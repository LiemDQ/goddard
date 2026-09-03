"""Matplotlib helpers for visualizing Method of Characteristics results.

Every function takes an optional ``ax`` and returns the axes it drew on, so plots
compose into larger figures. Matplotlib is imported lazily inside each function,
keeping ``import goddard`` fast and usable in headless contexts that never plot.

The helpers are deliberately tolerant of partial results: a solve that failed
partway still carries a net worth looking at, and seeing where the march stopped
is usually the whole point of plotting it.
"""

import numpy as np

from goddard._core import CharacteristicFamily, CharacteristicNet


__all__ = [
    "mesh_node_mask",
    "plot_characteristic_net",
    "plot_field",
    "plot_profile",
    "plot_exit_plane",
    "plot_front_diagnostics",
]


def _pyplot():
    """Import pyplot on demand, with a clearer message than a bare ImportError."""
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:  # pragma: no cover - depends on the environment
        raise ImportError(
            "goddard.plotting requires matplotlib. Install it with "
            "`pip install matplotlib`."
        ) from exc
    return plt


def _as_net(result_or_net):
    """Accept either a MocResult or a bare CharacteristicNet."""
    if isinstance(result_or_net, CharacteristicNet):
        return result_or_net
    return result_or_net.net


def plot_characteristic_net(result_or_net, ax=None, *, families=None, wall=True,
                            axis=True, linewidth=0.5, plus_color="tab:blue",
                            minus_color="tab:red", wall_color="k"):
    """Draw the characteristic mesh: every C+ and C- line in the net.

    Args:
        result_or_net: A MocResult or a CharacteristicNet.
        ax: Axes to draw on. A new figure is created when omitted.
        families: Iterable of CharacteristicFamily to draw. Defaults to both.
        wall: Draw the wall contour.
        axis: Draw the centerline.
        linewidth: Line width for the characteristics.
        plus_color: Colour for the C+ family.
        minus_color: Colour for the C- family.
        wall_color: Colour for the wall contour.

    Returns:
        The matplotlib Axes.
    """
    plt = _pyplot()
    net = _as_net(result_or_net)
    if ax is None:
        _, ax = plt.subplots()

    if families is None:
        families = (CharacteristicFamily.PLUS, CharacteristicFamily.MINUS)
    families = set(families)

    # Snapshot once: the columnar properties build a fresh array on every access,
    # so re-reading them inside the chain loop would be quadratic.
    x = np.asarray(net.x)
    y = np.asarray(net.y)
    chains = net.chains
    metadata = net.chain_metadata

    colors = {
        CharacteristicFamily.PLUS: plus_color,
        CharacteristicFamily.MINUS: minus_color,
    }
    labelled = set()
    for chain, meta in zip(chains, metadata):
        if meta.family not in families or len(chain) < 2:
            continue
        color = colors.get(meta.family, "0.5")
        label = None
        if meta.family not in labelled:
            label = "C+" if meta.family == CharacteristicFamily.PLUS else "C-"
            labelled.add(meta.family)
        ax.plot(x[chain], y[chain], color=color, linewidth=linewidth, label=label)

    if wall and len(net.wall_x) > 0:
        ax.plot(net.wall_x, net.wall_y, color=wall_color, linewidth=1.5, label="wall")
    if axis:
        ax.axhline(0.0, color="0.6", linewidth=0.8, linestyle="--")

    ax.set_xlabel("x")
    ax.set_ylabel("r")
    ax.set_aspect("equal", adjustable="datalim")
    if labelled or wall:
        ax.legend(loc="upper left", fontsize="small")
    return ax


def mesh_node_mask(net):
    """Boolean mask selecting the net points that are genuine solution nodes.

    A net carries a few points that were never produced by a unit process: the
    seeded throat-lip wall point, for one, which only bootstraps the wall march
    and holds placeholder zeros for Mach, pressure and temperature. Such points
    belong to no characteristic chain, which is what this tests for. Contouring
    over them puts a spurious cold spot at the throat and drags the colour scale
    down to zero.

    Args:
        net: A CharacteristicNet.

    Returns:
        A boolean numpy array with one entry per point in ``net.points``.
    """
    membership = net.membership
    return np.fromiter(
        (m.c_plus_chain_idx is not None or m.c_minus_chain_idx is not None
         for m in membership),
        dtype=bool,
        count=len(membership),
    )


def plot_field(result_or_net, ax=None, *, field="mach", levels=30, cmap="viridis",
               colorbar=True, wall=True, mesh_nodes_only=True):
    """Contour a flow quantity over the net's scattered points.

    Args:
        result_or_net: A MocResult or a CharacteristicNet.
        ax: Axes to draw on. A new figure is created when omitted.
        field: Name of a columnar net property ("mach", "theta", "pressure",
            "temperature", "nu", "mu", "gamma_s", "V").
        levels: Number of contour levels.
        cmap: Matplotlib colormap name.
        colorbar: Attach a colorbar to the axes.
        wall: Overlay the wall contour.
        mesh_nodes_only: Contour only genuine solution nodes, dropping bootstrap
            points that carry placeholder values. See :func:`mesh_node_mask`.
            Set False to contour every point in the net as-is.

    Returns:
        The matplotlib Axes.
    """
    plt = _pyplot()
    net = _as_net(result_or_net)
    if ax is None:
        _, ax = plt.subplots()

    values = np.asarray(getattr(net, field))
    x = np.asarray(net.x)
    y = np.asarray(net.y)
    if mesh_nodes_only:
        keep = mesh_node_mask(net)
        values, x, y = values[keep], x[keep], y[keep]
    if values.size < 3:
        raise ValueError(
            f"need at least 3 points to contour a field; got {values.size}"
        )

    contours = ax.tricontourf(x, y, values, levels=levels, cmap=cmap)
    if colorbar:
        ax.figure.colorbar(contours, ax=ax, label=field)
    if wall and len(net.wall_x) > 0:
        ax.plot(net.wall_x, net.wall_y, color="k", linewidth=1.5)

    ax.set_xlabel("x")
    ax.set_ylabel("r")
    ax.set_aspect("equal", adjustable="datalim")
    return ax


def plot_profile(profile, ax=None, *, mirror=False, label=None, **kwargs):
    """Draw a nozzle wall contour.

    Args:
        profile: A NozzleProfile, or a MocResult whose profile is used.
        ax: Axes to draw on. A new figure is created when omitted.
        mirror: Also draw the contour reflected about the centerline.
        label: Legend label for the contour.
        **kwargs: Forwarded to ``Axes.plot``.

    Returns:
        The matplotlib Axes.
    """
    plt = _pyplot()
    if hasattr(profile, "profile"):
        profile = profile.profile
    if ax is None:
        _, ax = plt.subplots()

    x = np.asarray(profile.x)
    y = np.asarray(profile.y)
    ax.plot(x, y, label=label, **kwargs)
    if mirror:
        ax.plot(x, -y, **kwargs)

    ax.set_xlabel("x")
    ax.set_ylabel("r")
    ax.set_aspect("equal", adjustable="datalim")
    if label is not None:
        ax.legend(loc="upper left", fontsize="small")
    return ax


def plot_exit_plane(result, ax=None, *, field="mach", **kwargs):
    """Plot a flow quantity across the exit plane, against radial station.

    Args:
        result: A MocResult.
        ax: Axes to draw on. A new figure is created when omitted.
        field: Name of an ExitPlane array ("mach", "theta", "pressure",
            "temperature", "gamma_s", "velocity").
        **kwargs: Forwarded to ``Axes.plot``.

    Returns:
        The matplotlib Axes.
    """
    plt = _pyplot()
    if ax is None:
        _, ax = plt.subplots()

    exit_plane = result.exit_plane
    values = np.asarray(getattr(exit_plane, field))
    ax.plot(values, np.asarray(exit_plane.y), **kwargs)
    ax.set_xlabel(field)
    ax.set_ylabel("r")
    return ax


def plot_front_diagnostics(result, axes=None):
    """Plot the marching front's per-pass geometry history.

    Three stacked panels: point spacing against its target, the worst mesh-cell
    aspect ratio, and the smallest spacelike margin. These are the quantities that
    expose a sampling void or a degenerating cell, neither of which is visible in
    the finished net.

    Args:
        result: A MocResult.
        axes: A sequence of three Axes. A new figure is created when omitted.

    Returns:
        The array of three Axes.

    Raises:
        ValueError: If the result carries no per-pass diagnostics.
    """
    plt = _pyplot()
    from goddard.convenience import pass_diagnostics_table

    table = pass_diagnostics_table(result)
    if table["pass_index"].size == 0:
        raise ValueError("result has no pass_diagnostics to plot")

    if axes is None:
        _, axes = plt.subplots(3, 1, sharex=True, figsize=(7, 8))
    passes = table["pass_index"]

    axes[0].plot(passes, table["min_spacing"], label="min")
    axes[0].plot(passes, table["mean_spacing"], label="mean")
    axes[0].plot(passes, table["max_spacing"], label="max")
    axes[0].plot(passes, table["target_spacing"], "k--", label="target")
    axes[0].set_ylabel("front spacing")
    axes[0].legend(fontsize="small")

    axes[1].plot(passes, table["max_cell_aspect"])
    axes[1].set_ylabel("max cell aspect")

    axes[2].plot(passes, table["min_spacelike_margin"])
    axes[2].axhline(0.0, color="tab:red", linewidth=0.8, linestyle="--")
    axes[2].set_ylabel("min spacelike margin")
    axes[2].set_xlabel("kernel pass")

    return axes

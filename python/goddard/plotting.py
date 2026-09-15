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
    "fan_apex_index",
    "mesh_node_mask",
    "trace_characteristics",
    "plot_characteristic_net",
    "plot_traced_characteristics",
    "plot_fronts",
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


def fan_apex_index(net):
    """Index of the throat lip a minimum-length design's expansion fan is centered on.

    Minimum-length design seeds every characteristic of the fan at one point -- the sharp
    throat lip -- and records each C- chain only from where it crosses the start line, so
    the segment of each ray between the lip and the start line belongs to the solution but
    is in no chain. That segment is what :func:`plot_characteristic_net` restores, and this
    is the point it draws back to.

    The lip is the net's first wall point. It is the bootstrap anchor rather than a solved
    node (see :func:`mesh_node_mask`), which is what distinguishes it here: a net whose
    first wall point carries a real Mach number was not built this way and has no fan apex.

    Args:
        net: A CharacteristicNet.

    Returns:
        The point index of the fan apex, or None if the net has no such point.
    """
    indices = np.asarray(net.wall_point_indices)
    if len(indices) == 0 or len(net.chains) == 0:
        return None
    apex = int(indices[0])
    return apex if net.points[apex].mach <= 0.0 else None


def _fan_ray_origins(net, apex):
    """Point indices where the fan's rays first meet the start line, one per ray.

    The start line is the first C+ chain the kernel creates. Every fan ray but the leading
    one goes on to spawn a C- chain from the point where it crosses that line, so those
    crossings are the C- chain origins lying on it; C- chains born later start at a wall
    reflection instead and are excluded. The leading ray spawns no C- chain -- it runs from
    the lip to the axis, where it reflects into the start line's own C+ -- so its crossing
    is the start line's axis end, taken here as the first point of that chain.

    Each origin is collinear with the apex along its own characteristic direction to
    machine precision, so a segment drawn to it is the ray itself, not an interpolation.
    """
    chain = [int(i) for i in net.chains[0]]
    start_line = set(chain)
    origins = {
        int(meta.origin_point_idx)
        for meta in net.chain_metadata
        if meta.family == CharacteristicFamily.MINUS
        and int(meta.origin_point_idx) in start_line
    }
    if chain:
        origins.add(chain[0])
    origins.discard(apex)
    return sorted(origins)


def plot_characteristic_net(result_or_net, ax=None, *, families=None, wall=True,
                            axis=True, linewidth=0.5, plus_color="tab:blue",
                            minus_color="tab:red", wall_color="k", fan_rays=True):
    """Draw the characteristic mesh: every C+ and C- line in the net.

    Analysis and Rao-design solves (the reference-plane march) prescribe a
    sequence of marching fronts instead of chain-paired characteristics, so
    CharacteristicNet.chains is empty for them; this draws CharacteristicNet.fronts
    (one polyline each) in that case instead of drawing nothing. Use
    :func:`plot_fronts` directly for more control over that drawing (e.g. thinning
    a fine march with ``every``). Fronts are near-vertical reference planes and show no
    wave structure; :func:`plot_traced_characteristics` reconstructs the characteristics
    themselves from such a net.

    Args:
        result_or_net: A MocResult or a CharacteristicNet.
        ax: Axes to draw on. A new figure is created when omitted.
        families: Iterable of CharacteristicFamily to draw. Defaults to both.
            Ignored when the net has no chains (fronts carry no family).
        wall: Draw the wall contour.
        axis: Draw the centerline.
        linewidth: Line width for the characteristics (or fronts).
        plus_color: Colour for the C+ family (also used for the fronts, when the
            net has no chains).
        minus_color: Colour for the C- family.
        wall_color: Colour for the wall contour.
        fan_rays: Extend a minimum-length design's fan characteristics back to the
            throat lip they are centered on. Without this the fan appears to emanate
            from the start line -- the locus where its rays first cross the leading
            C+ characteristic -- rather than from the lip. See :func:`fan_apex_index`.
            Ignored for nets that carry no such fan.

    Returns:
        The matplotlib Axes.
    """
    plt = _pyplot()
    net = _as_net(result_or_net)
    if ax is None:
        _, ax = plt.subplots()

    # Snapshot once: the columnar properties build a fresh array on every access,
    # so re-reading them inside the chain/front loop would be quadratic.
    x = np.asarray(net.x)
    y = np.asarray(net.y)
    chains = net.chains

    labelled = set()
    if len(chains) == 0 and len(net.fronts) > 0:
        for i, front in enumerate(net.fronts):
            idx = np.asarray(front)
            ax.plot(x[idx], y[idx], color=plus_color, linewidth=linewidth,
                    label="front" if i == 0 else None)
        labelled.add("front")
    else:
        if families is None:
            families = (CharacteristicFamily.PLUS, CharacteristicFamily.MINUS)
        families = set(families)
        metadata = net.chain_metadata

        colors = {
            CharacteristicFamily.PLUS: plus_color,
            CharacteristicFamily.MINUS: minus_color,
        }

        for chain, meta in zip(chains, metadata):
            if meta.family not in families or len(chain) < 2:
                continue
            color = colors.get(meta.family, "0.5")
            label = None
            if meta.family not in labelled:
                label = "C+" if meta.family == CharacteristicFamily.PLUS else "C-"
                labelled.add(meta.family)
            ax.plot(x[chain], y[chain], color=color, linewidth=linewidth, label=label)

        # The fan's rays run from the throat lip, but the net records each only from the
        # start line onward, so without this the fan appears to emanate from that line.
        # These segments are the missing upstream piece of rays already drawn above (and,
        # for the leading ray, the whole of it) -- each origin is collinear with the apex
        # along the ray, so this draws the characteristic rather than inventing geometry.
        apex = fan_apex_index(net) if fan_rays else None
        if apex is not None and CharacteristicFamily.MINUS in families:
            for origin in _fan_ray_origins(net, apex):
                label = None
                if CharacteristicFamily.MINUS not in labelled:
                    label = "C-"
                    labelled.add(CharacteristicFamily.MINUS)
                ax.plot(x[[apex, origin]], y[[apex, origin]], color=minus_color,
                        linewidth=linewidth, label=label)

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


def _ray_front_intersection(px, py, angle, fx, fy):
    """First forward intersection of a ray with a front polyline.

    Returns (x, y, segment, u) with u the position within that segment, or None when the
    ray misses -- which is how a traced characteristic learns it has left the domain
    through the wall or the axis rather than reaching the next front.
    """
    dx, dy = np.cos(angle), np.sin(angle)
    ax_, ay_ = fx[:-1], fy[:-1]
    ex, ey = fx[1:] - ax_, fy[1:] - ay_
    det = ex * dy - dx * ey
    usable = np.abs(det) > 1e-14
    wx, wy = ax_ - px, ay_ - py
    s = np.full(ax_.shape, -1.0)
    u = np.full(ax_.shape, -1.0)
    s[usable] = (-wx[usable] * ey[usable] + ex[usable] * wy[usable]) / det[usable]
    u[usable] = (dx * wy[usable] - dy * wx[usable]) / det[usable]
    hit = usable & (s > 1e-12) & (u >= -1e-9) & (u <= 1.0 + 1e-9)
    if not hit.any():
        return None
    # Nearest forward crossing: a characteristic that grazes a folded front must take the
    # first one, not whichever segment happens to come first in the array.
    i = int(np.argmin(np.where(hit, s, np.inf)))
    ui = min(max(u[i], 0.0), 1.0)
    return ax_[i] + ui * ex[i], ay_[i] + ui * ey[i], i, ui


def trace_characteristics(result_or_net, *, families=None, every=1, seed_start_front=True):
    """Reconstruct characteristic lines from an inverse-march net.

    Analysis and Rao-design solves march reference planes rather than pairing
    characteristics, so their mesh topology is the sequence of fronts in
    CharacteristicNet.fronts and CharacteristicNet.chains is empty. The fronts are very
    nearly vertical, which makes the wave structure invisible in them. This integrates the
    characteristic directions through the solved field to recover it: from a seed point it
    steps front to front along dy/dx = tan(theta +/- mu), taking theta and mu by linear
    interpolation along each front it crosses, with one corrector pass per step (which is
    enough -- the angle iteration converges immediately; the accuracy limit is the
    interpolation along the front).

    These curves are reconstructed after the fact, not the mesh the solver used, so treat
    them as a visual check on the field rather than as solution data. The reconstruction is
    convergent but only first order in the front count: in planar flow, where theta -/+ nu
    is exactly conserved along a C+/C- line, the drift along a traced line runs about 1.8,
    0.9, 0.43 and 0.18 degrees for 45, 97, 205 and 430 fronts.

    Args:
        result_or_net: A MocResult or a CharacteristicNet, from an inverse-march solve.
        families: Iterable of CharacteristicFamily to trace. Defaults to both.
        every: Seed from every Nth front. Raise it to thin a fine march.
        seed_start_front: Also seed every point of the first front, which fills in the
            throat region that the axis and wall seeds reach only further downstream.

    Returns:
        List of ``(family, points)`` pairs, each ``points`` an ``(n, 2)`` array of x, y.

    Raises:
        ValueError: If the net carries no fronts.
    """
    net = _as_net(result_or_net)
    fronts = net.fronts
    if len(fronts) == 0:
        raise ValueError(
            "net carries no fronts to trace through; a minimum-length design net already "
            "holds its characteristics in net.chains, so plot them with "
            "plot_characteristic_net instead"
        )
    if families is None:
        families = (CharacteristicFamily.PLUS, CharacteristicFamily.MINUS)
    families = set(families)

    # Snapshot the columnar properties once, then slice per front: each access rebuilds a
    # whole-net array, so reading them inside the trace loop would be quadratic.
    def by_front(field):
        values = np.asarray(getattr(net, field))
        return [values[np.asarray(f)] for f in fronts]

    X, Y = by_front("x"), by_front("y")
    TH, MU = by_front("theta"), by_front("mu")

    def trace(sign, k0, j):
        px, py = X[k0][j], Y[k0][j]
        theta, mu = TH[k0][j], MU[k0][j]
        points = [(px, py)]
        for k in range(k0 + 1, len(fronts)):
            angle = theta + sign * mu
            hit = _ray_front_intersection(px, py, angle, X[k], Y[k])
            if hit is None:
                break
            _, _, i, u = hit
            theta_q = TH[k][i] * (1 - u) + TH[k][i + 1] * u
            mu_q = MU[k][i] * (1 - u) + MU[k][i + 1] * u
            hit = _ray_front_intersection(
                px, py, 0.5 * (angle + theta_q + sign * mu_q), X[k], Y[k])
            if hit is None:
                break
            px, py, i, u = hit
            theta = TH[k][i] * (1 - u) + TH[k][i + 1] * u
            mu = MU[k][i] * (1 - u) + MU[k][i + 1] * u
            points.append((px, py))
        return np.asarray(points)

    # A C+ climbs away from the axis and a C- descends from the wall, so each family is
    # seeded on the boundary it leaves, and both on the first front to fill the throat.
    seeds = []
    for k in range(0, len(fronts) - 1, max(1, every)):
        seeds.append((CharacteristicFamily.PLUS, k, 0))
        seeds.append((CharacteristicFamily.MINUS, k, len(fronts[k]) - 1))
    if seed_start_front:
        for j in range(len(fronts[0])):
            seeds.append((CharacteristicFamily.PLUS, 0, j))
            seeds.append((CharacteristicFamily.MINUS, 0, j))

    traced = []
    for family, k, j in seeds:
        if family not in families:
            continue
        points = trace(+1 if family == CharacteristicFamily.PLUS else -1, k, j)
        if len(points) > 1:
            traced.append((family, points))
    return traced


def plot_traced_characteristics(result_or_net, ax=None, *, families=None, every=1,
                                seed_start_front=True, linewidth=0.5,
                                plus_color="tab:blue", minus_color="tab:red",
                                wall=True, axis=True, wall_color="k"):
    """Draw characteristic lines reconstructed from an inverse-march net.

    The counterpart of :func:`plot_characteristic_net` for analysis and Rao-design solves,
    whose nets hold fronts rather than chains. See :func:`trace_characteristics` for how
    the curves are recovered and how far to trust them.

    Args:
        result_or_net: A MocResult or a CharacteristicNet, from an inverse-march solve.
        ax: Axes to draw on. A new figure is created when omitted.
        families: Iterable of CharacteristicFamily to draw. Defaults to both.
        every: Seed from every Nth front.
        seed_start_front: Also seed every point of the first front.
        linewidth: Line width for the characteristics.
        plus_color: Colour for the C+ family.
        minus_color: Colour for the C- family.
        wall: Draw the wall contour.
        axis: Draw the centerline.
        wall_color: Colour for the wall contour.

    Returns:
        The matplotlib Axes.
    """
    plt = _pyplot()
    net = _as_net(result_or_net)
    if ax is None:
        _, ax = plt.subplots()

    colors = {
        CharacteristicFamily.PLUS: plus_color,
        CharacteristicFamily.MINUS: minus_color,
    }
    labels = {CharacteristicFamily.PLUS: "C+", CharacteristicFamily.MINUS: "C-"}
    labelled = set()
    for family, points in trace_characteristics(
            net, families=families, every=every, seed_start_front=seed_start_front):
        label = None
        if family not in labelled:
            label = labels[family]
            labelled.add(family)
        ax.plot(points[:, 0], points[:, 1], color=colors.get(family, "0.5"),
                linewidth=linewidth, label=label)

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


def plot_fronts(result_or_net, ax=None, *, every=1, **kwargs):
    """Draw every marching front of an inverse-march net as one polyline each.

    Meaningful for analysis and Rao-design solves (the reference-plane march),
    whose mesh topology is the sequence of fronts in CharacteristicNet.fronts
    (axis to wall, in march order) rather than the chain-paired characteristics
    that minimum-length design (the chain ladder) builds. A minimum-length design
    net carries no fronts, so this draws nothing for one.

    Args:
        result_or_net: A MocResult or a CharacteristicNet.
        ax: Axes to draw on. A new figure is created when omitted.
        every: Draw every Nth front only, to declutter a fine march. 1 draws all.
        **kwargs: Forwarded to ``Axes.plot`` (e.g. color, linewidth).

    Returns:
        The matplotlib Axes.
    """
    plt = _pyplot()
    net = _as_net(result_or_net)
    if ax is None:
        _, ax = plt.subplots()

    x = np.asarray(net.x)
    y = np.asarray(net.y)
    kwargs.setdefault("color", "tab:blue")
    kwargs.setdefault("linewidth", 0.5)
    for front in net.fronts[::every]:
        idx = np.asarray(front)
        ax.plot(x[idx], y[idx], **kwargs)

    ax.set_xlabel("x")
    ax.set_ylabel("r")
    ax.set_aspect("equal", adjustable="datalim")
    return ax


def mesh_node_mask(net):
    """Boolean mask selecting the net points that are genuine solution nodes.

    A minimum-length design net (the chain ladder) carries a few points that were
    never produced by a unit process: the seeded throat-lip wall point, for one,
    which only bootstraps the wall march and holds placeholder zeros for Mach,
    pressure and temperature. Such points belong to no characteristic chain, which
    is what this tests for. Contouring over them puts a spurious cold spot at the
    throat and drags the colour scale down to zero.

    Analysis and Rao-design solves (the reference-plane march) build no chains at
    all (CharacteristicNet.chains is empty; its topology lives in
    CharacteristicNet.fronts instead), so every point has empty membership
    regardless of whether it is a genuine solution node -- for those solves every
    point *is* one, so this treats them all as mesh nodes rather than reading
    "no chain membership" as "bootstrap placeholder".

    Args:
        net: A CharacteristicNet.

    Returns:
        A boolean numpy array with one entry per point in ``net.points``.
    """
    membership = net.membership
    has_any_chain_entry = any(
        m.c_plus_chain_idx is not None or m.c_minus_chain_idx is not None
        for m in membership
    )
    if not has_any_chain_entry:
        return np.ones(len(membership), dtype=bool)
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

    Three stacked panels: point spacing (min/mean/max), the axial step size, and
    the front's axis and wall stations. These are the quantities that expose a
    sampling void, a stalled march, or a degenerating front, none of which is
    visible in the finished net.

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
    axes[0].set_ylabel("front spacing")
    axes[0].legend(fontsize="small")

    axes[1].plot(passes, table["step_dx"])
    axes[1].set_ylabel("step dx")

    axes[2].plot(passes, table["front_axis_x"], label="axis")
    axes[2].plot(passes, table["front_wall_x"], label="wall")
    axes[2].set_ylabel("front station x")
    axes[2].set_xlabel("kernel pass")
    axes[2].legend(fontsize="small")

    return axes

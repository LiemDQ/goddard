"""Rocket solver front end: array views over a list of rocket results."""

from __future__ import annotations

from collections.abc import Iterable, Iterator, Sequence
from typing import Any, overload

import numpy as np

from goddard import _core
from goddard._core import RocketPoint, RocketResult, SolveStatus

_STATIONS = ("chamber", "stagnation", "combustion_end", "throat")
_PERFORMANCE_FIELDS = ("pressure_ratio", "area_ratio", "mach_number", "cstar", "CF", "isp",
                       "ivac", "isp_s", "ivac_s")


class RocketResults(Sequence[RocketResult]):
    """Results of a list of rocket points, in the order of the points.

    Indexing and iteration give the `RocketResult` of each point. The methods gather one
    quantity over every point into a NumPy array, so a sweep can be plotted without a loop.
    Exit quantities have one column per exit; points with fewer exits than the longest are
    padded with NaN, as are failed points.
    """

    def __init__(self, results: Iterable[RocketResult]):
        self._results = list(results)

    @overload
    def __getitem__(self, index: int) -> RocketResult: ...
    @overload
    def __getitem__(self, index: slice) -> RocketResults: ...

    def __getitem__(self, index):
        if isinstance(index, slice):
            return RocketResults(self._results[index])
        return self._results[index]

    def __len__(self) -> int:
        return len(self._results)

    def __iter__(self) -> Iterator[RocketResult]:
        return iter(self._results)

    def __repr__(self) -> str:
        failed = len(self.failed())
        return f"<RocketResults {len(self)} point(s), {failed} not fully solved>"

    @property
    def status(self) -> SolveStatus | None:
        """Status of the first point that was not fully solved, or an ok status.

        None for an empty list of results.
        """
        for result in self._results:
            if not result.status.ok:
                return result.status
        return self._results[0].status if self._results else None

    def failed(self) -> list[int]:
        """Indices of the points with a station that was not solved."""
        return [i for i, result in enumerate(self._results) if not result.status.ok]

    def _num_exits(self) -> int:
        return max((len(result.exits) for result in self._results), default=0)

    def performance(self, field: str) -> np.ndarray:
        """One performance quantity at every exit of every point.

        Args:
            field: A `RocketPerformance` attribute, e.g. "ivac", "isp_s" or "cstar".

        Returns:
            Array of shape (number of points, number of exits).
        """
        if field not in _PERFORMANCE_FIELDS:
            raise ValueError(f"unknown performance field {field!r}; "
                             f"expected one of {', '.join(_PERFORMANCE_FIELDS)}")
        values = np.full((len(self), self._num_exits()), np.nan)
        for i, result in enumerate(self._results):
            for j, performance in enumerate(result.performance):
                values[i, j] = getattr(performance, field)
        return values

    def station(self, station: str, field: str) -> np.ndarray:
        """One state quantity at a station of every point.

        Args:
            station: "chamber", "stagnation", "combustion_end", "throat" or "exits".
            field: A `ThermodynamicState` attribute, e.g. "temperature", or "velocity",
                "mach" or "area_ratio".

        Returns:
            Array of shape (number of points,), or (number of points, number of exits) for
            "exits".
        """
        def value(station_result: Any) -> float:
            if field in ("velocity", "mach", "area_ratio"):
                return getattr(station_result, field)
            return getattr(station_result.thermo, field)

        if station == "exits":
            values = np.full((len(self), self._num_exits()), np.nan)
            for i, result in enumerate(self._results):
                for j, exit_station in enumerate(result.exits):
                    values[i, j] = value(exit_station)
            return values
        if station not in _STATIONS:
            raise ValueError(f"unknown station {station!r}; expected one of "
                             f"{', '.join(_STATIONS)} or 'exits'")
        return np.array([value(getattr(result, station)) for result in self._results])

    def point(self, field: str) -> np.ndarray:
        """One input field of every point, e.g. "mixture_ratio" or "chamber_pressure"."""
        return np.array([getattr(result.point, field) for result in self._results])

    def report(self) -> str:
        """CEA-style text report, one page per point."""
        return _core.rocket_report(self._results)

    def to_dataframe(self) -> Any:
        """One row per exit of every point, with the point inputs and the exit performance.

        Requires pandas.
        """
        import pandas as pd

        rows = []
        for index, result in enumerate(self._results):
            point = result.point
            for exit_index, (exit_station, performance) in enumerate(
                    zip(result.exits, result.performance)):
                row = {
                    "point": index,
                    "label": point.label,
                    "mixture_ratio": point.mixture_ratio,
                    "chamber_pressure": point.chamber_pressure,
                    "exit": exit_index,
                    "expansion_ratio": point.expansion_ratios[exit_index],
                    "chamber_temperature": result.chamber.thermo.temperature,
                    "exit_temperature": exit_station.thermo.temperature,
                    "exit_pressure": exit_station.thermo.pressure,
                    "ok": performance.status.ok,
                }
                for field in _PERFORMANCE_FIELDS:
                    row[field] = getattr(performance, field)
                rows.append(row)
        return pd.DataFrame(rows)


class RocketSolver(_core.RocketSolver):
    __doc__ = _core.RocketSolver.__doc__

    @overload
    def solve(self, point: RocketPoint) -> RocketResult: ...
    @overload
    def solve(self, points: Iterable[RocketPoint]) -> RocketResults: ...

    def solve(self, points: RocketPoint | Iterable[RocketPoint]
              ) -> RocketResult | RocketResults:
        """Solve one point, or every point of an iterable.

        Args:
            points: A `RocketPoint`, or any iterable of them (a list, a generator).

        Returns:
            The `RocketResult` of a single point, or a `RocketResults` with one result per
            point, in order. Every point is checked before any is solved.
        """
        if isinstance(points, RocketPoint):
            return super().solve(points)
        return RocketResults(super().solve(list(points)))

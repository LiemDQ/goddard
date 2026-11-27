"""Method of Characteristics solver diagnostics.

Per-pass marching-front statistics, start-line diagnostics and the bookkeeping of the
characteristic net's chains. These are for studying the solver itself (grid convergence,
step limiting, front quality); a design or analysis needs only `MocResult`.
"""

from goddard._core import (
    ChainMetadata,
    MocFrontShear,
    MocInitDiagnostics,
    MocPassDiagnostics,
    MocStepLimiter,
    PointMembership,
    summarize_front_shear,
)
from goddard.convenience import pass_diagnostics_table

__all__ = [
    "ChainMetadata",
    "MocFrontShear",
    "MocInitDiagnostics",
    "MocPassDiagnostics",
    "MocStepLimiter",
    "PointMembership",
    "pass_diagnostics_table",
    "summarize_front_shear",
]

"""Public API for BioEqPy."""

from bioeqpy.analysis import analyze
from bioeqpy.core.types import ANOVAResult, BEStudy, BEResult, CIResult, DesignSpec

__all__ = [
    "ANOVAResult",
    "BEStudy",
    "BEResult",
    "CIResult",
    "DesignSpec",
    "analyze",
]

__version__ = "0.1.0"


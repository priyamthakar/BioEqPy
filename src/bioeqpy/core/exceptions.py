"""Custom exceptions raised by BioEqPy."""


class BioEqError(Exception):
    """Base exception for BioEqPy errors."""


class DesignError(BioEqError):
    """Raised when a study design cannot be detected or validated."""


class ValidationError(BioEqError):
    """Raised when input data fail schema or numeric validation."""


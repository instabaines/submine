"""Submine exception hierarchy.

These exceptions provide stable, semantically meaningful error types that
downstream applications can catch without parsing strings.
"""

from __future__ import annotations


class SubmineError(Exception):
    """Base class for all library-defined exceptions."""


class SubmineInputError(SubmineError, ValueError):
    """Raised when user-supplied inputs (files, graphs, parameters) are invalid."""


class ParameterValidationError(SubmineInputError):
    """Raised when algorithm parameters fail validation."""


class BackendUnavailableError(SubmineError, RuntimeError):
    """Raised when an optional external backend (binary, JVM, etc.) is unavailable."""


class BackendExecutionError(SubmineError, RuntimeError):
    """Raised when an external backend fails during execution."""


class ResourceLimitError(SubmineInputError):
    """Raised when an input exceeds configured resource limits."""

__all__ = [
    "SubmineError",
    "SubmineInputError",
    "ParameterValidationError",
    "BackendUnavailableError",
    "BackendExecutionError",
    "ResourceLimitError",
]


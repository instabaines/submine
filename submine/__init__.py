from __future__ import annotations

__version__ = "0.1.2"

from .registry import available_algorithms
from .errors import (
    SubmineError,
    SubmineInputError,
    ParameterValidationError,
    BackendUnavailableError,
    BackendExecutionError,
    ResourceLimitError,
)


# Import algorithms so they register themselves via @register
# (you can add more as you implement them)
from .algorithms import gspan, sopagrami




def get_mining_algorithm(name: str):
    key = name.lower()
    try:
        return available_algorithms[key]
    except KeyError:
        raise ValueError(
            f"Unknown algorithm '{name}'. "
            f"Available: {sorted(available_algorithms.keys())}"
        )

def list_algorithms():
    return sorted(available_algorithms.keys())

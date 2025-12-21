"""Environment checks for submine."""

from __future__ import annotations

import shutil
from typing import Optional


def is_tool_available(name: str) -> bool:
    """Return True if a given executable exists on the system PATH."""
    return shutil.which(name) is not None
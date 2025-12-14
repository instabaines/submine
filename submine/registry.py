# submine/registry.py
from __future__ import annotations

from typing import Dict, Type

# We intentionally don't import SubgraphMiner here to avoid cycles.
# We just store "type" and let base.py handle typing.
available_algorithms: Dict[str, Type] = {}

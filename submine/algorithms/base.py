# submine/algorithms/base.py
from __future__ import annotations

import subprocess
from abc import ABC, abstractmethod
from pathlib import Path
from typing import List, Optional

from ..registry import available_algorithms   # <-- IMPORTANT
from ..core.graph import Graph
from ..core.result import MiningResult
from ..utils.logging import get_logger
from typing import Iterable
__all__ = ["SubgraphMiner", "register"]


class SubgraphMiner(ABC):
    name: str = "base"

    def __init__(self, verbose: bool = False) -> None:
        self.verbose = verbose
        self.logger = get_logger(self.__class__.__name__)
        if self.verbose:
            self.logger.setLevel("DEBUG")

    @abstractmethod
    def mine(self, graphs: Iterable[Graph], min_support: int, **kwargs) -> MiningResult:
        raise NotImplementedError

    def check_availability(self) -> None:
        return None

    def run_external(self, cmd: List[str], *, cwd: Optional[Path] = None) -> str:
        self.logger.debug("Running external command: %s", " ".join(cmd))
        completed = subprocess.run(cmd, cwd=cwd, text=True, capture_output=True)
        self.logger.debug("Command stdout: %s", completed.stdout)
        if completed.returncode != 0:
            self.logger.error("Command failed with stderr: %s", completed.stderr)
            raise RuntimeError(
                f"Command '{' '.join(cmd)}' failed with exit code {completed.returncode}\n"
                f"stderr:\n{completed.stderr}"
            )
        return completed.stdout


def register(cls: type[SubgraphMiner]) -> type[SubgraphMiner]:
    if not issubclass(cls, SubgraphMiner):
        raise TypeError("Only subclasses of SubgraphMiner can be registered")

    name = getattr(cls, "name", None)
    if not isinstance(name, str):
        raise TypeError("Subgraph miner must define a string 'name' attribute")

    key = name.lower()
    if key in available_algorithms:
        raise ValueError(f"Algorithm '{name}' is already registered")

    available_algorithms[key] = cls
    return cls

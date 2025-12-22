# submine/algorithms/base.py
from __future__ import annotations

import subprocess
from abc import ABC, abstractmethod
from pathlib import Path
from typing import List, Optional


from ..registry import available_algorithms   
from ..core.graph import Graph
from ..core.result import MiningResult
from ..utils.logging import get_logger
from ..errors import BackendExecutionError, ParameterValidationError
from typing import Iterable
__all__ = ["SubgraphMiner", "register"]


class SubgraphMiner(ABC):
    name: str = "base"

    # Native file input contract
    # -------------------------
    # Some algorithms consume an on-disk dataset in a specific format.
    # If `expected_input_format` is set (e.g., "lg" or "gspan"), the high-level
    # API can transcode user-provided files to this format and call `mine_native`.
    #
    # By default, miners operate on in-memory Graph objects via `mine()`.
    expected_input_format: str | None = None
    multi_graph_policy: str = "reject"  # reject | batch | merge (reserved)

    # Weight handling
    # -------------
    # Most classical subgraph miners operate on *labeled* (unweighted) graphs.
    # We therefore make weight support explicit. If the input graph contains
    # weights and the algorithm does not support them, the weight_strategy
    # controls what happens.
    supports_weighted: bool = False
    weight_strategy: str = "ignore"  # one of: ignore | reject

    def __init__(self, verbose: bool = False) -> None:
        self.verbose = verbose
        self.logger = get_logger(self.__class__.__name__)
        if self.verbose:
            self.logger.setLevel("DEBUG")

    def _handle_weights(self, graphs: Iterable[Graph]) -> Iterable[Graph]:
        """Apply the configured weight strategy to input graphs.

        - If the algorithm supports weights: pass through.
        - If it does not and the graph is weighted:
            * ignore: drop weights (treat as unweighted)
            * reject: raise
        """
        for g in graphs:
            if getattr(g, "is_weighted", False) and g.is_weighted and not self.supports_weighted:
                if self.weight_strategy == "reject":
                    raise ValueError(
                        f"Algorithm '{self.name}' does not support weighted graphs; "
                        "set weight_strategy='ignore' to drop weights explicitly."
                    )
                # ignore: drop weights
                g.edge_weights = None
            yield g

    @abstractmethod
    def mine(self, graphs: Iterable[Graph], min_support: int, **kwargs) -> MiningResult:
        raise NotImplementedError

    def mine_native(self, path: str | Path, min_support: int, **kwargs) -> MiningResult:
        """Run the miner on a native on-disk dataset.

        Miners with `expected_input_format != None` should override this method.
        The default implementation indicates that the miner does not accept a
        native path entrypoint.
        """
        raise NotImplementedError(
            f"Algorithm '{self.name}' does not implement mine_native(); "
            "use mine(graphs=...) instead."
        )

    def check_availability(self) -> None:
        return None

    def run_external(
        self,
        cmd: List[str],
        *,
        cwd: Optional[Path] = None,
        timeout_s: int = 300,
        env: Optional[dict[str, str]] = None,
    ) -> str:
        """Run an external command defensively.

        - Uses ``shell=False`` implicitly (we pass a list).
        - Applies a default timeout to avoid hung processes.
        - Captures stdout/stderr for error reporting.
        """
        if not cmd or not isinstance(cmd[0], str):
            raise ParameterValidationError("cmd must be a non-empty list of strings")

        # Basic hardening against accidental injection via newlines/NULs.
        for part in cmd:
            if not isinstance(part, str):
                raise TypeError("All cmd parts must be strings")
            if "\x00" in part or "\n" in part or "\r" in part:
                raise ParameterValidationError("Unsafe characters in command argument")

        self.logger.debug("Running external command: %s", " ".join(cmd))
        completed = subprocess.run(
            cmd,
            cwd=cwd,
            text=True,
            capture_output=True,
            timeout=timeout_s,
            env=env,
            check=False,
            close_fds=True,
        )
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

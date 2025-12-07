"""Entry point for the submine command line interface.

Use this CLI to run frequent subgraph mining algorithms from the shell.
It supports selecting an algorithm, loading a dataset and specifying
common parameters such as the minimum support threshold. The results
are printed to standard output.

Example::

    python -m submine.cli.main --algorithm gspan --dataset toy --min-support 2

"""

from __future__ import annotations

import argparse
from typing import List

from .. import get_algorithm, load_dataset


def parse_args(argv: List[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run frequent subgraph mining algorithms")
    parser.add_argument(
        "--algorithm",
        "-a",
        required=True,
        help="Name of the algorithm to run (e.g., gspan, grami)"
    )
    parser.add_argument(
        "--dataset",
        "-d",
        default="toy",
        help="Dataset name to load (e.g., toy, mutag, enzymes)"
    )
    parser.add_argument(
        "--min-support",
        "-s",
        type=int,
        default=1,
        help="Minimum support threshold (positive integer)"
    )
    parser.add_argument(
        "--top-k",
        "-k",
        type=int,
        default=5,
        help="Print the top K subgraphs by support"
    )
    parser.add_argument(
        "--verbose",
        "-v",
        action="store_true",
        help="Enable verbose logging"
    )
    return parser.parse_args(argv)


def main(argv: List[str] | None = None) -> None:
    args = parse_args(argv)
    # Load dataset
    try:
        graphs = load_dataset(args.dataset)
    except Exception as e:
        raise SystemExit(f"Failed to load dataset '{args.dataset}': {e}")
    # Instantiate algorithm
    try:
        miner = get_algorithm(args.algorithm, verbose=args.verbose)
    except KeyError as e:
        raise SystemExit(str(e))
    # Run mining
    try:
        result = miner.mine(graphs, min_support=args.min_support)
    except NotImplementedError as e:
        raise SystemExit(str(e))
    except Exception as e:
        raise SystemExit(f"Error while running algorithm '{args.algorithm}': {e}")
    # Print results
    top = result.top_k(args.top_k)
    print(f"Found {len(result)} frequent subgraphs (displaying top {len(top)})")
    for idx, fs in enumerate(top, start=1):
        # Provide a simple textual representation
        print(f"#{idx}: support={fs.support}, nodes={fs.pattern.number_of_nodes()}, edges={fs.pattern.number_of_edges()}")


if __name__ == "__main__":
    main()
# submine

*submine* is a modular Python library for frequent subgraph mining that provides a unified, safe, and extensible interface over heterogeneous mining algorithms implemented in Python, C/C++, and Java.

The library is designed to support research-grade reproducibility and production-safe execution, while remaining lightweight and backend-agnostic.

## Core Design Principles

* Algorithm-agnostic API
Users select an algorithm; submine handles format adaptation and execution.

* No redundant graph rewrites
Input graphs are converted directly into the format required by the selected algorithm.

* Strict input validation & safety
Resource limits, parameter validation, and hardened subprocess execution are enforced by default.

* Clean extensibility model
New algorithms can be plugged in without modifying core logic.

## Package Layout
```
submine/
├── api.py                  # Public API entrypoint
├── errors.py               # Structured exception hierarchy
├── core/
│   ├── graph.py            # Canonical graph representation
│   └── result.py           # Mining result containers
├── algorithms/
│   ├── base.py             # Base miner abstraction
│   ├── gspan.py            # gSpan wrapper / implementation
│   ├── sopagrami.py        # SoPaGraMi backend wrapper
│   └── ...
├── io/
│   ├── common.py           # Shared readers / writers
│   ├── transcode.py        # Format detection & conversion
│   ├── gspan.py            # gSpan I/O
│   ├── sopagrami.py        # .lg I/O
│   └── ...
├── utils/
│   └── checks.py           # Input validation & resource limits
└── tests/
    ├── unit/
    └── functional/
```

## Supported Algorithms

| Algorithm | Graph Type         | Backend      | Notes                    |
| --------- | ------------------ | ------------ | ------------------------ |
| gSpan     | Multiple graphs    | Python / C++ | Frequent subgraph mining |
| SoPaGraMi | Single large graph | C++          | Social pattern mining    |

Each algorithm declares:

* required input format,
* parameter schema,
* execution strategy (in-process vs subprocess).

## Supported Input Formats

submine accepts graphs in multiple formats and converts them directly into the format required by the selected algorithm:

* Edge list (.txt, .edgelist)

* gSpan datasets (.data, .data.x, .data.N)

* SoPaGraMi .lg

* GEXF (.gexf)

* Format detection is automatic and deterministic.

## Installation
### Runtime installation


```bash
pip install submine

```
### Dev installation

```bash
pip install -e ".[dev]"
```

## Basic Usage

```python
from submine.api import mine_subgraphs

results = mine_subgraphs(
    data="graph.data",
    algorithm="gspan",
    min_support=5
)

```

For SoPaGraMi:
```python
results = mine_subgraphs("citeseer.lg", algorithm="sopagrami", 
                    min_support=100,
                    sorted_seeds=4,
                    dump_images_csv=True, #if you want to dump the images
                    dump_sample_embeddings=True, #if you want to dump the embeddings, not implemented yet
                    out_dir= "." # /path/to/dir to save the images and embeddings default is ./sopagrami_result
                    )
```
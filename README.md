# submine

**submine** is a research‑grade Python library for frequent subgraph mining that provides a unified, safe, and extensible interface over heterogeneous mining algorithms implemented in Python, C++, and Java.

The goal of _submine_ is to let users focus on **what** to mine rather than **how** each algorithm expects its input. Users select an algorithm and parameters; _submine_ automatically validates inputs, converts graph formats, and executes the backend in a controlled and reproducible manner.

---

## Key Features

- **Algorithm‑centric API**
  You specify the mining algorithm and parameters; _submine_ handles format adaptation and execution.

- **Direct format transcoding (no redundant rewrites)**
  Input graphs are converted directly into the native format required by the selected algorithm.

- **Multi‑format graph support**
  Edge lists, gSpan datasets, single‑graph `.lg` files, and GEXF are supported out of the box.

- **Safe and reproducible execution**
  Parameter validation, deterministic format detection, and hardened subprocess execution are enforced by default.

- **Extensible design**
  New algorithms can be added via a clean backend interface without modifying core logic.

---

## Supported Algorithms

### gSpan (Frequent Subgraph Mining)

- **Graph type:** Multiple graphs (transactional dataset)
- **Typical use case:** Discovering frequent substructures across many graphs
- **Backend:** C++

The gSpan backend in _submine_ is a C++ implementation adapted and extended from the widely used **gBoost / gSpan reference implementations**, with additional input validation, format handling, and Python bindings for safe integration.

### SoPaGraMi (Single‑Graph Pattern Mining)

- **Graph type:** Single large graph
- **Typical use case:** Social, biological, or information networks
- **Backend:** C++

SoPaGraMi is used for scalable subgraph mining on a single graph, where frequency is defined structurally rather than transactionally.

---

## Supported Input Formats

_submine_ automatically detects the input format and converts it to the format required by the chosen algorithm:

- **Edge lists**: `.txt`, `.edgelist`
- **gSpan datasets**: `.data`, `.data.x`, `.data.N`
- **SoPaGraMi graphs**: `.lg`
- **GEXF**: `.gexf`

Format detection is deterministic and does not rely on user‑supplied flags.

---

## Installation

### Standard installation

```bash
pip install submine
```

### Development installation

```bash
pip install -e ".[dev]"
```

---

## Basic Usage

### gSpan example

```python
from submine.api import mine_subgraphs

results = mine_subgraphs(
    data="graphs.data",
    algorithm="gspan",
    min_support=5
)
```

**Parameters**

- `data` (str or path): Path to the input graph dataset
- `algorithm` (str): Mining algorithm (`"gspan"`, `"sopagrami"`, …)
- `min_support` (int): Minimum support threshold (algorithm‑specific semantics)

---

### SoPaGraMi example

```python
results = mine_subgraphs(
    data="citeseer.lg",
    algorithm="sopagrami",
    min_support=100,
    sorted_seeds=4,
    dump_images_csv=True,
    dump_sample_embeddings=False,
    out_dir="."
)
```

**SoPaGraMi‑specific parameters**

- `min_support` (int): Minimum frequency threshold
- `max_edges` (int): maximum pattern size
- `sorted_seeds` (int): Seed sorting strategy (implementation‑specific)
- `dump_images_csv` (bool): Whether to dump pattern images as CSV metadata
- `dump_sample_embeddings` (bool): Whether to dump sample embeddings (experimental)
- `out_dir` (str or path): Output directory for results (default: `./sopagrami_result`)

---

### Optimization (Specific to SoPaGrami): Maximum Edge Limit (`max_edges`)

#### The Problem: Combinatorial Explosion
Even on small graphs (e.g., < 4,000 edges), subgraph mining can hang indefinitely if the graph contains **dense clusters** or **hub nodes** (nodes with high degrees of identical labels).

In these dense regions, the number of valid subgraphs grows **exponentially** with the number of edges. 
* A pattern with 5 edges might have 100 embeddings.
* A pattern with 20 edges might have $10^9$ (billions) of embeddings.

Without a limit, the algorithm attempts to enumerate every single one of these "super-patterns," causing the system to stall (mining for hours/days without progress).

## The Solution: `max_edges`
We have introduced a `max_edges` parameter to the configuration. This acts as a **hard depth limit** for the search tree.

* **Mechanism:** During the candidate generation phase (`SUBGRAPHEXTENSION`), if a pattern's edge count reaches `max_edges`, the algorithm **stops expanding** that branch.
* **Default Value:** 10 (Conservative start).

## Implications & Trade-offs

| Implication | Description |
| :--- | :--- |
| **Performance** | **Drastic Speedup.** Prevents the algorithm from entering "infinite" search spaces in dense cliques. |
| **Completeness** | **Partial Loss.** You will not discover frequent patterns that strictly require > `max_edges` to be defined. You will only find their sub-components (fragments of size `max_edges`). |
| **Relevance** | **High.** In most domain applications (bioinformatics, social networks), distinct functional motifs rarely exceed 10-15 edges. Larger patterns are often just "hairball" noise. |

## Tuning Guide
* **Start Small:** Set `max_edges = 5`. The code should finish in seconds.
* **Scale Up:** Increment to 6, 7, 8... until the runtime becomes unacceptable.
* **Disable:** Set `max_edges = -1` (or a very large number) to disable the limit, but be warned of potential hangs on dense graphs.
## Design Philosophy

- **No algorithm‑specific I/O burden on the user**
  Users never manually convert graph formats.

- **Minimal assumptions about graph structure**
  Directed/undirected and labeled/unlabeled graphs are handled at the backend level.

- **Research‑grade transparency**
  Backends are explicitly documented and citable.

---

## Citation

If you use **gSpan**, please cite:

```bibtex
@inproceedings{yan2002gspan,
  title={gspan: Graph-based substructure pattern mining},
  author={Yan, Xifeng and Han, Jiawei},
  booktitle={Proceedings of the IEEE International Conference on Data Mining},
  pages={721--724},
  year={2002}
}
```

If you use **SoPaGraMi**, please cite:

```bibtex
@article{nguyen2020fast,
  title={Fast and scalable algorithms for mining subgraphs in a single large graph},
  author={Nguyen, Lam BQ and Vo, Bay and Le, Ngoc-Thao and Snasel, Vaclav and Zelinka, Ivan},
  journal={Engineering Applications of Artificial Intelligence},
  volume={90},
  pages={103539},
  year={2020}
}
```

To cite this library:

```bibtex
@misc{amure_submine,
  title  = {submine: A Unified Subgraph Mining Library},
  author = {Amure, Ridwan},
  year   = {2025},
  url    = {https://github.com/instabaines/submine}
}
```

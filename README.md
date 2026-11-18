# The family of transformed MSTs

This repository implements balanced transformed minimum spanning trees
(T-MSTs) for incremental MST maintenance and path-max queries.
We will release the accompanying paper separately.

Every structure is a header-only class in `include/tmst/` with the same
interface (checked by the concept `tmst::IncrementalMST` in `common.h`):

```cpp
using weight_type = W;                           // template parameter, default int
explicit T(int n);                               // n vertices, no edges
void Insert(int u, int v, W w);
std::optional<W> PathMax(int u, int v);          // nullopt if u == v or disconnected
bool Check() const;                              // verifies invariants for debugging
static std::string Name();
```

## Data structures and complexity

Height is the maximum vertex depth in edges after an operation on $n$ input
vertices. The $n-1$ entries give the general bound for a forest. For the
link-cut baseline, depth counts edges of the represented minimum spanning
forest (MSF). **am.** means amortized, **w.c.** means worst case, and
**w.h.p.** means with high probability.

| Implementation | Height | `Insert` | `PathMax` | Search-tree analogue |
|---|---|---|---|---|
| `OECForest<W>`<br>[oec_forest.h](include/tmst/oec_forest.h) | $O(\log n)$ w.h.p. | $O(\log n)$ w.h.p. | $O(\log n)$ w.h.p. | Treap |
| `WeightedDSU<W>`<br>[weighted_dsu.h](include/tmst/weighted_dsu.h) | $O(\log n)$ w.h.p. | $O(\log n)$ w.h.p. | $O(\log n)$ w.h.p. | Treap (OEC variant) |
| `StrictAMTree<LinkMode, W, kHash>`<br>[am_tree_strict.h](include/tmst/am_tree_strict.h) | $\le \log_{3/2} n$ | $O(\log n)$ am.<br>$O(\log^2 n)$ w.c. | $O(\log n)$ w.c. | Weight-balanced tree (BB[α]) |
| `LazyAMTree<LinkMode, W>`<br>[am_tree_lazy.h](include/tmst/am_tree_lazy.h) | $\le n-1$ | $O(\log n)$ am. | $O(\log n)$ am. | — |
| `BottleneckTree<W>`<br>[bottleneck_tree.h](include/tmst/bottleneck_tree.h) | $\le n-1$ | $O(\log n)$ am. | $O(\log n)$ am. | Splay tree |
| `RivalTree<LinkMode, W, kHash>`<br>[rival_tree.h](include/tmst/rival_tree.h) | $\le \log_\varphi(n+1)$ | $O(\log n)$ am.<br>$O(\log^2 n)$ w.c. | $O(\log n)$ w.c. | AVL tree |
| `LinkCutTree<W>`<br>[link_cut_tree.h](include/tmst/link_cut_tree.h) · baseline | $\le n-1$ (MSF) | $O(\log n)$ am. | $O(\log n)$ am. | — |

With $\varphi=(1+\sqrt{5})/2$, the strict AM-tree and Rival-tree height bounds
are approximately $1.71\log_2 n$ and $1.44\log_2(n+1)$, respectively.

The w.h.p. bounds assume uniformly random vertex priorities independent of the
operation sequence. They hold simultaneously over polynomially many insertions
and the intervening queries. `OECForest` and `WeightedDSU` default to seed `0`
for reproducibility. For a fixed vertex in a component of size $s$, the
expected depth is $H_s-1$ in the OEC-forest and at most $H_s-1$ in the weighted DSU, where
$H_s=\sum_{i=1}^{s}1/i$. [OEC depth law](report.md#thm-oec-exact) · [Weighted DSU bound](report.md#lem-wdsu)

`LinkMode` selects the perch-based link (`kPerch`) or the stitch-based link
(`kStitch`). The OEC-forest allows ties. The bounds for the two-sided stitch,
including weighted DSU, and for BottleneckTree's access rule assume distinct
edge weights. The insertion bounds for `StrictAMTree` and `RivalTree` use the
default array of child-list heads, which takes $O(n\log n)$ space. Setting
`kHash = true` stores the nonempty groups in $O(n)$ space, and insertion time
depends on hash probe lengths.

Every structure implements its own link. The shared path-max walk and forest
check are in [common.h](include/tmst/common.h). The weight type `W` defaults to
`int`. Every inserted weight must be less than `std::numeric_limits<W>::max()`,
the weight reserved for the virtual parent edge of a root.

## Technical report

The [technical report](report.md) contains the algorithms, theorems, and proofs for the OEC-forest, priority-based stitch, weighted DSU, and Rival-tree.

## Build

Requires CMake >= 3.18 and a C++20 compiler.

```sh
cmake -B build
cmake --build build -j
./build/test/tmst_test [seed]              # randomized tests against brute force
./build/bench/tmst_bench [n] [m] [q] [seed]
```

Status: the OEC-forest, weighted DSU, AM-tree (both versions, both links),
BottleneckTree, Rival-tree (both links), and Link-Cut Tree are implemented and
tested.

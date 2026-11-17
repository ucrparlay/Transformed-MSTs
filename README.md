# TMST-code

Implementations of the family of balanced T-MSTs (transformed minimum spanning
trees) for incremental minimum spanning trees, accompanying the paper
*The Family of Balanced T-MSTs for Incremental Minimum Spanning Trees*.

Every structure is a header-only class in `include/tmst/` with the same
interface (checked by the concept `tmst::IncrementalMST` in `common.h`):

```cpp
using weight_type = W;                           // template parameter, default int
explicit T(int n);                               // n vertices, no edges
void Insert(int u, int v, W w);
std::optional<W> PathMax(int u, int v);          // nullopt if u == v or disconnected
bool Check() const;                              // invariants, O(n), for debugging
static std::string Name();
```

| Header | Class | Paper |
|---|---|---|
| `oec_forest.h` | `OECForest<W>` | Section 6 |
| `weighted_dsu.h` | `WeightedDSU<W>` | Section 6 (remark) |
| `am_tree_lazy.h` | `LazyAMTree<LinkMode, W>` | Section 3 |
| `am_tree_strict.h` | `StrictAMTree<LinkMode, W, kHash>` | Section 3 |
| `bottleneck_tree.h` | `BottleneckTree<W>` | Section 4 |
| `rival_tree.h` | `RivalTree<LinkMode, W, kHash>` | Section 5 |
| `link_cut_tree.h` | `LinkCutTree<W>` | baseline |

`LinkMode` is `kPerch` (Algorithm 1) or `kStitch` (Algorithm 2).
`StrictAMTree` and `RivalTree` keep the children of every vertex grouped by
size or height, in an array of list heads by default, or in a hash table with
`kHash = true` (O(n) space instead of O(n log n)). Every structure implements its own link; only the
path-max walk and the forest check used by `Check()` are shared, in
`common.h`. The weight
type `W` defaults to `int`; `std::numeric_limits<W>::max()` is reserved for the
virtual parent edge of a root, so every inserted weight must be smaller than it.

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

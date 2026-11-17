#ifndef TMST_WEIGHTED_DSU_H_
#define TMST_WEIGHTED_DSU_H_

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <limits>
#include <numeric>
#include <optional>
#include <random>
#include <string>
#include <utility>
#include <vector>

#include "tmst/common.h"

namespace tmst {

// The weighted DSU [Codeforces blog 130107] (Remark in Section 6).
// The two-sided priority stitch, plus: the heaviest edge of the cycle is cut
// first, and pointers are shortcut during accesses.
template <class W = int>
class WeightedDSU {
 public:
  using weight_type = W;

  explicit WeightedDSU(int n, uint64_t seed = 0)
      : p(n, -1), wt(n, kInf), priority(n) {
    std::mt19937_64 rng(seed);
    std::iota(priority.begin(), priority.end(), 0);
    std::shuffle(priority.begin(), priority.end(), rng);
  }

  static std::string Name() { return "weighted-dsu"; }

  void Insert(int u, int v, W w) {
    assert(w < kInf);
    if (u == v) return;
    int heaviest = MainEdge(u, v);
    if (heaviest != -1) {
      if (w >= wt[heaviest]) return;
      p[heaviest] = -1;
      wt[heaviest] = kInf;
    }
    Stitch(u, v, w);
  }

  std::optional<W> PathMax(int u, int v) {
    if (u == v) return std::nullopt;
    int heaviest = MainEdge(u, v);
    if (heaviest == -1) return std::nullopt;
    return wt[heaviest];
  }

  bool Check() const {
    if (!CheckForest(p, wt)) return false;
    for (int v = 0; v < static_cast<int>(p.size()); v++) {
      if (p[v] != -1 && priority[p[v]] <= priority[v]) return false;
    }
    return true;
  }

 private:
  static constexpr W kInf = std::numeric_limits<W>::max();

  // Walks to the root in the forest restricted to edges of weight <= limit.
  int GetRoot(int v, W limit = kInf) {
    while (p[v] != -1 && wt[v] <= limit) {
      // The edge from v can skip any following edge of no greater weight.
      while (wt[p[v]] <= wt[v]) p[v] = p[p[v]];
      v = p[v];
    }
    return v;
  }

  // Returns the child endpoint of a maximum-weight edge on the path.
  int MainEdge(int u, int v) {
    int u_root = GetRoot(u);
    int v_root = GetRoot(v);
    if (u_root != v_root) return -1;

    // The shortcuts make weights strictly increase on both root paths.
    // Advance along the lighter edge until only the heaviest edge remains.
    while (p[u] != v && p[v] != u) {
      if (wt[u] < wt[v]) {
        u = p[u];
      } else {
        v = p[v];
      }
    }
    return p[u] == v ? u : v;
  }

  // Each pending edge joins different trees. Replacing a root's parent
  // finishes the stitch, since that root has no old edge to carry forward.
  void Stitch(int u, int v, W w) {
    while (u != -1) {
      u = GetRoot(u, w);
      v = GetRoot(v, w);
      if (priority[u] < priority[v]) std::swap(u, v);

      int parent = p[v];
      p[v] = u;
      std::swap(w, wt[v]);
      u = parent;
    }
  }

  std::vector<int> p;         // parent, -1 for a root
  std::vector<W> wt;          // parent edge weight, max() for a root
  std::vector<int> priority;  // permutation of 0, ..., n - 1
};

}  // namespace tmst

#endif  // TMST_WEIGHTED_DSU_H_

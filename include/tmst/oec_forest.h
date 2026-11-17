#ifndef TMST_OEC_FOREST_H_
#define TMST_OEC_FOREST_H_

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

// The OEC-forest [Song et al., SIGMOD'24] (Section 6 of the paper).
// Random vertex priorities, heap-ordered forest, one-sided stitch that only
// reads and writes the endpoint of lower priority.
template <class W = int>
class OECForest {
 public:
  using weight_type = W;

  explicit OECForest(int n, uint64_t seed = 0)
      : p(n, -1), wt(n, kInf), priority(n) {
    std::mt19937_64 rng(seed);
    std::iota(priority.begin(), priority.end(), 0);
    std::shuffle(priority.begin(), priority.end(), rng);
  }

  static std::string Name() { return "oec-forest"; }

  void Insert(int u, int v, W w) {
    assert(w < kInf);
    // Each step leaves only endpoints above v in priority. Writing v now
    // cannot affect the parent edges read by the remaining steps.
    while (v != -1 && u != v) {
      if (priority[u] < priority[v]) std::swap(u, v);
      int parent = p[v];
      if (w <= wt[v]) {
        p[v] = u;
        std::swap(w, wt[v]);
      }
      v = parent;
    }
  }

  std::optional<W> PathMax(int u, int v) const {
    std::optional<W> max_weight;
    while (u != v) {
      // The lower endpoint cannot be an ancestor of the higher endpoint.
      if (priority[u] < priority[v]) std::swap(u, v);
      if (p[v] == -1) return std::nullopt;
      if (!max_weight || *max_weight < wt[v]) max_weight = wt[v];
      v = p[v];
    }
    return max_weight;
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

  std::vector<int> p;         // parent, -1 for a root
  std::vector<W> wt;          // parent edge weight, max() for a root
  std::vector<int> priority;  // permutation of 0, ..., n - 1
};

}  // namespace tmst

#endif  // TMST_OEC_FOREST_H_

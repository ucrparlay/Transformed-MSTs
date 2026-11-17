#ifndef TMST_BOTTLENECK_TREE_H_
#define TMST_BOTTLENECK_TREE_H_

#include <cassert>
#include <limits>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "tmst/common.h"

namespace tmst {

// The BottleneckTree of Bicsi and Budau [arXiv:2611.05713] (Section 4 of the
// paper). Self-adjusting: Access(x) moves x to the root with a two-level rule,
// and the insertion is the perch-based link with Access as Perch. Stores only
// p and wt.
template <class W = int>
class BottleneckTree {
 public:
  using weight_type = W;

  explicit BottleneckTree(int n) : p(n, -1), wt(n, kInf) {}

  static std::string Name() { return "bottleneck-tree"; }

  void Insert(int u, int v, W w) {
    assert(w < kInf);
    if (u == v) return;
    Access(u);
    Access(v);
    // Now u is either a root with weight kInf or a direct child of v.
    if (w < wt[u]) {
      p[u] = v;
      wt[u] = w;
    }
  }

  std::optional<W> PathMax(int u, int v) {
    if (u == v) return std::nullopt;
    Access(u);
    Access(v);
    if (p[u] != v) return std::nullopt;
    return wt[u];
  }

  bool Check() const { return CheckForest(p, wt); }

 private:
  static constexpr W kInf = std::numeric_limits<W>::max();

  void Rotate(int x) {
    int y = p[x];
    p[x] = p[y];
    p[y] = x;
    std::swap(wt[x], wt[y]);
  }

  void Access(int x) {
    while (p[x] != -1) {
      int y = p[x];
      if (wt[x] <= wt[y]) {
        Rotate(x);
        continue;
      }

      // Here wt[y] < wt[x] < kInf, so y has a parent.
      int z = p[y];
      if (wt[y] > wt[z]) {
        p[y] = p[z];
      } else if (wt[z] <= wt[x]) {
        Rotate(y);
      }
      // y's edge is now light enough for x to skip it.
      p[x] = p[y];
    }
  }

  std::vector<int> p;  // parent, -1 for a root
  std::vector<W> wt;   // parent edge weight, max() for a root
};

}  // namespace tmst

#endif  // TMST_BOTTLENECK_TREE_H_

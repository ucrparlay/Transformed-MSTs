#ifndef TMST_AM_TREE_LAZY_H_
#define TMST_AM_TREE_LAZY_H_

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "tmst/common.h"

namespace tmst {

// The lazy AM-tree: UpwardCalibrate on both endpoints before every operation.
template <LinkMode kLink, class W = int>
class LazyAMTree {
 public:
  using weight_type = W;

  explicit LazyAMTree(int n) : p(n, -1), wt(n, kInf), s(n, 1) {}

  static std::string Name() {
    return std::string("am-lazy-") + ToString(kLink);
  }

  void Insert(int u, int v, W w) {
    assert(w < kInf);
    if (u == v) return;
    UpwardCalibrate(u);
    UpwardCalibrate(v);
    if constexpr (kLink == LinkMode::kPerch) {
      PerchLink(u, v, w);
    } else {
      StitchLink(u, v, w);
    }
  }

  std::optional<W> PathMax(int u, int v) {
    if (u == v) return std::nullopt;
    UpwardCalibrate(u);
    UpwardCalibrate(v);
    return PathMaxWalk(p, wt, u, v);
  }

  bool Check() const { return CheckForest(p, wt, &s); }

 private:
  static constexpr W kInf = std::numeric_limits<W>::max();

  // s(x) > 2/3 s(y)
  bool IsHeavy(int x, int y) const {
    return 3 * static_cast<int64_t>(s[x]) > 2 * static_cast<int64_t>(s[y]);
  }

  void Promote(int x) {
    int y = p[x];
    s[y] -= s[x];
    if (p[y] != -1 && wt[x] >= wt[y]) {  // shortcut
      p[x] = p[y];
    } else {  // rotate
      s[x] += s[y];
      std::swap(wt[x], wt[y]);
      p[x] = p[y];
      p[y] = x;
    }
  }

  void UpwardCalibrate(int x) {
    while (p[x] != -1) {
      int y = p[x];
      if (IsHeavy(x, y)) {
        Promote(x);
      } else {
        x = y;
      }
    }
  }

  // Promotes x until it is a root or a child of stop.
  void Perch(int x, int stop) {
    while (p[x] != -1 && p[x] != stop) Promote(x);
  }

  void PerchLink(int u, int v, W w) {
    Perch(v, -1);
    Perch(u, v);
    if (p[u] == v) {
      if (w < wt[u]) wt[u] = w;
    } else {
      if (s[u] > s[v]) std::swap(u, v);
      p[u] = v;
      wt[u] = w;
      s[v] += s[u];
    }
  }

  // Algorithm 2. s is not touched in the loop, so the rule sees the old
  // sizes; the left vertices are recorded and the sizes fixed afterwards.
  void StitchLink(int u, int v, W w) {
    left.clear();
    int a = u, b = v;
    bool attached = false;
    while (a != b && a != -1 && b != -1) {
      if (w > wt[a]) {
        if (attached) left.push_back({a, p[a], 0});
        a = p[a];
      } else if (w > wt[b]) {
        if (attached) left.push_back({b, p[b], 0});
        b = p[b];
      } else {
        if (s[b] < s[a] || (s[b] == s[a] && b < a)) std::swap(a, b);
        attached = true;
        int pa = p[a];
        W wa = wt[a];
        left.push_back({a, pa, 0});
        p[a] = b;
        wt[a] = w;
        a = pa;
        w = wa;
      }
    }
    if (!attached) return;

    // Disconnected: top's ancestors gain the attached tree.
    int top = a == b ? -1 : (a == -1 ? b : a);
    int old_top = top == -1 ? 0 : s[top];
    KeepLastLeft();
    RecomputeSizes();
    if (top != -1) {
      int delta = s[top] - old_top;
      for (int x = p[top]; x != -1; x = p[x]) s[x] += delta;
    }
  }

  struct Left {
    int v, p, s;  // vertex, old parent, old size
  };

  // Dedups left (keeping the last occurrence) and fills in the old sizes.
  void KeepLastLeft() {
    for (std::size_t i = left.size(); i-- > 0;) {
      Left& e = left[i];
      if (s[e.v] > 0) {
        e.s = s[e.v];
        s[e.v] = -e.s;
      } else {
        e.v = -1;
      }
    }
    std::size_t k = 0;
    for (std::size_t i = 0; i < left.size(); i++) {
      if (left[i].v == -1) continue;
      s[left[i].v] = left[i].s;
      left[k++] = left[i];
    }
    left.resize(k);
  }

  // left is in children-first order, so one pass suffices.
  void RecomputeSizes() {
    for (const Left& e : left) {
      if (e.p != -1) s[e.p] -= e.s;
    }
    for (const Left& e : left) {
      if (p[e.v] != -1) s[p[e.v]] += s[e.v];
    }
  }

  std::vector<int> p;  // parent, -1 for a root
  std::vector<W> wt;   // parent edge weight, max() for a root
  std::vector<int> s;  // subtree size
  std::vector<Left> left;
};

}  // namespace tmst

#endif  // TMST_AM_TREE_LAZY_H_

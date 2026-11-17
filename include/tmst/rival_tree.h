#ifndef TMST_RIVAL_TREE_H_
#define TMST_RIVAL_TREE_H_

#include <algorithm>
#include <bit>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "tmst/common.h"

namespace tmst {

// The children of every vertex, grouped by stored height. Group k of y is a
// list from head[y * keys + k] via next/prev; key[c] is c's group, and bit k
// of mask[y] is set iff group k of y is nonempty.
class ChildrenByHeight {
 public:
  ChildrenByHeight(int n, int keys)
      : keys(keys),
        mask(n, 0),
        head(static_cast<std::size_t>(n) * keys, -1),
        next(n, -1),
        prev(n, -1),
        key(n, 0) {
    assert(keys <= 64);
  }

  void Add(int y, int c, int k) {
    assert(0 <= k && k < keys);
    key[c] = static_cast<uint8_t>(k);
    int& h = Head(y, k);
    prev[c] = -1;
    next[c] = h;
    if (h != -1) prev[h] = c;
    h = c;
    mask[y] |= uint64_t{1} << k;
  }

  void Remove(int y, int c) {
    int k = key[c];
    if (prev[c] != -1) {
      next[prev[c]] = next[c];
    } else {
      Head(y, k) = next[c];
      if (next[c] == -1) mask[y] &= ~(uint64_t{1} << k);
    }
    if (next[c] != -1) prev[next[c]] = prev[c];
  }

  // One plus the largest key of a child of y, or 0.
  int Height(int y) const { return std::bit_width(mask[y]); }

  // The heavy child of y if y is unbalanced, or -1: y is unbalanced iff its
  // largest key k >= 2 has a single child and group k - 1 is empty.
  int Heavy(int y) const {
    int k = std::bit_width(mask[y]) - 1;
    if (k < 2 || (mask[y] >> (k - 1) & 1) != 0) return -1;
    int c = Head(y, k);
    return next[c] == -1 ? c : -1;
  }

  // Checks that the lists hold exactly the children given by p, grouped by h.
  bool Check(const std::vector<int>& p, const std::vector<int>& h) const {
    const int n = static_cast<int>(p.size());
    std::vector<int> num_children(n, 0), listed(n, 0);
    for (int x = 0; x < n; x++) {
      if (p[x] != -1) num_children[p[x]]++;
    }
    for (int y = 0; y < n; y++) {
      uint64_t m = 0;
      for (int k = 0; k < keys; k++) {
        int last = -1;
        for (int c = Head(y, k); c != -1; c = next[c]) {
          if (p[c] != y || key[c] != k || h[c] != k) return false;
          if (prev[c] != last) return false;
          if (++listed[y] > num_children[y]) return false;
          last = c;
        }
        if (Head(y, k) != -1) m |= uint64_t{1} << k;
      }
      if (m != mask[y] || listed[y] != num_children[y]) return false;
    }
    return true;
  }

 private:
  int& Head(int y, int k) {
    return head[static_cast<std::size_t>(y) * keys + k];
  }
  int Head(int y, int k) const {
    return head[static_cast<std::size_t>(y) * keys + k];
  }

  int keys;
  std::vector<uint64_t> mask;
  std::vector<int> head;
  std::vector<int> next, prev;
  std::vector<uint8_t> key;
};

// The same as ChildrenByHeight, but the nonempty heads are kept in a hash
// table keyed by (y, k), so the space is O(n) instead of O(n log n). Open
// addressing with linear probing and backward-shift deletion; the table has
// at least 2n slots and at most n - 1 entries.
class HashedChildrenByHeight {
 public:
  HashedChildrenByHeight(int n, int keys)
      : keys(keys),
        bits(std::bit_width(std::bit_ceil(2 * static_cast<std::size_t>(
                 std::max(n, 1)))) - 1),
        mod((std::size_t{1} << bits) - 1),
        slots(mod + 1),
        mask(n, 0),
        next(n, -1),
        prev(n, -1),
        key(n, 0) {
    assert(keys <= 64);
  }

  void Add(int y, int c, int k) {
    assert(0 <= k && k < keys);
    key[c] = static_cast<uint8_t>(k);
    Slot& e = slots[Find(y, k)];
    prev[c] = -1;
    if (e.y == -1) {
      e = {y, c, static_cast<uint8_t>(k)};
      next[c] = -1;
      mask[y] |= uint64_t{1} << k;
    } else {
      next[c] = e.head;
      prev[e.head] = c;
      e.head = c;
    }
  }

  void Remove(int y, int c) {
    int k = key[c];
    if (prev[c] != -1) {
      next[prev[c]] = next[c];
    } else if (next[c] != -1) {
      slots[Find(y, k)].head = next[c];
    } else {
      Erase(Find(y, k));
      mask[y] &= ~(uint64_t{1} << k);
    }
    if (next[c] != -1) prev[next[c]] = prev[c];
  }

  int Height(int y) const { return std::bit_width(mask[y]); }

  // As ChildrenByHeight::Heavy; looks up the table only if group k - 1 is
  // empty.
  int Heavy(int y) const {
    int k = std::bit_width(mask[y]) - 1;
    if (k < 2 || (mask[y] >> (k - 1) & 1) != 0) return -1;
    int c = slots[Find(y, k)].head;
    return next[c] == -1 ? c : -1;
  }

  // Checks that the lists hold exactly the children given by p, grouped by h.
  bool Check(const std::vector<int>& p, const std::vector<int>& h) const {
    const int n = static_cast<int>(p.size());
    std::vector<int> num_children(n, 0), listed(n, 0);
    for (int x = 0; x < n; x++) {
      if (p[x] != -1) num_children[p[x]]++;
    }
    std::size_t entries = 0, nonempty = 0;
    for (const Slot& e : slots) entries += e.y != -1;
    for (int y = 0; y < n; y++) {
      for (int k = 0; k < keys; k++) {
        const Slot& e = slots[Find(y, k)];
        bool listed_here = (mask[y] >> k & 1) != 0;
        if ((e.y != -1) != listed_here) return false;
        if (!listed_here) continue;
        nonempty++;
        int last = -1;
        for (int c = e.head; c != -1; c = next[c]) {
          if (p[c] != y || key[c] != k || h[c] != k) return false;
          if (prev[c] != last) return false;
          if (++listed[y] > num_children[y]) return false;
          last = c;
        }
      }
      if ((keys < 64 && mask[y] >> keys != 0) ||
          listed[y] != num_children[y]) {
        return false;
      }
    }
    return entries == nonempty;
  }

 private:
  struct Slot {
    int y = -1;  // -1 for an empty slot
    int head = -1;
    uint8_t k = 0;
  };

  std::size_t Home(int y, int k) const {
    uint64_t x = static_cast<uint64_t>(y) * keys + k;
    return static_cast<std::size_t>((x * 0x9E3779B97F4A7C15ull) >>
                                    (64 - bits));
  }

  // The slot of (y, k), or the empty slot where it would go.
  std::size_t Find(int y, int k) const {
    std::size_t i = Home(y, k);
    while (slots[i].y != -1 && (slots[i].y != y || slots[i].k != k)) {
      i = (i + 1) & mod;
    }
    return i;
  }

  void Erase(std::size_t i) {
    for (std::size_t j = (i + 1) & mod; slots[j].y != -1; j = (j + 1) & mod) {
      std::size_t h = Home(slots[j].y, slots[j].k);
      if (((j - h) & mod) >= ((j - i) & mod)) {
        slots[i] = slots[j];
        i = j;
      }
    }
    slots[i] = Slot{};
  }

  int keys;
  int bits;
  std::size_t mod;
  std::vector<Slot> slots;
  std::vector<uint64_t> mask;
  std::vector<int> next, prev;
  std::vector<uint8_t> key;
};

// The Rival-tree (Algorithm 3 of the paper): every vertex of height k >= 3
// has a child of height k - 1 and a rival of height at least k - 2. After the
// link, Rebalance fixes the vertices on the old root paths of u and v.
// kHash selects HashedChildrenByHeight.
template <LinkMode kLink, class W = int, bool kHash = false>
class RivalTree {
 public:
  using weight_type = W;

  explicit RivalTree(int n)
      : p(n, -1),
        wt(n, kInf),
        h(n, 0),
        children(n, MaxHeight(n) + 2),
        mark(n, 0),
        depth(n, 0) {}

  static std::string Name() {
    return std::string("rival-") + ToString(kLink) + (kHash ? "-hash" : "");
  }

  void Insert(int u, int v, W w) {
    assert(w < kInf);
    if (u == v) return;
    CollectPaths(u, v);
    if constexpr (kLink == LinkMode::kPerch) {
      PerchLink(u, v, w);
    } else {
      if (!StitchLink(u, v, w)) return;
    }
    Rebalance();
  }

  std::optional<W> PathMax(int u, int v) const {
    return PathMaxWalk(p, wt, u, v);
  }

  bool Check() const {
    std::vector<int> d;
    if (!CheckForest(p, wt, nullptr, &d)) return false;
    const int n = static_cast<int>(p.size());
    // Heights from the deepest vertices up, and the two tallest children.
    std::vector<int> order(n);
    for (int x = 0; x < n; x++) order[x] = x;
    std::sort(order.begin(), order.end(),
              [&](int a, int b) { return d[a] > d[b]; });
    std::vector<int> height(n, 0), c1(n, -1), c2(n, -1);
    for (int x : order) {
      height[x] = c1[x] + 1;
      if (h[x] != height[x]) return false;
      int y = p[x];
      if (y == -1) continue;
      if (height[x] > c1[y]) {
        c2[y] = c1[y];
        c1[y] = height[x];
      } else if (height[x] > c2[y]) {
        c2[y] = height[x];
      }
    }
    for (int y = 0; y < n; y++) {
      if (c1[y] >= 2 && c2[y] < c1[y] - 1) return false;  // rival rule
    }
    return children.Check(p, h);
  }

 private:
  static constexpr W kInf = std::numeric_limits<W>::max();

  // The largest height of a balanced tree with at most n vertices.
  static int MaxHeight(int n) {
    int64_t a = 2, b = 3;  // minimum sizes of heights k - 1 and k
    int k = n >= 3 ? 2 : std::max(n - 1, 0);
    while (n >= 3 && a + b + 1 <= n) {
      a = std::exchange(b, a + b + 1);
      k++;
    }
    return k;
  }

  void Promote(int x) {
    int y = p[x], z = p[y];
    children.Remove(y, x);
    if (z != -1 && wt[x] >= wt[y]) {  // shortcut
      p[x] = z;
      children.Add(z, x, h[x]);
    } else {  // rotate
      if (z != -1) children.Remove(z, y);
      std::swap(wt[x], wt[y]);
      p[x] = z;
      p[y] = x;
      if (z != -1) children.Add(z, x, h[x]);
      children.Add(x, y, h[y]);
    }
  }

  // Recomputes h(y) from its children and moves y to the right group.
  void UpdateHeight(int y) {
    int k = children.Height(y);
    if (k == h[y]) return;
    if (p[y] != -1) children.Remove(p[y], y);
    h[y] = k;
    if (p[y] != -1) children.Add(p[y], y, k);
  }

  void Fix(int y) {
    for (int c; (c = children.Heavy(y)) != -1;) {
      Promote(c);
      UpdateHeight(y);
    }
  }

  // path = u, v and their ancestors, without repeats.
  void CollectPaths(int u, int v) {
    path.clear();
    stamp++;
    for (int x = u; x != -1; x = p[x]) {
      mark[x] = stamp;
      path.push_back(x);
    }
    for (int x = v; x != -1 && mark[x] != stamp; x = p[x]) {
      mark[x] = stamp;
      path.push_back(x);
    }
  }

  void PerchLink(int u, int v, W w) {
    while (p[u] != -1) Promote(u);
    while (p[v] != -1) Promote(v);
    if (p[u] == v) {
      if (w < wt[u]) wt[u] = w;
    } else {
      p[u] = v;
      wt[u] = w;
      children.Add(v, u, h[u]);
    }
  }

  // Algorithm 2 with the height rule; h is not updated in the loop, so the
  // rule sees the old heights. Returns false if the forest did not change.
  bool StitchLink(int u, int v, W w) {
    int a = u, b = v;
    bool attached = false;
    while (a != b && a != -1 && b != -1) {
      if (w > wt[a]) {
        a = p[a];
      } else if (w > wt[b]) {
        b = p[b];
      } else {
        if (h[b] < h[a]) std::swap(a, b);
        attached = true;
        int pa = p[a];
        W wa = wt[a];
        if (pa != -1) children.Remove(pa, a);
        p[a] = b;
        wt[a] = w;
        children.Add(b, a, h[a]);
        a = pa;
        w = wa;
      }
    }
    return attached;
  }

  // Fixes the vertices of path in order of non-increasing depth after the
  // link. path is closed under parents, so the depths take O(|path|) time.
  void Rebalance() {
    int max_depth = 0;
    for (int x : path) {
      if (mark[x] == -stamp) continue;
      int y = x;
      stack.clear();
      while (y != -1 && mark[y] != -stamp) {
        stack.push_back(y);
        y = p[y];
      }
      int d = y == -1 ? -1 : depth[y];
      while (!stack.empty()) {
        depth[stack.back()] = ++d;
        mark[stack.back()] = -stamp;
        stack.pop_back();
      }
      max_depth = std::max(max_depth, d);
    }
    count.assign(max_depth + 2, 0);
    for (int x : path) count[max_depth - depth[x] + 1]++;
    for (int i = 1; i <= max_depth + 1; i++) count[i] += count[i - 1];
    order.resize(path.size());
    for (int x : path) order[count[max_depth - depth[x]]++] = x;
    for (int y : order) {
      UpdateHeight(y);
      Fix(y);
    }
  }

  std::vector<int> p;  // parent, -1 for a root
  std::vector<W> wt;   // parent edge weight, max() for a root
  std::vector<int> h;  // stored height, 0 for a leaf
  std::conditional_t<kHash, HashedChildrenByHeight, ChildrenByHeight>
      children;

  // Scratch space for Insert. mark[x] is stamp if x is in path, and -stamp
  // once depth[x] is known.
  std::vector<int> mark, depth;
  int stamp = 0;
  std::vector<int> path, stack, count, order;
};

}  // namespace tmst

#endif  // TMST_RIVAL_TREE_H_

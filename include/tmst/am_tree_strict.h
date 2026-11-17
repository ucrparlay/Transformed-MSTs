#ifndef TMST_AM_TREE_STRICT_H_
#define TMST_AM_TREE_STRICT_H_

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

// The children of every vertex, grouped by floor(log2 size). Group g of y is
// a list from head[y * groups + g] via next/prev; group[c] is c's group, and
// bit g of mask[y] is set iff group g of y is nonempty.
class ChildrenBySize {
 public:
  explicit ChildrenBySize(int n)
      : groups(std::max(1, static_cast<int>(std::bit_width(
                               static_cast<unsigned>(std::max(n, 0)))))),
        mask(n, 0),
        head(static_cast<std::size_t>(n) * groups, -1),
        next(n, -1),
        prev(n, -1),
        group(n, 0) {}

  void Add(int y, int c, int size) {
    int g = group[c] = static_cast<uint8_t>(Group(size));
    int& h = Head(y, g);
    prev[c] = -1;
    next[c] = h;
    if (h != -1) prev[h] = c;
    h = c;
    mask[y] |= uint32_t{1} << g;
  }

  void Remove(int y, int c) {
    int g = group[c];
    if (prev[c] != -1) {
      next[prev[c]] = next[c];
    } else {
      Head(y, g) = next[c];
      if (next[c] == -1) mask[y] &= ~(uint32_t{1} << g);
    }
    if (next[c] != -1) prev[next[c]] = prev[c];
  }

  // The child alone in the highest group, or -1; no other child can be heavy.
  int Candidate(int y) const {
    if (mask[y] == 0) return -1;
    int c = Head(y, std::bit_width(mask[y]) - 1);
    return next[c] == -1 ? c : -1;
  }

  // Checks that the lists hold exactly the children given by p, grouped by s.
  bool Check(const std::vector<int>& p, const std::vector<int>& s) const {
    const int n = static_cast<int>(p.size());
    std::vector<int> num_children(n, 0), listed(n, 0);
    for (int x = 0; x < n; x++) {
      if (p[x] != -1) num_children[p[x]]++;
    }
    for (int y = 0; y < n; y++) {
      uint32_t m = 0;
      for (int g = 0; g < groups; g++) {
        int last = -1;
        for (int c = Head(y, g); c != -1; c = next[c]) {
          if (p[c] != y || group[c] != g || Group(s[c]) != g) return false;
          if (prev[c] != last) return false;
          if (++listed[y] > num_children[y]) return false;
          last = c;
        }
        if (Head(y, g) != -1) m |= uint32_t{1} << g;
      }
      if (m != mask[y] || listed[y] != num_children[y]) return false;
    }
    return true;
  }

 private:
  static int Group(int size) {
    return std::bit_width(static_cast<unsigned>(size)) - 1;
  }

  int& Head(int y, int g) {
    return head[static_cast<std::size_t>(y) * groups + g];
  }
  int Head(int y, int g) const {
    return head[static_cast<std::size_t>(y) * groups + g];
  }

  int groups;
  std::vector<uint32_t> mask;
  std::vector<int> head;
  std::vector<int> next, prev;
  std::vector<uint8_t> group;
};

// The same as ChildrenBySize, but the nonempty heads are kept in a hash table
// keyed by (y, g), so the space is O(n) instead of O(n log n). Open
// addressing with linear probing and backward-shift deletion; the table has
// at least 2n slots and at most n - 1 entries.
class HashedChildrenBySize {
 public:
  explicit HashedChildrenBySize(int n)
      : groups(std::max(1, static_cast<int>(std::bit_width(
                               static_cast<unsigned>(std::max(n, 0)))))),
        bits(std::bit_width(std::bit_ceil(2 * static_cast<std::size_t>(
                 std::max(n, 1)))) - 1),
        mod((std::size_t{1} << bits) - 1),
        slots(mod + 1),
        mask(n, 0),
        next(n, -1),
        prev(n, -1),
        group(n, 0) {}

  void Add(int y, int c, int size) {
    int g = group[c] = static_cast<uint8_t>(Group(size));
    Slot& e = slots[Find(y, g)];
    prev[c] = -1;
    if (e.y == -1) {
      e = {y, c, static_cast<uint8_t>(g)};
      next[c] = -1;
      mask[y] |= uint32_t{1} << g;
    } else {
      next[c] = e.head;
      prev[e.head] = c;
      e.head = c;
    }
  }

  void Remove(int y, int c) {
    int g = group[c];
    if (prev[c] != -1) {
      next[prev[c]] = next[c];
    } else if (next[c] != -1) {
      slots[Find(y, g)].head = next[c];
    } else {
      Erase(Find(y, g));
      mask[y] &= ~(uint32_t{1} << g);
    }
    if (next[c] != -1) prev[next[c]] = prev[c];
  }

  // The child alone in the highest group, or -1; no other child can be heavy.
  int Candidate(int y) const {
    if (mask[y] == 0) return -1;
    int c = slots[Find(y, std::bit_width(mask[y]) - 1)].head;
    return next[c] == -1 ? c : -1;
  }

  // Checks that the lists hold exactly the children given by p, grouped by s.
  bool Check(const std::vector<int>& p, const std::vector<int>& s) const {
    const int n = static_cast<int>(p.size());
    std::vector<int> num_children(n, 0), listed(n, 0);
    for (int x = 0; x < n; x++) {
      if (p[x] != -1) num_children[p[x]]++;
    }
    std::size_t entries = 0, nonempty = 0;
    for (const Slot& e : slots) entries += e.y != -1;
    for (int y = 0; y < n; y++) {
      for (int g = 0; g < groups; g++) {
        const Slot& e = slots[Find(y, g)];
        bool listed_here = (mask[y] >> g & 1) != 0;
        if ((e.y != -1) != listed_here) return false;
        if (!listed_here) continue;
        nonempty++;
        int last = -1;
        for (int c = e.head; c != -1; c = next[c]) {
          if (p[c] != y || group[c] != g || Group(s[c]) != g) return false;
          if (prev[c] != last) return false;
          if (++listed[y] > num_children[y]) return false;
          last = c;
        }
      }
      if (mask[y] >> groups != 0 || listed[y] != num_children[y]) return false;
    }
    return entries == nonempty;
  }

 private:
  struct Slot {
    int y = -1;  // -1 for an empty slot
    int head = -1;
    uint8_t g = 0;
  };

  static int Group(int size) {
    return std::bit_width(static_cast<unsigned>(size)) - 1;
  }

  std::size_t Home(int y, int g) const {
    uint64_t key = static_cast<uint64_t>(y) * groups + g;
    return static_cast<std::size_t>((key * 0x9E3779B97F4A7C15ull) >>
                                    (64 - bits));
  }

  // The slot of (y, g), or the empty slot where it would go.
  std::size_t Find(int y, int g) const {
    std::size_t i = Home(y, g);
    while (slots[i].y != -1 && (slots[i].y != y || slots[i].g != g)) {
      i = (i + 1) & mod;
    }
    return i;
  }

  void Erase(std::size_t i) {
    for (std::size_t j = (i + 1) & mod; slots[j].y != -1; j = (j + 1) & mod) {
      std::size_t h = Home(slots[j].y, slots[j].g);
      if (((j - h) & mod) >= ((j - i) & mod)) {
        slots[i] = slots[j];
        i = j;
      }
    }
    slots[i] = Slot{};
  }

  int groups;
  int bits;
  std::size_t mod;
  std::vector<Slot> slots;
  std::vector<uint32_t> mask;
  std::vector<int> next, prev;
  std::vector<uint8_t> group;
};

// The strict AM-tree: DownwardCalibrate after every link, so the rule always
// holds. kHash selects HashedChildrenBySize.
template <LinkMode kLink, class W = int, bool kHash = false>
class StrictAMTree {
 public:
  using weight_type = W;

  explicit StrictAMTree(int n) : p(n, -1), wt(n, kInf), s(n, 1), children(n) {}

  static std::string Name() {
    return std::string("am-strict-") + ToString(kLink) +
           (kHash ? "-hash" : "");
  }

  void Insert(int u, int v, W w) {
    assert(w < kInf);
    if (u == v) return;
    if constexpr (kLink == LinkMode::kPerch) {
      PerchLink(u, v, w);
    } else {
      StitchLink(u, v, w);
    }
  }

  std::optional<W> PathMax(int u, int v) const {
    return PathMaxWalk(p, wt, u, v);
  }

  bool Check() const {
    if (!CheckForest(p, wt, &s)) return false;
    for (int x = 0; x < static_cast<int>(p.size()); x++) {
      if (p[x] != -1 && IsHeavy(x, p[x])) return false;
    }
    return children.Check(p, s);
  }

 private:
  static constexpr W kInf = std::numeric_limits<W>::max();

  void Link(int c) {
    if (p[c] != -1) children.Add(p[c], c, s[c]);
  }
  void Unlink(int c) {
    if (p[c] != -1) children.Remove(p[c], c);
  }

  // s(x) > 2/3 s(y)
  bool IsHeavy(int x, int y) const {
    return 3 * static_cast<int64_t>(s[x]) > 2 * static_cast<int64_t>(s[y]);
  }

  int HeavyChild(int y) const {
    int c = children.Candidate(y);
    return c != -1 && IsHeavy(c, y) ? c : -1;
  }

  void Promote(int x) {
    int y = p[x], z = p[y];
    Unlink(x);
    Unlink(y);
    s[y] -= s[x];
    if (z != -1 && wt[x] >= wt[y]) {  // shortcut
      p[x] = z;
    } else {  // rotate
      s[x] += s[y];
      std::swap(wt[x], wt[y]);
      p[x] = z;
      p[y] = x;
    }
    Link(x);
    Link(y);
  }

  void DownwardCalibrate(int y) {
    for (int c; (c = HeavyChild(y)) != -1;) Promote(c);
  }

  // Promotes x until it is a root or a child of stop.
  void Perch(int x, int stop) {
    while (p[x] != -1 && p[x] != stop) {
      int y = p[x];
      Promote(x);
      DownwardCalibrate(y);
    }
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
      Link(u);
    }
    for (int x = u; x != -1; x = p[x]) DownwardCalibrate(x);
  }

  // Algorithm 2, as in LazyAMTree. The children lists are not touched in the
  // loop, so every vertex is still listed under its old parent afterwards.
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

    // Connected: the pointers met at z. Disconnected: top's ancestors gain
    // the attached tree.
    int z = a == b ? a : -1;
    int top = a == b ? -1 : (a == -1 ? b : a);
    KeepLastLeft();

    for (const Left& e : left) {
      if (e.p != -1) children.Remove(e.p, e.v);
    }
    for (int x = top; x != -1; x = p[x]) Unlink(x);
    int old_top = top == -1 ? 0 : s[top];
    RecomputeSizes();
    if (top != -1) {
      int delta = s[top] - old_top;
      for (int x = p[top]; x != -1; x = p[x]) s[x] += delta;
    }
    for (const Left& e : left) Link(e.v);
    for (int x = top; x != -1; x = p[x]) Link(x);

    // Children first.
    for (const Left& e : left) DownwardCalibrate(e.v);
    if (z != -1) DownwardCalibrate(z);
    for (int x = top; x != -1; x = p[x]) DownwardCalibrate(x);
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
  std::conditional_t<kHash, HashedChildrenBySize, ChildrenBySize> children;
  std::vector<Left> left;
};

}  // namespace tmst

#endif  // TMST_AM_TREE_STRICT_H_

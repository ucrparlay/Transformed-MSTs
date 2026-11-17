#ifndef TMST_COMMON_H_
#define TMST_COMMON_H_

#include <algorithm>
#include <concepts>
#include <cstdint>
#include <limits>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace tmst {

// The two link procedures of the paper (Algorithms 1 and 2).
enum class LinkMode { kPerch, kStitch };

constexpr const char* ToString(LinkMode mode) {
  switch (mode) {
    case LinkMode::kPerch:
      return "perch";
    case LinkMode::kStitch:
      return "stitch";
  }
  return "?";
}

// The interface of every structure; PathMax returns nullopt if u == v or
// u and v are disconnected.
template <class T>
concept IncrementalMST =
    std::constructible_from<T, int> &&
    requires(T t, const T ct, int u, int v, typename T::weight_type w) {
      { t.Insert(u, v, w) } -> std::same_as<void>;
      {
        t.PathMax(u, v)
      } -> std::same_as<std::optional<typename T::weight_type>>;
      { ct.Check() } -> std::same_as<bool>;
      { T::Name() } -> std::convertible_to<std::string>;
    };

// Walks from u and v up to their LCA.
template <class W>
std::optional<W> PathMaxWalk(const std::vector<int>& p,
                             const std::vector<W>& wt, int u, int v) {
  if (u == v) return std::nullopt;
  int du = 0, dv = 0, ru = u, rv = v;
  while (p[ru] != -1) ru = p[ru], du++;
  while (p[rv] != -1) rv = p[rv], dv++;
  if (ru != rv) return std::nullopt;
  W res = std::numeric_limits<W>::lowest();
  for (; du > dv; du--) res = std::max(res, wt[u]), u = p[u];
  for (; dv > du; dv--) res = std::max(res, wt[v]), v = p[v];
  while (u != v) {
    res = std::max(res, std::max(wt[u], wt[v]));
    u = p[u];
    v = p[v];
  }
  return res;
}

// Checks that p is a forest whose roots are exactly the vertices of weight
// max(), and that s (if given) holds the subtree sizes.
template <class W>
bool CheckForest(const std::vector<int>& p, const std::vector<W>& wt,
                 const std::vector<int>* s = nullptr,
                 std::vector<int>* depth = nullptr) {
  const int n = static_cast<int>(p.size());
  if (static_cast<int>(wt.size()) != n) return false;
  if (s != nullptr && static_cast<int>(s->size()) != n) return false;
  for (int v = 0; v < n; v++) {
    if (p[v] < -1 || p[v] >= n) return false;
    if ((p[v] == -1) != (wt[v] == std::numeric_limits<W>::max())) {
      return false;
    }
  }

  // state: 0 = new, 1 = on the stack, 2 = done
  std::vector<int> d(n, 0);
  std::vector<uint8_t> state(n, 0);
  std::vector<int> stack;
  for (int v = 0; v < n; v++) {
    int x = v;
    while (state[x] == 0) {
      state[x] = 1;
      stack.push_back(x);
      if (p[x] == -1) break;
      x = p[x];
    }
    int dx;
    if (state[x] == 1 && p[x] == -1) {  // reached a new root
      stack.pop_back();
      state[x] = 2;
      d[x] = dx = 0;
    } else if (state[x] == 1) {
      return false;  // cycle
    } else {
      dx = d[x];
    }
    while (!stack.empty()) {
      int y = stack.back();
      stack.pop_back();
      state[y] = 2;
      d[y] = ++dx;
    }
  }

  if (s != nullptr) {
    int max_d = 0;
    for (int v = 0; v < n; v++) max_d = std::max(max_d, d[v]);
    std::vector<int> start(max_d + 2, 0), order(n);
    for (int v = 0; v < n; v++) start[d[v] + 1]++;
    for (int i = 1; i <= max_d + 1; i++) start[i] += start[i - 1];
    for (int v = 0; v < n; v++) order[start[d[v]]++] = v;
    std::vector<int64_t> size(n, 1);
    for (int i = n - 1; i >= 0; i--) {
      int v = order[i];
      if (p[v] != -1) size[p[v]] += size[v];
    }
    for (int v = 0; v < n; v++) {
      if (size[v] != (*s)[v]) return false;
    }
  }
  if (depth != nullptr) *depth = std::move(d);
  return true;
}

}  // namespace tmst

#endif  // TMST_COMMON_H_

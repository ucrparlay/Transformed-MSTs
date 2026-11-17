// Randomized tests: every implemented structure against a brute-force
// incremental MST, with Check() after every operation.
//
// Usage: tmst_test [seed]

#include <algorithm>
#include <climits>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <optional>
#include <random>
#include <string>
#include <tuple>
#include <vector>

#include "tmst/tmst.h"

namespace {

// Keeps the edges of an MST explicitly; every operation walks the forest.
template <class W>
class Reference {
 public:
  explicit Reference(int n) : adj(n) {}

  void Insert(int u, int v, W w) {
    if (u == v) return;
    auto path = Path(u, v);
    if (!path) {
      AddEdge(u, v, w);
      return;
    }
    int heaviest = (*path)[0];
    for (int id : *path) {
      if (edges[id].w > edges[heaviest].w) heaviest = id;
    }
    if (w < edges[heaviest].w) {
      RemoveEdge(heaviest);
      AddEdge(u, v, w);
    }
  }

  std::optional<W> PathMax(int u, int v) const {
    if (u == v) return std::nullopt;
    auto path = Path(u, v);
    if (!path) return std::nullopt;
    W res = edges[(*path)[0]].w;
    for (int id : *path) res = std::max(res, edges[id].w);
    return res;
  }

 private:
  struct Edge {
    int u, v;
    W w;
  };

  void AddEdge(int u, int v, W w) {
    int id = static_cast<int>(edges.size());
    edges.push_back({u, v, w});
    adj[u].push_back(id);
    adj[v].push_back(id);
  }

  void RemoveEdge(int id) {
    for (int x : {edges[id].u, edges[id].v}) {
      adj[x].erase(std::find(adj[x].begin(), adj[x].end(), id));
    }
  }

  // The edge ids on the tree path from u to v, or nullopt if disconnected.
  std::optional<std::vector<int>> Path(int u, int v) const {
    std::vector<int> via(adj.size(), -2);
    std::vector<int> queue{u};
    via[u] = -1;
    for (std::size_t i = 0; i < queue.size(); i++) {
      int x = queue[i];
      for (int id : adj[x]) {
        int y = edges[id].u ^ edges[id].v ^ x;
        if (via[y] == -2) via[y] = id, queue.push_back(y);
      }
    }
    if (via[v] == -2) return std::nullopt;
    std::vector<int> path;
    for (int x = v; x != u;) {
      int id = via[x];
      path.push_back(id);
      x = edges[id].u ^ edges[id].v ^ x;
    }
    return path;
  }

  std::vector<std::vector<int>> adj;
  std::vector<Edge> edges;
};

int failures = 0;

std::string Show(const std::optional<int64_t>& x) {
  return x ? std::to_string(*x) : "none";
}

// Runs a sequence of operations; an operation (u, v, w) with w < 0 is a
// query. Compares all pairs at the end if n is small.
template <class T>
bool Run(const char* label, int n,
         const std::vector<std::tuple<int, int, int64_t>>& ops) {
  using W = typename T::weight_type;
  T t(n);
  Reference<W> ref(n);
  auto fail = [&](std::size_t i, const std::string& what) {
    std::fprintf(stderr, "FAIL %s [%s] n=%d op %zu: %s\n", T::Name().c_str(),
                 label, n, i, what.c_str());
    failures++;
    return false;
  };
  auto compare = [&](std::size_t i, int u, int v) {
    std::optional<int64_t> got, want;
    if (auto r = t.PathMax(u, v)) got = static_cast<int64_t>(*r);
    if (auto r = ref.PathMax(u, v)) want = static_cast<int64_t>(*r);
    if (got != want) {
      return fail(i, "PathMax(" + std::to_string(u) + "," + std::to_string(v) +
                         ") = " + Show(got) + ", want " + Show(want));
    }
    return true;
  };
  for (std::size_t i = 0; i < ops.size(); i++) {
    auto [u, v, w] = ops[i];
    if (w < 0) {
      if (!compare(i, u, v)) return false;
    } else {
      t.Insert(u, v, static_cast<W>(w));
      ref.Insert(u, v, static_cast<W>(w));
    }
    if (!t.Check()) return fail(i, "Check() failed");
  }
  if (n <= 64) {
    for (int u = 0; u < n; u++) {
      for (int v = 0; v < n; v++) {
        if (!compare(ops.size(), u, v)) return false;
      }
    }
  }
  return true;
}

using Ops = std::vector<std::tuple<int, int, int64_t>>;

Ops RandomOps(int n, int m, int64_t max_w, double query_ratio,
              std::mt19937_64& rng) {
  std::uniform_int_distribution<int> vertex(0, n - 1);
  std::uniform_int_distribution<int64_t> weight(0, max_w);
  std::bernoulli_distribution query(query_ratio);
  Ops ops;
  for (int i = 0; i < m; i++) {
    int u = vertex(rng), v = vertex(rng);
    ops.emplace_back(u, v, query(rng) ? -1 : weight(rng));
  }
  return ops;
}

// Structured inputs: paths and stars with monotone weights, in sorted and
// random orders, followed by random edges.
std::vector<std::pair<std::string, Ops>> StructuredOps(int n,
                                                       std::mt19937_64& rng) {
  std::vector<std::pair<std::string, Ops>> cases;
  std::vector<int> perm(n);
  for (int i = 0; i < n; i++) perm[i] = i;
  std::shuffle(perm.begin(), perm.end(), rng);
  auto with_queries = [&](Ops ops) {
    Ops out;
    std::uniform_int_distribution<int> vertex(0, n - 1);
    for (auto& op : ops) {
      out.push_back(op);
      out.emplace_back(vertex(rng), vertex(rng), -1);
    }
    auto tail = RandomOps(n, 2 * n, 1000, 0.3, rng);
    out.insert(out.end(), tail.begin(), tail.end());
    return out;
  };
  Ops inc, dec, star, rnd;
  for (int i = 0; i + 1 < n; i++) {
    inc.emplace_back(perm[i], perm[i + 1], i + 1);
    dec.emplace_back(perm[i], perm[i + 1], n - i);
    star.emplace_back(perm[0], perm[i + 1], i % 3 == 0 ? n - i : i + 1);
  }
  rnd = inc;
  std::shuffle(rnd.begin(), rnd.end(), rng);
  Ops rev(inc.rbegin(), inc.rend());
  cases.emplace_back("path-increasing", with_queries(inc));
  cases.emplace_back("path-decreasing", with_queries(dec));
  cases.emplace_back("path-reversed", with_queries(rev));
  cases.emplace_back("path-shuffled", with_queries(rnd));
  cases.emplace_back("star", with_queries(star));
  return cases;
}

template <class T>
void Test(uint64_t seed) {
  std::mt19937_64 rng(seed);
  int cases = 0;
  int before = failures;
  for (int n : {1, 2, 3, 5, 8, 20, 64, 300}) {
    for (int rep = 0; rep < (n <= 64 ? 30 : 5); rep++) {
      int m = 4 * n + 10;
      // Small weight ranges force many equal weights.
      for (int64_t max_w : {int64_t{3}, int64_t{n}, int64_t{1'000'000'000}}) {
        for (double qr : {0.0, 0.5}) {
          cases++;
          if (!Run<T>("random", n, RandomOps(n, m, max_w, qr, rng))) return;
        }
      }
    }
    for (auto& [label, ops] : StructuredOps(n, rng)) {
      cases++;
      if (!Run<T>(label.c_str(), n, ops)) return;
    }
  }
  if (failures == before) {
    std::printf("ok   %-24s %d cases\n", T::Name().c_str(), cases);
  }
}

}  // namespace

int main(int argc, char** argv) {
  std::setvbuf(stdout, nullptr, _IOLBF, 0);  // show progress if it hangs
  uint64_t seed = argc > 1 ? std::strtoull(argv[1], nullptr, 10) : 1;
  using tmst::LinkMode;
  Test<tmst::OECForest<>>(seed);
  Test<tmst::WeightedDSU<>>(seed);
  Test<tmst::LazyAMTree<LinkMode::kPerch>>(seed);
  Test<tmst::LazyAMTree<LinkMode::kStitch>>(seed);
  Test<tmst::StrictAMTree<LinkMode::kPerch>>(seed);
  Test<tmst::StrictAMTree<LinkMode::kStitch>>(seed);
  Test<tmst::StrictAMTree<LinkMode::kPerch, int, true>>(seed);
  Test<tmst::StrictAMTree<LinkMode::kStitch, int, true>>(seed);
  Test<tmst::BottleneckTree<>>(seed);
  Test<tmst::RivalTree<LinkMode::kPerch>>(seed);
  Test<tmst::RivalTree<LinkMode::kStitch>>(seed);
  Test<tmst::RivalTree<LinkMode::kPerch, int, true>>(seed);
  Test<tmst::RivalTree<LinkMode::kStitch, int, true>>(seed);
  Test<tmst::LinkCutTree<>>(seed);
  // Other weight types.
  Test<tmst::LazyAMTree<LinkMode::kStitch, int64_t>>(seed);
  Test<tmst::StrictAMTree<LinkMode::kPerch, double>>(seed);
  Test<tmst::RivalTree<LinkMode::kStitch, int64_t>>(seed);
  if (failures > 0) {
    std::printf("%d failure(s)\n", failures);
    return 1;
  }
  return 0;
}

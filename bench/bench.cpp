// Times every structure on the same random workload.
//
// Usage: tmst_bench [n] [m] [q] [seed]
//   n vertices, m random edge insertions with random weights, then q random
//   path-max queries.

#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <vector>

#include "tmst/tmst.h"

namespace {

struct Workload {
  int n;
  struct Edge {
    int u, v;
    int w;
  };
  std::vector<Edge> edges;
  std::vector<std::pair<int, int>> queries;
};

Workload MakeRandomWorkload(int n, int m, int q, uint64_t seed) {
  std::mt19937_64 rng(seed);
  std::uniform_int_distribution<int> vertex(0, n - 1);
  std::uniform_int_distribution<int> weight(0, 1'000'000'000);
  Workload wl{n, {}, {}};
  wl.edges.reserve(m);
  for (int i = 0; i < m; i++) {
    wl.edges.push_back({vertex(rng), vertex(rng), weight(rng)});
  }
  wl.queries.reserve(q);
  for (int i = 0; i < q; i++) wl.queries.push_back({vertex(rng), vertex(rng)});
  return wl;
}

template <tmst::IncrementalMST T>
void Run(const Workload& wl) {
  using Clock = std::chrono::steady_clock;
  T t(wl.n);
  auto t0 = Clock::now();
  for (const auto& e : wl.edges) t.Insert(e.u, e.v, e.w);
  auto t1 = Clock::now();
  uint64_t checksum = 0;
  for (auto [u, v] : wl.queries) {
    auto r = t.PathMax(u, v).value_or(-1);
    checksum = checksum * 31 + static_cast<uint64_t>(r);
  }
  auto t2 = Clock::now();
  auto sec = [](auto d) { return std::chrono::duration<double>(d).count(); };
  std::printf("%-28s insert %8.3fs  query %8.3fs  checksum %016llx\n",
              T::Name().c_str(), sec(t1 - t0), sec(t2 - t1),
              static_cast<unsigned long long>(checksum));
}

template <class... Ts>
void RunAll(const Workload& wl) {
  (Run<Ts>(wl), ...);
}

}  // namespace

int main(int argc, char** argv) {
  int n = argc > 1 ? std::atoi(argv[1]) : 1'000'000;
  int m = argc > 2 ? std::atoi(argv[2]) : 10'000'000;
  int q = argc > 3 ? std::atoi(argv[3]) : 10'000'000;
  uint64_t seed = argc > 4 ? std::strtoull(argv[4], nullptr, 10) : 0;
  std::printf("n = %d, m = %d, q = %d, seed = %llu\n", n, m, q,
              static_cast<unsigned long long>(seed));

  auto wl = MakeRandomWorkload(n, m, q, seed);

  using tmst::LinkMode;
  RunAll<tmst::OECForest<>, tmst::WeightedDSU<>,
         tmst::LazyAMTree<LinkMode::kPerch>,
         tmst::LazyAMTree<LinkMode::kStitch>,
         tmst::StrictAMTree<LinkMode::kPerch>,
         tmst::StrictAMTree<LinkMode::kStitch>,
         tmst::StrictAMTree<LinkMode::kPerch, int, true>,
         tmst::StrictAMTree<LinkMode::kStitch, int, true>,
         tmst::BottleneckTree<>,
         tmst::RivalTree<LinkMode::kPerch>, tmst::RivalTree<LinkMode::kStitch>,
         tmst::RivalTree<LinkMode::kPerch, int, true>,
         tmst::RivalTree<LinkMode::kStitch, int, true>,
         tmst::LinkCutTree<>>(wl);
  return 0;
}

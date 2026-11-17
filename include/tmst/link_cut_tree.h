#ifndef TMST_LINK_CUT_TREE_H_
#define TMST_LINK_CUT_TREE_H_

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

// Baseline: a link-cut tree [Sleator and Tarjan, 1983] on the MST, with the
// edges as extra nodes. O(log n) amortized per operation.
template <class W = int>
class LinkCutTree {
 public:
  using weight_type = W;

  explicit LinkCutTree(int n) : n(n) {
    const std::size_t max_edges = n == 0 ? 0 : n - 1;
    nodes.reserve(static_cast<std::size_t>(n) + max_edges);
    nodes.resize(n);
    endpoints.reserve(max_edges);
  }

  static std::string Name() { return "link-cut-tree"; }

  void Insert(int u, int v, W w) {
    assert(w < std::numeric_limits<W>::max());
    if (u == v) return;
    NodeId edge = ExposePath(u, v);
    if (edge == kNone) {
      edge = static_cast<NodeId>(nodes.size());
      nodes.emplace_back();
      endpoints.emplace_back(u, v);
    } else {
      if (!(w < nodes[edge].weight)) return;
      auto [a, b] = endpoints[edge - n];
      Cut(edge, a);
      Cut(edge, b);
      endpoints[edge - n] = {u, v};
    }

    // Both components can now be linked to the isolated edge node.
    nodes[edge].weight = w;
    nodes[edge].max_edge = edge;
    MakeRoot(u);
    MakeRoot(v);
    nodes[u].parent = edge;
    nodes[v].parent = edge;
  }

  std::optional<W> PathMax(int u, int v) {
    if (u == v) return std::nullopt;
    NodeId edge = ExposePath(u, v);
    if (edge == kNone) return std::nullopt;
    return nodes[edge].weight;
  }

  bool Check() const {
    if (nodes.size() != n + endpoints.size()) return false;
    for (auto [u, v] : endpoints) {
      if (u >= n || v >= n || u == v) return false;
    }
    for (NodeId x = 0; x < nodes.size(); x++) {
      const auto& node = nodes[x];
      if (node.parent != kNone && node.parent >= nodes.size()) return false;
      if (node.max_edge != kNone &&
          (node.max_edge < n || node.max_edge >= nodes.size())) {
        return false;
      }
      for (NodeId child : node.child) {
        if (child != kNone &&
            (child >= nodes.size() || nodes[child].parent != x)) {
          return false;
        }
      }
      if (node.child[0] != kNone && node.child[0] == node.child[1]) {
        return false;
      }
    }

    // The auxiliary parents and the path parents together form a forest.
    std::vector<uint8_t> state(nodes.size(), 0);
    for (NodeId x = 0; x < nodes.size(); x++) {
      if (nodes[x].max_edge != SubtreeMax(x)) return false;
      NodeId y = x;
      while (y != kNone && state[y] == 0) {
        state[y] = 1;
        y = nodes[y].parent;
      }
      if (y != kNone && state[y] == 1) return false;
      for (y = x; y != kNone && state[y] == 1; y = nodes[y].parent) {
        state[y] = 2;
      }
    }
    return true;
  }

 private:
  using NodeId = uint32_t;
  static constexpr NodeId kNone = std::numeric_limits<NodeId>::max();

  struct Node {
    W weight{};
    // Parent in the auxiliary splay tree, or path parent at an auxiliary root.
    NodeId parent = kNone;
    NodeId child[2] = {kNone, kNone};
    NodeId max_edge = kNone;
    bool reversed = false;
  };

  NodeId Heavier(NodeId a, NodeId b) const {
    if (a == kNone) return b;
    if (b == kNone) return a;
    // Break ties by node ID so reversing a path leaves its maximum unchanged.
    if (b < a) std::swap(a, b);
    return nodes[a].weight < nodes[b].weight ? b : a;
  }

  NodeId SubtreeMax(NodeId x) const {
    NodeId best = x < n ? kNone : x;
    for (NodeId child : nodes[x].child) {
      if (child != kNone) best = Heavier(best, nodes[child].max_edge);
    }
    return best;
  }

  void Pull(NodeId x) { nodes[x].max_edge = SubtreeMax(x); }

  bool IsAuxRoot(NodeId x) const {
    NodeId p = nodes[x].parent;
    return p == kNone || (nodes[p].child[0] != x && nodes[p].child[1] != x);
  }

  void Reverse(NodeId x) {
    auto& node = nodes[x];
    std::swap(node.child[0], node.child[1]);
    node.reversed = !node.reversed;
  }

  void Push(NodeId x) {
    if (!nodes[x].reversed) return;
    for (NodeId child : nodes[x].child) {
      if (child != kNone) Reverse(child);
    }
    nodes[x].reversed = false;
  }

  void Rotate(NodeId x) {
    NodeId y = nodes[x].parent;
    NodeId z = nodes[y].parent;
    int side = nodes[y].child[1] == x;
    NodeId middle = nodes[x].child[1 - side];
    if (!IsAuxRoot(y)) nodes[z].child[nodes[z].child[1] == y] = x;
    nodes[x].parent = z;
    nodes[y].parent = x;
    nodes[x].child[1 - side] = y;
    nodes[y].child[side] = middle;
    if (middle != kNone) nodes[middle].parent = y;
    Pull(y);
    Pull(x);
  }

  void Splay(NodeId x) {
    ancestors.clear();
    ancestors.push_back(x);
    for (NodeId y = x; !IsAuxRoot(y);) {
      y = nodes[y].parent;
      ancestors.push_back(y);
    }
    for (auto it = ancestors.rbegin(); it != ancestors.rend(); ++it) Push(*it);

    while (!IsAuxRoot(x)) {
      NodeId y = nodes[x].parent;
      if (!IsAuxRoot(y)) {
        NodeId z = nodes[y].parent;
        bool same_side = (nodes[y].child[1] == x) == (nodes[z].child[1] == y);
        Rotate(same_side ? y : x);
      }
      Rotate(x);
    }
  }

  void Access(NodeId x) {
    NodeId preferred = kNone;
    for (NodeId y = x; y != kNone; y = nodes[y].parent) {
      Splay(y);
      // The old and new right children both point to y, as a path parent
      // or an auxiliary parent depending on which child is preferred.
      nodes[y].child[1] = preferred;
      Pull(y);
      preferred = y;
    }
    Splay(x);
  }

  void MakeRoot(NodeId x) {
    Access(x);
    Reverse(x);
  }

  NodeId ExposePath(NodeId u, NodeId v) {
    MakeRoot(u);
    Access(v);
    // For distinct connected vertices, u is now a proper descendant of v
    // in the auxiliary tree. Otherwise u stays an auxiliary root.
    return nodes[u].parent == kNone ? kNone : nodes[v].max_edge;
  }

  void Cut(NodeId x, NodeId y) {
    MakeRoot(x);
    Access(y);
    assert(nodes[y].child[0] == x && nodes[x].child[1] == kNone);
    nodes[y].child[0] = kNone;
    nodes[x].parent = kNone;
    Pull(y);
  }

  NodeId n;
  std::vector<Node> nodes;
  std::vector<std::pair<NodeId, NodeId>> endpoints;
  std::vector<NodeId> ancestors;
};

}  // namespace tmst

#endif  // TMST_LINK_CUT_TREE_H_

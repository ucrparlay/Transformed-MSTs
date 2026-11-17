# OEC-forest and Rival-tree technical report

This report accompanies *The Family of Balanced T-MSTs for Incremental Minimum Spanning Trees*. The [OEC-forest analysis](#sec-oec) follows Section 6 of the paper and covers the exact depth law and its extensions to priority-based stitch and weighted DSU. The [Rival-tree analysis](#sec-hb) follows Section 5 and covers the balance rule, both link procedures, and the worst-case stitch construction. Algorithm, theorem, and lemma numbers are local to each analysis.

## Preliminaries

Following Section 2 of the paper, let $G_t$ be the undirected graph on $n$ vertices after the first $t$ insertions, with real edge weights. For distinct vertices $u,v$ in the same component, the *path-max value* is

```math
\mathrm{PM}_{G_t}(u,v)=\min_{P:u\leadsto v}\ \max_{e\in P} w(e).
```

We set $\mathrm{PM}_{G_t}(u,v)=+\infty$ for vertices in different components and $\mathrm{PM}_{G_t}(u,u)=-\infty$ for the empty path. The C++ method `PathMax` returns `std::nullopt` in both cases.

The maintained forest $T_t$ stores a parent $p(v)$ and the weight $\mathit{wt}(v)$ of the parent edge. At a root, $p(v)=\bot$ (`null` in the pseudocode) and $\mathit{wt}(v)=+\infty$. The C++ implementations require input weights strictly below `std::numeric_limits<W>::max()`, which represents the root weight. Depth counts the edges from a vertex to its root, and height counts the edges on its longest downward path.

A transformed minimum spanning tree (T-MST) of $G_t$ is a rooted forest with the same path-max values as $G_t$ and the same multiset of edge weights as a minimum spanning forest of $G_t$. For a forest, path-max equivalence already implies the weight condition. A minimum spanning forest realizes the path-max values of $G_t$ on its tree paths. Thus, at every threshold $\tau$, both forests have the same components using edges of weight at most $\tau$, and in each forest the number of such edges is $n$ minus the number of these components.

<a name="tw-transformations"></a>

A *TW transformation* on adjacent edges $(a,b,w_1)$ and $(b,c,w_2)$ with $w_1\le w_2$ replaces the latter edge by $(a,c,w_2)$. At every threshold that includes the moved edge, the lighter edge already joins $a$ and $b$. The threshold components, and hence the path-max values, are preserved. Deleting a self-loop or a maximum-weight edge of a cycle also preserves them.

## The OEC-Forest Has Logarithmic Height

<a name="sec-oec"></a>

### Model and algorithm

<a name="sec-oec-model"></a>

**Random priorities.**

Each vertex $v$ receives a priority $\pi(v)$, and $\pi$ is a uniformly random bijection from $V$ to $`\{1,\dots,n\}`$. We say that $u$ is *higher* than $v$ if $\pi(u)>\pi(v)$, and *lower* otherwise. A rooted forest is *heap-ordered* if $\pi(p(v))>\pi(v)$ for every non-root $v$. For a set $S\subseteq V$, we write $\max_\pi S$ for the highest vertex in $S$.

**The input model.**

The insertion sequence $\sigma=(e_1,\dots,e_M)$ is fixed in advance, independently of the random priorities $\pi$. Path-max queries may be interleaved with insertions and are fixed in advance as well. Queries leave the OEC-forest and priority-based stitch unchanged. Weighted DSU shortcuts pointers during queries, as handled by [Lemma 7](#lem-wdsu). The OEC-forest allows arbitrary weights, including ties. The analysis of [priority-based stitch](#sec-two-sided) and weighted DSU assumes distinct edge weights.

The depth distributions and height probabilities in this section are taken over the uniformly random permutation $\pi$ for the fixed input. The C++ constructors for `OECForest` and `WeightedDSU` shuffle the vertices with `std::mt19937_64`. Both default to seed `0`, so the same $n$ and seed reproduce the same priority order in a given build. [OECForest implementation](include/tmst/oec_forest.h) · [WeightedDSU implementation](include/tmst/weighted_dsu.h)

**The OEC-forest.**

The OEC-forest \[[Song et al., 2024](#ref-song2024oec); [Xu et al., 2026](#ref-xu2026oecjournal)\] maintains a heap-ordered rooted forest $T_t$ that is path-max equivalent to $G_t$. Initially every vertex is a root. To insert $e_t=(u_t,v_t,w_t)$, the algorithm calls $\mathrm{Update}(u_t,v_t,w_t)$, shown in [Algorithm 1](#alg-update). The original description uses a maximum spanning forest over arrival times; [Algorithm 1](#alg-update) is the same procedure stated for minimum spanning forests and arbitrary weights (negating all weights translates one into the other).

<a name="alg-update"></a>

**Algorithm 1. The OEC-forest update.**

```
Update(u, v, w):
    if u == v or u == null or v == null: return
    if priority(u) < priority(v): swap(u, v)
    # v is the lower endpoint
    if w <= wt(v):
        Update(u, p(v), wt(v))
        p(v) = u; wt(v) = w
    else:
        Update(u, p(v), w)
```

Each step applies a TW transformation at the pending edge's lower endpoint $v$. When $w\le\mathit{wt}(v)$, attach $v$ below $u$ and slide its old parent edge to $(u,p(v),\mathit{wt}(v))$ using the old fields. Otherwise slide the pending edge along the parent edge of $v$. The update stops when the pending edge becomes a self-loop or a root is attached. Since $v$ is always attached below the higher endpoint, the forest stays heap-ordered. Path-max equivalence gives the same connected components as $G_t$.

**Fact 1 (\[[Song et al., 2024](#ref-song2024oec)\]).**

<a name="fact-correct"></a>

For every $t$, the forest $T_t$ is heap-ordered and path-max equivalent to $G_t$.

**Fact 2.**

<a name="fact-heap"></a>

For every $t$, the trees of $T_t$ are exactly the connected components of $G_t$. In particular, the root of the tree containing $x$ is $\max_\pi C$, where $C$ is the connected component of $x$ in $G_t$.

**Fact 3.**

<a name="fact-cost"></a>

The call $\mathrm{Update}(u,v,w)$ makes at most $\mathrm{depth}(u)+\mathrm{depth}(v)+2$ recursive calls, and a path-max query on $u$ and $v$ takes $O(1+\mathrm{depth}(u)+\mathrm{depth}(v))$ time, where depths are measured in the current forest.

**Proof.**

Every call of [Algorithm 1](#alg-update) makes its recursive call before writing any field, so all fields read during an update are those of the forest before the update. Each nonterminal call moves one endpoint one step up its original root path or past its root.

A query advances the lower-priority endpoint and records the largest edge weight traversed. By heap order, a proper ancestor has higher priority than its descendant, so on connected inputs the pointers meet at their lowest common ancestor. Reaching a root at the lower-priority endpoint while the pointers differ certifies disconnection. Each advance costs constant time and traverses one parent edge, proving the query bound. $\square$

Hence it suffices to bound the depth of every vertex at every time. Throughout, $H_m=\sum_{i=1}^{m}1/i$ denotes the $m$-th harmonic number, with $H_0=0$.

### The depth of the OEC-forest

<a name="sec-proof"></a>

Throughout this subsection, fix the insertion sequence $\sigma$ and a time $t$, and write $T=T_t$, $p=p_t$ and $\mathit{wt}=\mathit{wt}_t$. All probabilities are over the random priorities $\pi$.

#### Revealing the priorities from below

<a name="sec-reveal"></a>

Let $\ell_1,\ell_2,\dots,\ell_n$ be the vertices listed in increasing order of priority, so $\pi(\ell_i)=i$. For $`k\in\{0,\dots,n\}`$, the *prefix* of length $k$ is the random sequence $\Lambda_k=(\ell_1,\dots,\ell_k)$.

Revealing $\Lambda_k$ gives the relative order of the $k$ lowest vertices and leaves the order of the suffix hidden.

**Lemma 1 (Uniform suffix).**

<a name="lem-suffix"></a>

For every $k$ and every sequence $\lambda$ of $k$ distinct vertices with $\Pr[\Lambda_k=\lambda]>0$, conditioned on $\Lambda_k=\lambda$, the sequence $(\ell_{k+1},\dots,\ell_n)$ is a uniformly random ordering of $V\setminus\lambda$.

The next lemma is the heart of the analysis. It states that the parents of the revealed vertices are already determined by the prefix, even when those parents are unrevealed.

**Lemma 2 (Determinacy).**

<a name="lem-determinacy"></a>

For every $k$, the values $p(v)$ and $\mathit{wt}(v)$ for all $v\in \Lambda_k$ are determined by $\sigma$ and $\Lambda_k$. That is, there are functions $P_\sigma$ and $W_\sigma$ such that $p(v)=P_\sigma(\Lambda_k,v)$ and $\mathit{wt}(v)=W_\sigma(\Lambda_k,v)$ for all $v\in \Lambda_k$.

**Proof.**

Consider the execution of [Algorithm 1](#alg-update) on $\sigma$, viewed as a sequence of calls $\mathrm{Update}(a,b,\omega)$ in execution order. In each call that does not return immediately, let $v$ denote the lower endpoint and $u$ the higher endpoint after the swap. We use three properties of the algorithm.

- **(P1)** A call reads and writes only the two fields $p(v)$ and $\mathit{wt}(v)$ of its lower endpoint $v$. The other endpoint $u$ is used only as a name: it is passed to the child call and possibly written into $p(v)$.

- **(P2)** If both endpoints of a call are outside $\Lambda_k$ (or equal to $\bot$), then neither this call nor any of its descendant calls reads or writes a field of a vertex in $\Lambda_k$. Indeed, the child call has endpoints $u$ and $p(v)$; the vertex $u$ is outside $\Lambda_k$, and $p(v)$ is either $\bot$ or higher than $v$ by heap order, hence also outside $\Lambda_k$. By induction, every descendant call has both endpoints outside $\Lambda_k$, and by (P1) it only touches fields of vertices outside $\Lambda_k$.

- **(P3)** If at least one endpoint of a call lies in $\Lambda_k$, then the lower endpoint lies in $\Lambda_k$, and $\Lambda_k$ determines which endpoint is the lower one: if both endpoints lie in $\Lambda_k$ this is given by their order in $\Lambda_k$, and otherwise the lower endpoint is the unique one in $\Lambda_k$.

Now run [Algorithm 1](#alg-update) on $\sigma$, storing fields only for the vertices of $\Lambda_k$ and skipping every call whose two endpoints both lie outside $\Lambda_k$. By (P2), the skipped calls and their descendants do not touch the fields of $\Lambda_k$. Every other call can identify its lower endpoint $v\in\Lambda_k$ by (P3), and by (P1) it computes its writes and its child call from its arguments and the fields of $v$ alone. By induction over the execution, this procedure, which uses only $\sigma$ and $\Lambda_k$, makes the same calls with an endpoint in $\Lambda_k$ as the true execution and writes the same values, which proves the lemma. $\square$

#### An exact recursion for the depth

<a name="sec-recursion"></a>

<a name="eq-recursion"></a>

Define random variables $D_0,D_1,D_2,\dots$ by

```math
D_0 = 0,
  \qquad
  D_m \;\stackrel{d}{=}\; 1 + D_{J_m}\quad (m\ge 1),
\tag{1}
```
where $J_m$ is uniformly distributed on $`\{0,1,\dots,m-1\}`$ and independent of $D_0,\dots,D_{m-1}$. The number of records (left-to-right maxima) of a uniform permutation satisfies this recursion. The maximum is the final record, and the elements preceding it form a uniformly ordered prefix of length $J_m$. [Knuth, *The Art of Computer Programming*, Vol. 1, 3rd ed. (1997), Section 1.2.10](#ref-knuth1997taocp1)

Throughout, $B_1,B_2,\dots$ are independent Bernoulli random variables with $\mathbb{E}[B_i]=1/i$, so $B_1=1$.

**Lemma 3 (Law of $D_m$).**

<a name="lem-moments"></a>

For every $m\ge 0$, the random variable $D_m$ has the same distribution as $B_1+\dots+B_m$. In particular, $\mathbb{E}[D_m]=H_m$, where $H_m$ is the $m$-th harmonic number, and $\mathbb{E}[a^{D_m}]=\prod_{i=1}^{m}\bigl(1+\frac{a-1}{i}\bigr)$ for every $a>0$; for $a=2$ this equals $m+1$.

**Proof.**

Since $D_m$ takes values in $`\{0,\dots,m\}`$, its distribution is determined by its generating polynomial $g_m(z)=\mathbb{E}[z^{D_m}]$. By [(1)](#eq-recursion), $g_0(z)=1$ and $g_m(z)=\frac{z}{m}\sum_{j=0}^{m-1}g_j(z)$ for $m\ge1$. The generating polynomial of $B_1+\dots+B_m$ is $P_m(z)=\prod_{i=1}^{m}\frac{z+i-1}{i}$, because $\mathbb{E}[z^{B_i}]=1-\frac1i+\frac zi$. It satisfies $jP_j(z)=(z+j-1)P_{j-1}(z)$, that is, $zP_{j-1}(z)=jP_j(z)-(j-1)P_{j-1}(z)$, for every $j\ge 1$. Summing over $j=1,\dots,m$, the right-hand side telescopes and gives $z\sum_{j=0}^{m-1}P_j(z)=mP_m(z)$. Hence the polynomials $P_m$ satisfy the same recursion as $g_m$, with $P_0=1$, so $g_m=P_m$ for every $m$ by induction. The moments follow from independence: $\mathbb{E}[D_m]=\sum_{i=1}^m\frac1i=H_m$ and $\mathbb{E}[a^{D_m}]=\prod_{i=1}^m\mathbb{E}[a^{B_i}]=\prod_{i=1}^m\bigl(1+\frac{a-1}{i}\bigr)$, and for $a=2$ the product $\prod_{i=1}^m\frac{i+1}{i}$ telescopes to $m+1$. $\square$

For a vertex $z$, let $C(z)$ denote the connected component of $z$ in $G_t$; it depends only on $\sigma$ and not on $\pi$. The next lemma shows that, conditioned on the prefix up to $z$, the depth of $z$ follows the recursion [(1)](#eq-recursion).

**Lemma 4 (Law of the depth).**

<a name="lem-law"></a>

Let $k\ge 1$, and let $\lambda=(\lambda_1,\dots,\lambda_k)$ be a sequence of distinct vertices with $\Pr[\Lambda_k=\lambda]>0$. Let $z=\lambda_k$ and $m=|C(z)\setminus\lambda|$, the number of vertices of $C(z)$ that are higher than $z$. Then, conditioned on $\Lambda_k=\lambda$, the depth $\mathrm{depth}_T(z)$ has the same distribution as $D_m$.

**Proof.**

We use induction on $m$. If $m=0$, then $z=\max_\pi C(z)$, so $z$ is a root by [Fact 2](#fact-heap) and its depth is $0=D_0$.

Let $m\ge 1$. Then $z$ is not the highest vertex of $C(z)$, so it is not a root, and its parent $y=p(z)$ is a vertex of $C(z)$ that is higher than $z$, that is, $y\in C(z)\setminus\lambda$. By [Lemma 2](#lem-determinacy), conditioned on $\Lambda_k=\lambda$, the vertex $y=P_\sigma(\lambda,z)$ is a fixed vertex. By [Lemma 1](#lem-suffix), the vertices of $V\setminus\lambda$ appear in uniformly random order in the rest of the permutation. Hence the number $J$ of vertices of $C(z)$ that are higher than $y$ is uniformly distributed on $`\{0,1,\dots,m-1\}`$.

To apply the induction hypothesis to $y$, we refine the conditioning to the prefix that ends with $y$. Let $\mathcal{E}$ be the set of sequences $\lambda'$ of distinct vertices that extend $\lambda$ and end with $y$; each of them has $\Pr[\Lambda_{|\lambda'|}=\lambda']>0$. For $\lambda'\in\mathcal{E}$, let $A_{\lambda'}$ denote the event $\Lambda_{|\lambda'|}=\lambda'$. Since $y$ is a fixed vertex outside $\lambda$, conditioned on $\Lambda_k=\lambda$, exactly one of the events $A_{\lambda'}$ with $\lambda'\in\mathcal{E}$ occurs, namely the one with $|\lambda'|$ equal to the position of $y$. Moreover, each $A_{\lambda'}$ implies $\Lambda_k=\lambda$. Let $j(\lambda')=|C(y)\setminus\lambda'|$; since $C(y)=C(z)$, this is the value of $J$ on $A_{\lambda'}$, and $j(\lambda')<m$. By the induction hypothesis applied to $y$ and $\lambda'$, conditioned on $A_{\lambda'}$, the depth of $y$ is distributed as $D_{j(\lambda')}$. As $\mathrm{depth}_T(z)=1+\mathrm{depth}_T(y)$, the law of total probability gives, for every integer $d$,
```math
\begin{aligned}
  \Pr\bigl[\mathrm{depth}_T(z)=d \mid \Lambda_k=\lambda\bigr]
  &= \sum_{\lambda'\in\mathcal{E}}\Pr\bigl[A_{\lambda'} \mid \Lambda_k=\lambda\bigr]\cdot\Pr\bigl[\mathrm{depth}_T(y)=d-1 \mid A_{\lambda'}\bigr]\\
  &= \sum_{\lambda'\in\mathcal{E}}\Pr\bigl[A_{\lambda'} \mid \Lambda_k=\lambda\bigr]\cdot\Pr\bigl[D_{j(\lambda')}=d-1\bigr]\\
  &= \sum_{j=0}^{m-1}\Pr\bigl[J=j \mid \Lambda_k=\lambda\bigr]\cdot\Pr[D_j=d-1]\\
  &= \frac1m\sum_{j=0}^{m-1}\Pr[D_j=d-1]
  = \Pr[1+D_{J_m}=d],
\end{aligned}
```
where the third equality groups the sequences $\lambda'$ according to the value $j(\lambda')$. By [(1)](#eq-recursion), this is $\Pr[D_m=d]$. $\square$

#### Main results

<a name="sec-main"></a>

**Theorem 1 (Exact law of the depth).**

<a name="thm-oec-exact"></a>

Fix any insertion sequence $\sigma$, any time $t$, and any vertex $x$, and let $s=|C(x)|$. Then $`\mathrm{depth}_{T_t}(x)`$ has the same distribution as $D_s-1$, that is, as $B_2+B_3+\dots+B_s$. In particular, $`\mathbb{E}[\mathrm{depth}_{T_t}(x)] = H_s-1`$, $`\mathbb{E}\bigl[2^{\mathrm{depth}_{T_t}(x)}\bigr] = (s+1)/2`$, and $`\mathbb{E}\bigl[a^{\mathrm{depth}_{T_t}(x)}\bigr]=\prod_{i=2}^{s}\bigl(1+\frac{a-1}{i}\bigr)`$ for every $a>0$.

**Proof.**

Let $k=\pi(x)$ be the position of $x$, so that $\Lambda_k$ is the prefix ending with $x$. The number $R$ of vertices of $C(x)$ higher than $x$ is determined by $\Lambda_k$, and since the rank of $x$ within $C(x)$ is uniform, $R$ is uniformly distributed on $`\{0,\dots,s-1\}`$. By [Lemma 4](#lem-law) applied to every possible value of $\Lambda_k$, the depth of $x$ conditioned on $R=m$ is distributed as $D_m$. Hence, for every $d$,
```math
\Pr[\mathrm{depth}(x)=d] = \frac1s\sum_{m=0}^{s-1}\Pr[D_m=d] = \Pr[D_{J_s}=d] = \Pr[D_s-1=d],
```
where the last equality is the recursion [(1)](#eq-recursion) for $D_s$. By [Lemma 3](#lem-moments) and $B_1=1$, the variable $D_s-1$ has the same distribution as $B_2+\dots+B_s$, and the moments follow as in [Lemma 3](#lem-moments); for $a=2$ the product $\prod_{i=2}^{s}\frac{i+1}{i}$ telescopes to $(s+1)/2$. $\square$

The exact law concerns one vertex at a time. The joint distribution of the depths and the shape of $T_t$ can depend on the insertion history.

**Theorem 2 (Height).**

<a name="thm-oec-tail"></a>

Fix any insertion sequence $\sigma$ and any time $t$.

- **(i)** For every $h\ge 0$, $\Pr\bigl[\mathrm{height}(T_t)\ge h\bigr] \le n(n+1)/2^{h+1}$. In particular, for every $c>0$, $\mathrm{height}(T_t)<(c+2)\log_2 n$ with probability at least $1-n^{-c}$.

- **(ii)** For every $\beta\ge 1$, $\Pr\bigl[\mathrm{height}(T_t)\ge \beta\ln n\bigr] \le n^{-\beta(\ln\beta-1)}$. In particular, for every constant $\beta>e$, $\mathrm{height}(T_t)<\beta\ln n$ with probability $1-n^{-\Omega(1)}$.

**Proof.**

Let $x$ be a vertex, $s=|C(x)|\le n$, and $a\ge 1$. By Markov’s inequality applied to $a^{\mathrm{depth}(x)}$, we have $\Pr[\mathrm{depth}(x)\ge h]\le\mathbb{E}[a^{\mathrm{depth}(x)}]/a^h$ for every $h\ge 0$. For (i), take $a=2$: by [Theorem 1](#thm-oec-exact), $\Pr[\mathrm{depth}(x)\ge h]\le (s+1)/2^{h+1}\le(n+1)/2^{h+1}$. A union bound over the $n$ vertices gives the first claim, and for $h=(c+2)\log_2 n$ the bound is $n(n+1)/(2n^{c+2})\le n^{-c}$. For (ii), [Theorem 1](#thm-oec-exact), the inequality $1+y\le e^y$, and $H_s-1\le\ln s\le\ln n$ give
```math
\mathbb{E}\bigl[a^{\mathrm{depth}(x)}\bigr]=\prod_{i=2}^{s}\Bigl(1+\frac{a-1}{i}\Bigr)\le e^{(a-1)(H_s-1)}\le n^{a-1}.
```
Take $a=\beta$ and $h=\beta\ln n$, so that $a^h=n^{\beta\ln\beta}$. Then $\Pr[\mathrm{depth}(x)\ge\beta\ln n]\le n^{\beta-1-\beta\ln\beta}$, and a union bound over the $n$ vertices gives $n^{\beta-\beta\ln\beta}=n^{-\beta(\ln\beta-1)}$. For $\beta>e$ the exponent $\beta(\ln\beta-1)$ is a positive constant. $\square$

For example, part (ii) gives height less than $3.6\ln n\approx2.5\log_2 n$ with probability at least $1-1/n$, compared with $3\log_2 n$ from part (i). When $G_t$ is connected, a fixed vertex has expected depth $H_n-1$. Applying a Chernoff lower-tail bound to its Bernoulli sum also gives depth $\Omega(\log n)$ with high probability.

**Corollary 1.**

<a name="cor-ops"></a>

For an insertion sequence of length $M$, with probability at least $`1-M\,n^{-c}`$, the forest has height less than $(c+2)\log_2 n$ at all times. In particular, if $M$ is polynomial in $n$, then every update and every path-max query takes $O(\log n)$ time with high probability.

**Proof.**

Apply [Theorem 2](#thm-oec-tail) at each of the $M$ times and take a union bound; then use [Fact 3](#fact-cost). $\square$

### The priority-based stitch

<a name="sec-two-sided"></a>

The *priority-based stitch* uses [Algorithm 2](#alg-priority-stitch) with the rule that attaches the endpoint of lower priority. Either endpoint can move up past lighter parent edges. A revealed prefix therefore determines a lower bound on its vertices' parent priorities rather than their exact parents. This is the stitch-based AM-tree with priorities in place of subtree sizes. Weighted DSU adds the steps in [Remark 1](#rem-wdsu). [Ding et al. (2025)](#ref-ding2025amtree) · [not_amir (2024)](#ref-notamir2024wdsu)

The two endpoints are called *pointers*. A climb moves one pointer to its parent, and attaching $x$ below $y$ moves the pointer at $x$ to its old parent. With the priority rule, each pointer advances to a higher vertex or reaches $\bot$.

<a name="alg-priority-stitch"></a>

**Algorithm 2. The priority-based stitch (distinct weights).**

```
Stitch(u, v, w):
    if u == v or u == null or v == null: return
    if w > wt(u): Stitch(p(u), v, w)
    else if w > wt(v): Stitch(u, p(v), w)
    else:
        x = the endpoint of lower priority
        y = the other endpoint
        old_parent = p(x); old_weight = wt(x)
        p(x) = y; wt(x) = w
        Stitch(old_parent, y, old_weight)
```

<a name="stitch-correctness"></a>

**Stitch correctness.** Assume distinct edge weights and a selection rule that chooses the descendant whenever one endpoint is an ancestor of the other. At the start of each call, the pointers lie on the two original root paths, and each pointer's path to its root still follows the original parent edges. A climb advances one pointer along that path. An attachment changes the parent pointer and edge weight of $x$ alone. The path from $y$ to its root avoids $x$ by the selection rule, so attaching $x$ below $y$ preserves the forest. If $x$ had a parent $x'$, the path from $x'$ to its root also avoids $x$. The next call therefore satisfies the same invariant with pointers $x'$ and $y$. Attaching a root ends the stitch.

Each nonterminal call advances one pointer along its original path or past its root. Including the terminal call, the stitch makes at most $\mathrm{depth}(u)+\mathrm{depth}(v)+2$ calls, with both depths measured before the insertion. Every vertex whose parent or child set changes lies on those two paths. A climb slides the pending edge along a lighter parent edge. An attachment incorporates the pending edge and, if $x$ had a parent, slides its old parent edge along the new one to become the next pending edge. When the pointers meet, the pending edge is a self-loop. These transformations preserve the path-max values of the forest together with the pending edge, so the final forest is a T-MST of the updated graph. [TW transformations](#tw-transformations)

The priority rule preserves heap order at each attachment. In a heap-ordered forest the descendant has lower priority, so the rule meets the selection condition. This proves [Facts 1](#fact-correct), [2](#fact-heap), and [3](#fact-cost) for priority-based stitch.

#### Lower bounds on the parents

Fix $k$ and a prefix $\lambda$ of length $k$, and write $\bar\lambda=V\setminus\lambda$. The *frozen simulation* runs the priority-based stitch on $\sigma$, maintaining fields $\hat p(v)$ and $\widehat{\mathit{wt}}(v)$ only for $v\in \lambda$, with the following changes.

- A pointer located at a vertex of $\bar\lambda$ never moves; that is, the tests in [the first climb](#alg-priority-stitch) and [the second climb](#alg-priority-stitch) are treated as false for such a pointer.

- When both pointers are located at vertices of $\bar\lambda$, the simulation stops processing the current insertion.

Whenever at least one pointer lies in $\lambda$, the lower of the two pointers lies in $\lambda$ and is determined by $\lambda$, as in (P3). Hence the frozen simulation is a function of $\sigma$ and $\lambda$ only.

**Lemma 5 (Simulation).**

<a name="lem-simulation"></a>

For every time $t$ and every $v\in \lambda$:

- **(i)** $\widehat{\mathit{wt}}_t(v)=\mathit{wt}_t(v)$;

- **(ii)** if $`\hat p_t(v)\in \lambda\cup\{\bot\}`$, then $p_t(v)=\hat p_t(v)$;

- **(iii)** if $\hat p_t(v)\in \bar\lambda$, then $p_t(v)\in \bar\lambda$, $\pi(p_t(v))\ge\pi(\hat p_t(v))$, and $\hat p_t(v)$ lies in the connected component of $v$ in $G_t$.

In particular, $p_t(v)\in \lambda$ if and only if $\hat p_t(v)\in \lambda$, and in that case the two coincide.

**Proof.**

Induct over insertions, starting from the isolated vertices. An attachment writes only the fields of $x$, and both pointers then lie above $x$. Since pointers only move up in priority, fields read at their positions still have their pre-insertion values.

Compare corresponding pointers in the two executions. They coincide inside $\lambda$. Outside $\lambda$ the true pointer is at least as high as the simulated pointer. The pending weights agree while either pointer is revealed. This holds at the start and is preserved as follows.

- For a revealed pointer, the two executions compare the same pending weight with the same vertex weight. A climb preserves the comparison of the pointers by (ii) and (iii).
- An unrevealed pointer climbs only in the true execution. Its priority increases and the pending weight stays the same.
- If both pointers are revealed, the executions take the same step. If they meet, both stop.
- At an attachment with one revealed endpoint $x$, the extra climbs of the unrevealed true pointer leave the decisions at $x$ unchanged. Both executions write the same weight at $x$ and attach it below their respective unrevealed pointers. The true parent is at least as high as the simulated one. Both runs then carry the old weight of $x$ and continue from its respective old parent.
- Once both pointers are unrevealed, only unrevealed fields can change.

In the mixed attachment, the simulated pointer started in the component of one endpoint of the inserted edge and followed parent edges within that component, by (ii) and (iii). The inserted edge places both endpoints in one component of $G_t$, so the simulated parent belongs to $C(x)$. Components can only grow after that attachment.

This preserves (i)–(iii). A true parent is revealed exactly when its simulated parent is revealed, and in that case they coincide. $\square$

#### Stochastic domination

For random variables $X$ and $Y$ with values in $\mathbb{Z}_{\ge0}$, we write $X\preceq Y$ if $\Pr[X\ge d]\le\Pr[Y\ge d]$ for all $d$, or equivalently, if $\mathbb{E}[f(X)]\le\mathbb{E}[f(Y)]$ for every non-decreasing function $f$.

**Lemma 6 (Monotonicity).**

<a name="lem-monotone"></a>

$D_0\preceq D_1\preceq D_2\preceq\cdots$.

**Proof.**

Use the Bernoulli representation from [Lemma 3](#lem-moments). Couple all $D_m$ as the partial sums of the same sequence of independent $B_i$. Each additional summand is nonnegative. $\square$

**Theorem 3.**

<a name="thm-two-sided"></a>

For the priority-based stitch, with distinct edge weights and under the assumptions of [Theorem 1](#thm-oec-exact), we have $\mathrm{depth}_{T_t}(x)\preceq D_s-1$. In particular,
```math
\mathbb{E}[\mathrm{depth}_{T_t}(x)] \le H_s-1
  \qquad\text{and}\qquad
  \mathbb{E}\bigl[a^{\mathrm{depth}_{T_t}(x)}\bigr] \le \prod_{i=2}^{s}\Bigl(1+\frac{a-1}{i}\Bigr)
```
for every $a\ge1$; for $a=2$ the right-hand side is $(s+1)/2$. Consequently, [Theorem 2](#thm-oec-tail) and [Corollary 1](#cor-ops) hold for this variant as well.

**Proof.**

We show the analogue of [Lemma 4](#lem-law) with equality in distribution replaced by $\preceq$: for every prefix $\lambda$ of length $k\ge 1$ ending with $z$, with $m=|C(z)\setminus\lambda|$, conditioned on $\Lambda_k=\lambda$ we have $\mathrm{depth}(z)\preceq D_m$. The theorem then follows as in the proof of [Theorem 1](#thm-oec-exact): with $R$ as there, for every $d$,
```math
\Pr[\mathrm{depth}(x)\ge d]=\frac1s\sum_{m=0}^{s-1}\Pr[\mathrm{depth}(x)\ge d\mid R=m]\le\frac1s\sum_{m=0}^{s-1}\Pr[D_m\ge d]=\Pr[D_s-1\ge d].
```
The bounds on the moments follow because $f(d)=d$ and $f(d)=a^d$ for $a\ge1$ are non-decreasing, and the proof of [Theorem 2](#thm-oec-tail) only uses upper bounds on $\mathbb{E}[a^{\mathrm{depth}(x)}]$ for $a\ge1$.

We use induction on $m$; the case $m=0$ is as before. Let $m\ge 1$, so $z$ is not a root and its parent $y=p(z)$ is a vertex of $C(z)$ higher than $z$, hence $y\in \bar\lambda$. By [Lemma 5](#lem-simulation) applied to the prefix $\lambda$, the simulated parent $a=\hat p(z)$ lies in $\bar\lambda$, it lies in $C(z)$, and $\pi(y)\ge\pi(a)$. Moreover, $a$ is determined by $\lambda$. By [Lemma 1](#lem-suffix), the number $U$ of vertices of $C(z)$ higher than $a$ is uniformly distributed on $`\{0,\dots,m-1\}`$, conditioned on $\Lambda_k=\lambda$. The number $K$ of vertices of $C(z)$ higher than $y$ satisfies $K\le U$.

Next, we refine the conditioning to the prefix that ends with $y$, as in the proof of [Lemma 4](#lem-law). Here $y$ is not determined by $\lambda$, but for every $j>k$, the prefix $\Lambda_j$ determines whether the position $k'=\pi(y)$ of $y$ equals $j$. Indeed, by the last statement of [Lemma 5](#lem-simulation), for every $i\ge k$, the parent of $z$ belongs to $\Lambda_i$ if and only if the frozen simulation for $\Lambda_i$ reports a parent in $\Lambda_i$; and $k'=j$ means that this holds for $i=j$ but not for $i=j-1$. Let $\mathcal{E}$ be the set of sequences $\lambda'$ of distinct vertices that extend $\lambda$ and for which $\Lambda_{|\lambda'|}=\lambda'$ implies $k'=|\lambda'|$, and let $A_{\lambda'}$ denote the event $\Lambda_{|\lambda'|}=\lambda'$. Conditioned on $\Lambda_k=\lambda$, exactly one of the events $A_{\lambda'}$ with $\lambda'\in\mathcal{E}$ occurs, namely the one with $\lambda'=\Lambda_{k'}$, and each of them implies $\Lambda_k=\lambda$. On $A_{\lambda'}$, the vertex $y$ is the last vertex of $\lambda'$, and $K=j(\lambda')$, where $j(\lambda')=|C(z)\setminus\lambda'|<m$. By the induction hypothesis applied to $\lambda'$, conditioned on $A_{\lambda'}$, we have $\mathrm{depth}(y)\preceq D_{j(\lambda')}$.

Let $f$ be non-decreasing, and let $g(j)=\mathbb{E}[f(1+D_j)]$, which is non-decreasing by [Lemma 6](#lem-monotone). Since $\mathrm{depth}(z)=1+\mathrm{depth}(y)$, the law of total expectation and $K\le U$ give
```math
\begin{aligned}
  \mathbb{E}\bigl[f(\mathrm{depth}(z)) \mid \Lambda_k=\lambda\bigr]
  &= \sum_{\lambda'\in\mathcal{E}}\Pr\bigl[A_{\lambda'} \mid \Lambda_k=\lambda\bigr]\cdot\mathbb{E}\bigl[f(1+\mathrm{depth}(y)) \mid A_{\lambda'}\bigr] \\
  &\le \sum_{\lambda'\in\mathcal{E}}\Pr\bigl[A_{\lambda'} \mid \Lambda_k=\lambda\bigr]\cdot g\bigl(j(\lambda')\bigr)
  = \mathbb{E}\bigl[g(K) \mid \Lambda_k=\lambda\bigr] \\
  &\le \mathbb{E}\bigl[g(U) \mid \Lambda_k=\lambda\bigr]
  = \mathbb{E}[f(D_m)],
\end{aligned}
```
where the last equality holds because $U$ is uniform on $`\{0,\dots,m-1\}`$. Hence $\mathrm{depth}(z)\preceq D_m$ conditioned on $\Lambda_k=\lambda$. $\square$

#### The weighted DSU

<a name="sec-wdsu"></a>

**Remark 1 (The weighted DSU).**

<a name="rem-wdsu"></a>

The weighted DSU attaches the endpoint with the smaller random index. We model the indices as the uniformly random order $\pi$. Both operations first test $u=v$, before either root walk. In that case an insertion ignores the self-loop and a query returns `std::nullopt`. For distinct endpoints and distinct edge weights, insertion uses the following three steps in addition to the priority rule. [not_amir (2024)](#ref-notamir2024wdsu)

- **(a)** Walk from both insertion endpoints to their roots. If the roots differ, $\theta=\mathrm{PM}_{G_{t-1}}(u,v)=+\infty$ and the stitch starts directly on the two trees. If the roots agree, the shortcuts in (c) leave strictly increasing weights along both root paths. Advance the pointer with the lighter parent edge until the pointers are adjacent. Their connecting edge is the heaviest edge on the tree path, of weight $\theta=\mathrm{PM}_{G_{t-1}}(u,v)$. If $w>\theta$, discard the new edge. If $w<\theta$, cut the heaviest edge at its child $z^*$ by setting $p(z^*)\gets\bot$ and $\mathit{wt}(z^*)\gets+\infty$, then stitch the two resulting trees. In both connected cases the discarded edge is heaviest on the cycle created by the insertion.
- **(b)** Attaching a root completes the stitch. If the attached vertex $x$ had a parent $x'$, the implementation continues with pointers at $x'$ and $x$. The next climb from $x$ reaches $y$, since $\mathit{wt}(x)=w<w'$ for the new pending weight $w'$. The other pointer climbs in a different tree, so the two descriptions make the same sequence of attachments.
- **(c)** When a root walk or climb leaves $z$, repeatedly set $p(z)\gets p(p(z))$ while $\mathit{wt}(p(z))\le\mathit{wt}(z)$. Each shortcut slides the edge from $z$ along a parent edge of at most its weight and replaces its parent by a higher vertex in the same tree. The climb with pending weight $w$ stopped at $y$ before the attachment, so afterwards $\mathit{wt}(y)>w=\mathit{wt}(x)$. The extra step from $x$ in (b) therefore reaches $y$ directly. [TW transformations](#tw-transformations)

A query with distinct endpoints performs the two root walks with the shortcuts in (c). It returns `std::nullopt` if the roots differ and otherwise returns the weight found by the path-max scan in (a). [Facts 1](#fact-correct) and [2](#fact-heap) hold after every operation. The root walks, the path-max scan, and the stitch each advance along the original ancestor paths. Each shortcut skips a vertex and each attachment incurs at most one extra step as in (b). Thus an operation on $u$ and $v$ takes $O(1+\mathrm{depth}(u)+\mathrm{depth}(v))$ time.

For the weighted DSU, $\sigma$ denotes the whole sequence of operations, including queries. A *time* is a moment between two operations, while the index $t$ in $e_t$ and $G_t$ counts insertions. The frozen simulation keeps the fields $\hat p(v)$ and $\widehat{\mathit{wt}}(v)$ only for $v\in\lambda$. It skips every operation with $u=v$ before its first walk. For distinct endpoints it applies the following rules.

- *Walks and climbs.* A pointer located at a vertex of $\bar\lambda$ never moves, and a walk to the root ends when its pointer reaches a root or a vertex of $\bar\lambda$. Whenever a walk to the root or a climb leaves a vertex $z\in\lambda$, the simulation first sets $\hat p(z)\gets\hat p(\hat p(z))$ repeatedly as long as $\hat p(z)\in\lambda$ and $\widehat{\mathit{wt}}(\hat p(z))\le\widehat{\mathit{wt}}(z)$, and then moves the pointer to $\hat p(z)$. The walk to the lowest common ancestor in (a) changes no field and is not simulated.

- *Cut (insertions only).* The value $`\theta=\mathrm{PM}_{G_{t-1}}(u,v)`$ depends only on $\sigma$. If $w<\theta<+\infty$ and some $z\in\lambda$ has $\widehat{\mathit{wt}}(z)=\theta$, the simulation sets $\hat p(z)\gets\bot$ and $\widehat{\mathit{wt}}(z)\gets+\infty$.

- *Stitch (insertions only).* If $w<\theta$, the simulation proceeds as the frozen simulation of the priority-based stitch, with the climbs as above. This includes distinct components, for which $\theta=+\infty$.

As before, the lower of the two pointers of a stitch is determined by $\lambda$ whenever at least one of them lies in $\lambda$, so this simulation is a function of $\sigma$ and $\lambda$ only.

**Lemma 7 (Simulation for the weighted DSU).**

<a name="lem-wdsu"></a>

For the weighted DSU and the frozen simulation above, (i)–(iii) of [Lemma 5](#lem-simulation) hold at every time and for every $v\in\lambda$. Consequently, [Theorem 3](#thm-two-sided) holds for the weighted DSU as well.

**Proof.**

We compare the two executions on the current fields, including during an operation, because the walks in (a) can change fields that the stitch visits again. During insertion $e_t$ the graph in (iii) is $G_t$. For $u=v$ both executions return before reading a parent pointer. A self-loop preserves the graph components, so all three invariants hold through a same-endpoint operation.

For distinct endpoints we use two properties of the true execution. First, the forest is heap-ordered, and each parent lies in the same graph component as its child. Second, each finite weight occurs at most once in the forest. The inserted weights are distinct, a shortcut preserves the weights, a cut removes one, and an attachment exchanges $\mathit{wt}(x)$ with the pending weight.

*Walks and climbs.* At a revealed vertex $z$ with a revealed parent $q$, both executions see the same two weights, so both make the same shortcut decision. After a shortcut, the new parents of $z$ are the respective parents of $q$, so they satisfy (ii) or (iii) by the invariant for $q$. In particular, a new simulated parent stays in the component of $q$ and $z$. Once the simulated parent is unrevealed, only the true execution can apply another shortcut at $z$, which can only raise the true parent's priority. After advancing both pointers, they agree inside $\lambda$ or are both unrevealed with the true pointer at least as high. A true pointer outside $\lambda$ can only visit other unrevealed vertices.

*Cut.* When $w<\theta<+\infty$, the true cut identifies $z^*$ by its unique weight $\theta$. If $z^*$ is revealed, (i) gives the same weight in the simulation, and both runs make $z^*$ a root. If it is unrevealed, no simulated weight equals $\theta$, so the revealed fields stay intact.

*Stitch.* Every stitch starts in distinct trees, either because the roots differed in (a) or because of the cut. When $x$ is attached below $y$, its old parent, if present, stays in the other tree. Thus each pending edge again joins distinct trees. By (b) and (c), the extra step from $x$ to $y$ preserves the fields and the order of attachments. The attachment argument of [Lemma 5](#lem-simulation) applies to the current fields. By (i), both runs make the same climb decision at a revealed vertex. After advancing, the pointers coincide inside $\lambda$ or are both outside it with the true pointer at least as high. Hence every attachment preserves (i)–(iii).

This proves (i)–(iii) at every time, and the last statement of [Lemma 5](#lem-simulation) follows from them as before. With these simulation properties, [Facts 1](#fact-correct) and [2](#fact-heap), the uniform suffix in [Lemma 1](#lem-suffix), and the monotonicity in [Lemma 6](#lem-monotone), the proof of [Theorem 3](#thm-two-sided) gives the same depth bound for the weighted DSU. Between insertions, queries can only shorten root paths by shortcutting to ancestors. Hence the union bound over the $M$ insertion times in [Corollary 1](#cor-ops) covers the intervening queries as well. $\square$

### The oblivious adversary is necessary

<a name="sec-adaptive"></a>

All our bounds assume that the insertion sequence does not depend on the priorities. This assumption cannot be dropped.

**Proposition 1.**

<a name="prop-adaptive"></a>

There exists an insertion sequence of length $n-1$, depending on $\pi$, after which the OEC-forest (and also the priority-based stitch) is a path of length $n-1$.

**Proof.**

Let $\ell_1,\dots,\ell_n$ be the vertices in increasing priority order. Insert $(\ell_i,\ell_{i+1})$ for $i=1,\dots,n-1$ with distinct weights. At step $i$, $\ell_i$ is the root of the component $\{\ell_1,\dots,\ell_i\}$ and $\ell_{i+1}$ is an isolated root. Both algorithms attach $\ell_i$ directly below $\ell_{i+1}$. No other parent changes. In the resulting path, a query or further insertion between $\ell_1$ and $\ell_n$ takes $\Theta(n)$ time. $\square$

Here the adversary uses $\pi$ to choose the edges. In contrast, if the $n-1$ edges of a path are chosen independently of $\pi$ and inserted in any fixed order, then by [Theorem 1](#thm-oec-exact) every vertex has expected depth $H_n-1$ after the last insertion.

## The Rival-Tree

<a name="sec-hb"></a>

The strict AM-tree of Section 3 of the paper bounds the subtree size of every child relative to its parent. The Rival-tree balances the heights of the children, as in an AVL tree. A tall child needs a rival of almost the same height. It guarantees height at most $\log_\varphi(n+1)<1.441\log_2(n+1)$, where $\varphi=(1+\sqrt5)/2$. With either link, a path-max query takes $O(\log n)$ worst-case time. With array-indexed child lists, an insertion takes $O(\log n)$ amortized and $O(\log^2 n)$ worst-case time. The stitch analysis assumes distinct weights.

### Algorithms

<a name="sec-rival-algorithms"></a>

**Promotions.**

A *promotion* moves $x$ one level up past its parent $y$ using a TW transformation. [Ding et al. (2025)](#ref-ding2025amtree)

- *Shortcut*. If $y$ has a parent and $\mathit{wt}(x)\ge\mathit{wt}(y)$, make $x$ a child of $p(y)$, keeping $\mathit{wt}(x)$.
- *Rotate*. Otherwise $x$ takes the place of $y$ with weight $\mathit{wt}(y)$, and $y$ becomes a child of $x$ with weight $\mathit{wt}(x)$.

Only $x$, $y$, and the old parent of $y$ can change their children. $\mathrm{Perch}(x)$ promotes $x$ until it becomes a root.

<a name="alg-perch"></a>

**Algorithm 1. The perch-based link.**

```
Perch(x):
    while p(x) != null: Promote(x)

PerchLink(u, v, w):
    if u == v: return
    Perch(u)
    Perch(v)
    if p(u) == v: wt(u) = min(wt(u), w)
    else: p(u) = v; wt(u) = w
```

<a name="alg-stitch"></a>

**Algorithm 2. The stitch-based link (distinct weights).**

```
Stitch(u, v, w):
    if u == v or u == null or v == null: return
    if w > wt(u): Stitch(p(u), v, w)
    else if w > wt(v): Stitch(u, p(v), w)
    else:
        x = the endpoint of smaller stored height (ties arbitrary)
        y = the other endpoint
        old_parent = p(x); old_weight = wt(x)
        p(x) = y; wt(x) = w
        Stitch(old_parent, y, old_weight)
```

<a name="alg-hb"></a>

**Algorithm 3. The Rival-tree.**

```
Insert(u, v, w):
    if u == v: return
    P = {u, v} together with their ancestors before the link
    PerchLink(u, v, w) or Stitch(u, v, w) with the height rule
    Rebalance(P)

Rebalance(P):
    for y in P, in non-increasing depth immediately after the link:
        UpdateHeight(y)
        Fix(y)

Fix(y):
    while y is unbalanced:
        promote the heavy child of y
        UpdateHeight(y)

UpdateHeight(y):
    k = 0 if y has no children, else 1 + the largest stored child height
    if k == h(y): return
    if p(y) != null: remove y from the list for (p(y), h(y))
    h(y) = k
    if p(y) != null: add y to the list for (p(y), k)
```

A path-max query on $u$ and $v$ walks from $u$ and $v$ to their lowest common ancestor and takes $O(1+\mathrm{depth}(u)+\mathrm{depth}(v))$ time.

### Balance and height

For a vertex $y$, let $c_1(y)\ge c_2(y)$ be the largest and the second largest height of a child of $y$, where a missing value is $-1$; thus $h(y)=c_1(y)+1$.

**Definition 1.**

<a name="def-hb"></a>

A vertex $y$ is *balanced* if $c_1(y)\le 1$ or $c_2(y)\ge c_1(y)-1$. Otherwise $y$ is *unbalanced*; then $y$ has a unique tallest child, called its *heavy child*, whose height is at least $2$, and every other child of $y$ has height at most $h(y)-3$. A forest is *balanced* if all its vertices are balanced. The *Rival-tree* is a T-MST that is kept balanced.

In other words, a vertex of height $k\ge3$ is balanced if it has a child of height $k-1$ and another child, a *rival*, of height at least $k-2$; we call this condition the *rival rule*. Whether $y$ is balanced, and its height, depend only on the multiset of the heights of its children. The following observation is immediate.

**Fact 1.**

<a name="fact-hb-extra"></a>

Let $O$ be a nonempty set of children of $y$ whose heights, as a multiset, satisfy the condition of [Definition 1](#def-hb), and let $k$ be the largest height in $O$. If every other child of $y$ has height at most $k$, then $y$ is balanced and $h(y)=k+1$.

**Lemma 1.**

<a name="lem-hb-height"></a>

Every tree of a balanced forest with $s$ vertices has height at most $\log_\varphi(s+1)$.

**Proof.**

Let $S(k)$ be the minimum number of vertices of a balanced tree of height $k$. Clearly $S(0)=1$, $S(1)=2$, $S(2)=3$, and $S$ is increasing, since the subtree of a tallest child of the root is a balanced tree of height $k-1$. For $k\ge3$, the root has a child of height $k-1$ and another child of height at least $k-2$, so $S(k)\ge S(k-1)+S(k-2)+1$. Let $F_i$ denote the Fibonacci numbers, with $F_1=F_2=1$. Then $S(k)+1\ge F_{k+2}$ for every $k\ge1$: this holds for $k=1,2$, since $3\ge F_3$ and $4\ge F_4$, and for $k\ge3$ we have $S(k)+1\ge (S(k-1)+1)+(S(k-2)+1)\ge F_{k+1}+F_k=F_{k+2}$. Since $F_{k+2}\ge\varphi^k$, a balanced tree of height $k\ge1$ with $s$ vertices satisfies $s+1\ge\varphi^k$. $\square$

### Rebalancing

**Data layout.**

<a name="rival-layout"></a>

Every vertex $v$ stores $p(v)$, $\mathit{wt}(v)$, $h(v)$, and a doubly linked list for each stored child height. A bit mask identifies the nonempty lists. When a child changes parents, move its list entry from the old parent to the new one, using its current stored height. Roots have no such entry. When only a child's stored height changes, `UpdateHeight` moves its entry between lists of the same parent. Each list operation updates the corresponding bit when the list gains its first child or loses its last. [Algorithm 3](#alg-hb)

All heights are $O(\log n)$, so each mask fits in $O(1)$ machine words. With an array of list heads per vertex, finding a tallest child, computing the height, testing balance, and adding or removing a child each take $O(1)$ worst-case time. This layout uses $O(n\log n)$ space and gives the insertion time bounds in the two theorems below.

There are at most $n-1$ nonempty (parent, height) lists, so storing their heads in a hash table takes $O(n)$ space. If table operations take expected $O(1)$ time for inputs independent of the random hash function, the expected time is $O(\log n)$ amortized per insertion and $O(\log^2 n)$ for any single insertion. The C++ option `kHash=true` uses linear probing with a fixed multiplicative hash. A time guarantee for that implementation requires a bound on its probe lengths for the keys encountered. Path-max queries use only parent pointers and edge weights, so their worst-case bound applies to both layouts. [RivalTree implementation](include/tmst/rival_tree.h)

**The rebalancing procedure.**

The procedure $\mathrm{Fix}(y)$ in [Algorithm 3](#alg-hb) promotes the heavy child of $y$ as long as $y$ is unbalanced, and $\mathrm{Rebalance}(P)$ calls $\mathrm{Fix}$ on every vertex of a set $P$, from the bottom up. In the analysis we use the potential
```math
\Psi=\sum_{v\in V}h(v).
```

**Lemma 2.**

<a name="lem-hb-fix"></a>

Let $y$ be a vertex whose proper descendants are all balanced, and let $z=p(y)$ before the call $\mathrm{Fix}(y)$.

- **(i)** Every promotion in [the promotion step of Fix](#alg-hb) decreases $h(y)$ by at least $2$ and does not increase the height of any vertex; hence it decreases $\Psi$ by at least $2$, and $\mathrm{Fix}(y)$ performs at most $h(y)/2$ promotions.

- **(ii)** After the call, $y$ and all vertices that were descendants of $y$ before the call are balanced.

- **(iii)** The call does not change the parent of any vertex that was not a descendant of $y$, and among these vertices, only $z$ may change its set of children.

**Proof.**

Let $R$ be the children of $y$ promoted by rotation so far. The current parent $q$ of $y$ is either its original parent $z$ or a vertex in $R$. Maintain that every proper descendant of $y$ is balanced and each $r\in R$ keeps its original children and height. Its additional children have height below $h(r)$.

In one iteration the heavy child $c$ has height $k\ge2$, so $h(y)=k+1$ and all other children have height at most $k-2$. Removing $c$ drops $h(y)$ to at most $k-1$. A shortcut leaves the subtree of $c$ intact. A rotation also gives $c$ the child $y$ of height at most $k-1$, so its original children still witness its height and balance by [Fact 1](#fact-hb-extra).

The parent $q$, if present, replaces a child of height $k+1$ by children of height at most $k$. When $q\in R$, those children are all shorter than $q$ and its original children still witness balance. When $q=z$, the maximum child height can only decrease, so the same is true of $h(z)$ and the heights of its ancestors. This establishes (i) and preserves the invariant until $y$ is balanced, giving (ii). Only $y$ and the promoted children change parents. The changed child sets belong to those vertices and $z$, giving (iii). $\square$

**Lemma 3.**

<a name="lem-hb-rebalance"></a>

Suppose that $P$ is closed under taking parents (the parent of every vertex of $P$ is in $P$ or is $\bot$), and that every unbalanced vertex and every vertex whose stored height is wrong lies in $P$. Then after $\mathrm{Rebalance}(P)$, the forest is balanced and all stored heights are correct.

**Proof.**

The root paths of unprocessed vertices in $P$ stay as they were just after the link. Indeed, a fix at $y$ changes parents only in its subtree. If a remaining vertex's root path entered that subtree, it would have been deeper than $y$ after the link and processed earlier.

When $y$ is reached, its proper descendants are therefore outside $P$ or already processed. Inductively they are balanced and have correct stored heights. Recomputing $h(y)$ and applying [Lemma 2](#lem-hb-fix) fixes $y$ and preserves those descendants' heights and balance. The only vertices above the original subtree whose heights or balance can change lie on the original ancestor path of $y$. They belong to $P$ and are still unprocessed. This proves the invariant and the result. $\square$

### Insertion with the perch-based link

The next lemma bounds the effect of a perch on the heights.

**Lemma 4.**

<a name="lem-hb-perch"></a>

Suppose that $\mathrm{Perch}(x)$ promotes $x$ itself until it becomes a root. Let $a_1,\dots,a_k$ be the proper ancestors of $x$ before $\mathrm{Perch}(x)$, from the bottom up, and let $`h^*`$ be the height of the tree of $x$. After $\mathrm{Perch}(x)$, the height of every $a_i$ is at most its height before the call, every other vertex except $x$ has the same subtree as before, and $`h(x)\le h^*+1`$.

**Proof.**

A promotion of $x$ with parent $y$ changes only the parents of $x$ and, in the case of a rotation, of $y$. Hence $p(a_{i})=a_{i+1}$ until the $i$-th promotion, which therefore promotes $x$ over $a_i$, and after it $a_i$ never changes again. Let $a_0=x$. Right before the $i$-th promotion, the children of $a_i$ are its original children, except that $a_{i-1}$ is replaced by $x$ if the $(i-1)$-th promotion was a rotation, and $x$ is added if it was a shortcut (for $i=1$, the child $a_0=x$ is simply present). The $i$-th promotion removes $x$ from the children of $a_i$. So the final children of $a_i$ are its original children other than $a_{i-1}$, together with $a_{i-1}$ if $i\ge2$ and the $(i-1)$-th promotion was a shortcut. By induction on $i$, the final height of $a_{i-1}$ is at most its original height, and hence so is the final height of $a_i$. The final children of $x$ are its original children and the vertices $a_i$ over which $x$ was rotated, all of which have height at most $`h^*`$; so $`h(x)\le h^*+1`$. Every other vertex keeps its children, and its subtree contains neither $x$ nor any $a_i$, so its subtree does not change. $\square$

**Theorem 1.**

<a name="thm-hb-perch"></a>

Consider the Rival-tree with the perch-based link, starting from $n$ isolated vertices. After every insertion, the forest is a balanced T-MST of the current graph, so its height is at most $\log_\varphi(n+1)$, and every path-max query takes $O(\log n)$ worst-case time. With the array of list heads, every insertion takes $O(\log n)$ amortized time and $O(\log^2 n)$ worst-case time. [Data layout](#rival-layout)

**Proof.**

All promotions preserve path-max values. If $u$ and $v$ were connected, the last promotion of $v$ rotates it over the root $u$, so $u$ is now its child. The link keeps the lighter of the two parallel edges between them. If they were disconnected, it joins the two perched roots by the new edge. [TW transformations](#tw-transformations)

Initially $P$ is closed under taking parents. Promoting $x\in P$ changes only children sets and parent pointers involving $x$, its parent, and its grandparent, which are in $P$ or null. The final attachment also stays within $P$. Thus subtrees outside $P$ retain their original balance and heights until the fix pass, and [Lemma 3](#lem-hb-rebalance) restores balance everywhere.

Let $L=\log_\varphi(n+1)$. Both old root paths have length at most $L$ by [Lemma 1](#lem-hb-height), so $|P|\le2L+2$. $\mathrm{Perch}(u)$ takes at most $L$ promotions and increases the tree height to at most $L+1$ by [Lemma 4](#lem-hb-perch). The second perch takes at most $L+1$ promotions. After both perches and the final attachment, all heights are $O(L)$. Compute post-link depths by memoization and bucket-sort $P$ in $O(|P|)$ time. Each fix takes $O(L)$ promotions by [Lemma 2](#lem-hb-fix), giving the $O(L^2)$ worst-case insertion bound.

By [Lemma 4](#lem-hb-perch), each perch contributes $O(L)$ to $\Psi$ through its endpoint, and the final attachment can increase only the new parent's height. Thus the link takes $O(L)$ time and adds $O(L)$ potential. Each promotion in a fix decreases $\Psi$ by at least two, and the other rebalancing work takes $O(|P|)$. Scaling the nonnegative potential, which starts at zero, gives $O(L)$ amortized insertion time. $\square$

### Insertion with the stitch-based link

Use the smaller pre-link height at the selection step of [Algorithm 2](#alg-stitch), breaking ties arbitrarily. At the start of a call, suppose both pointers still have their original root paths. If one is an ancestor of the other, it was also an ancestor before the link and has greater stored height. The rule therefore chooses the descendant, which establishes the induction step in the shared proof. [Stitch correctness](#stitch-correctness)

A stitch can increase many heights, so we analyze the link and fixes together. Write $r(v)$ for the height before insertion, $T'$ for the forest right after the link, and $T'[y]$ for the subtree of $y$ in $T'$. A vertex $y\in P$ is *processed* when `Rebalance` recomputes its height just before $\mathrm{Fix}(y)$.

**Lemma 5.**

<a name="lem-hb-stitch-height"></a>

In an insertion with the stitch-based link, every vertex $y\in P$ satisfies $h(y)\le r(y)+7$ when it is processed.

**Proof.**

We first collect four properties of the insertion.

- **(a)** Every edge of $T'$ satisfies $r(\text{child})\le r(\text{parent})$: an edge of the forest before the insertion joins a vertex to a taller parent, and an attached vertex is the endpoint of smaller stored height, which is its height before the insertion. Since heights strictly increase along every root path before the insertion, every value of $r$ occurs at most twice in $P$.

- **(b)** A vertex outside $P$ keeps its height and changes its parent only when promoted by $\mathrm{Fix}$ at that parent. The link leaves its subtree intact. In a fix it can only be a proper descendant of the fixed vertex, in which case its height stays the same, or be outside the fix altogether. The original ancestors of every fixed vertex lie in $P$ by the root-path invariant of [Lemma 3](#lem-hb-rebalance).

- **(c)** When $y$ is processed, its subtree has the same vertex set as in $T'$. A fix can change subtree membership only at the vertex it fixes and that vertex's proper descendants. By the bottom-up invariant, no earlier fix was at an ancestor of $y$.

- **(d)** Let $a\in P$, and let $C_a$ be the set of children of $a$ before the insertion that lie outside $P$ and are promoted in the call $\mathrm{Fix}(a)$. Then $\sum_{c\in C_a}\varphi^{h(c)}\le\varphi^{r(a)}$. Indeed, every promotion in $\mathrm{Fix}(a)$ promotes a child of height $h(a)-1$ and then decreases $h(a)$ by at least $2$ ([Lemma 2](#lem-hb-fix)(i)), so the heights of the promoted children differ pairwise by at least $2$. By (b), every $c\in C_a$ has height $r(c)\le r(a)-1$, so the sum is at most $\varphi^{r(a)-1}\sum_{j\ge0}\varphi^{-2j}=\varphi^{r(a)-1}/(1-\varphi^{-2})=\varphi^{r(a)}$.

Let $A=P\cap V(T'[y])$ and $W=\sum_{a\in A}\varphi^{r(a)}$. By (a), every $a\in A$ satisfies $r(a)\le r(y)$, and every value occurs at most twice, so
```math
W\le 2\sum_{j=0}^{r(y)}\varphi^j<2\varphi^{r(y)+2},
```
using $\sum_{j=0}^{k}\varphi^j=(\varphi^{k+1}-1)/(\varphi-1)<\varphi^{k+2}$.

If $y$ has no children, then $h(y)=0$. Otherwise, let $c_0$ be a tallest child of $y$ and $H=h(c_0)$, so $h(y)=H+1$. When $y$ is processed, all its proper descendants are balanced (proof of [Lemma 3](#lem-hb-rebalance)). We build a *witness* tree $Q$ rooted at $c_0$: a vertex $x$ of $Q$ is *expanded* if $x\in P$ and $h(x)\ge3$, and then its two tallest children are its children in $Q$; otherwise $x$ is a leaf of $Q$. An expanded vertex $x$ of height $k\ge3$ is balanced, so its two tallest children have heights $k-1$ and at least $k-2$, and $\varphi^{k}=\varphi^{k-1}+\varphi^{k-2}$ is at most the sum of $\varphi^{h}$ over these two children. By induction over $Q$,
```math
\varphi^H\le\sum_{\text{leaves }x\text{ of }Q}\varphi^{h(x)}.
```
By (c), all vertices of $Q$ lie in $V(T'[y])$, so the expanded vertices lie in $A$, and the current parent of every leaf is $y$ or an expanded vertex, hence lies in $A$. We split the leaves $x$ into three classes.

- If $x\in P$, then $h(x)\le2$ and $x\in A$, so these leaves contribute at most $\varphi^2W$.

- If $x\notin P$ and its current parent $a$ is its parent before the insertion, then $h(x)=r(x)\le r(a)-1$ by (b), and every $a\in A$ is the parent of at most two leaves of $Q$; these leaves contribute at most $2\varphi^{-1}W$.

- If $x\notin P$ and its current parent differs from its parent $a$ before the insertion, then $a$ is also the parent of $x$ in $T'$, so $a\in V(T'[y])$, since $x\ne y$ lies in $V(T'[y])$. By (b), $x$ first left $a$ when it was promoted in $\mathrm{Fix}(a)$, so $a\in P$, hence $a\in A$, and $x\in C_a$. By (d), these leaves contribute at most $\sum_{a\in A}\varphi^{r(a)}=W$.

Since $\varphi^2+2\varphi^{-1}+1=(\varphi+1)+2(\varphi-1)+1=3\varphi$, we get
```math
\varphi^H\le3\varphi W<6\varphi^{r(y)+3}<\varphi^{r(y)+7},
```
because $\varphi^4>6$. Hence $H\le r(y)+6$ and $h(y)\le r(y)+7$. $\square$

The bound holds only when $y$ is processed: right after the link, before the descendants of $y$ are fixed, the height of $y$ can be much larger.

**Theorem 2.**

<a name="thm-hb-stitch"></a>

With the stitch-based link and distinct weights, the forest is a balanced T-MST after every insertion, so its height is at most $\log_\varphi(n+1)$, and every path-max query takes $O(\log n)$ worst-case time. With the array of list heads, every insertion takes $O(\log n)$ amortized time and $O(\log^2 n)$ worst-case time. [Data layout](#rival-layout)

**Proof.**

The stitch advances along the original paths in $O(|P|)$ time. It changes only parent and children sets in $P$, keeps $P$ closed under taking parents, and preserves path-max values. The post-link paths have at most $|P|$ vertices in $P$ and $O(L)$ vertices outside it. Thus the balance and worst-case time arguments of [Theorem 1](#thm-hb-perch) apply.

For the amortized time, let $k_y$ be the number of promotions in $\mathrm{Fix}(y)$ for $y\in P$, and $K=\sum_{y\in P}k_y$. By [Lemma 5](#lem-hb-stitch-height) and [Lemma 2](#lem-hb-fix)(i), $h(y)\le r(y)+7-2k_y$ after the call $\mathrm{Fix}(y)$. The height of $y$ does not change afterwards: by the proof of [Lemma 3](#lem-hb-rebalance), a later call $\mathrm{Fix}(a)$ changes, among the descendants of $a$, only the height of $a$, and otherwise only the heights of ancestors of $a$, which are unprocessed. Vertices outside $P$ keep their heights ([Lemma 5](#lem-hb-stitch-height), property (b) in its proof). Hence the insertion changes $\Psi$ by
```math
\Delta\Psi=\sum_{y\in P}\bigl(h(y)-r(y)\bigr)\le 7|P|-2K,
```
where $h(y)$ is the height after the insertion. The insertion takes $O(|P|+K)$ time: the link takes $O(L)$ time, the depths and the order of `Rebalance` take $O(|P|)$ time, and $\mathrm{Fix}(y)$ takes $O(1+k_y)$ time. With the potential $c\Psi$ for a suitable constant $c$, the amortized time is $O(|P|)=O(L)$, since $|P|\le2L+2$. Since $\Psi\ge0$ and $\Psi=0$ initially, every insertion takes $O(\log n)$ amortized time. $\square$

The $O(\log^2 n)$ worst-case bound of [Theorem 2](#thm-hb-stitch) is tight for a single insertion.

**Proposition 1.**

<a name="prop-hb-staircase"></a>

For infinitely many $n$, there is a balanced forest with $n$ vertices and an edge whose stitch-based insertion performs $\Omega(\log^2 n)$ promotions.

**Proof.**

For $k\ge0$, let $B_k$ be the balanced tree of height $k$ defined by $B_0$ being a single vertex, $B_1$ a path with two vertices, $B_2$ a path with three vertices, and $B_k$ for $k\ge3$ a root with two children whose subtrees are $B_{k-1}$ and $B_{k-2}$. Its size is $|B_k|=O(\varphi^k)$.

Fix $\ell\ge2$ and let $h_i=\ell+i$ for $i=0,\dots,\ell$. The first tree has a path $x_0,x_1,\dots,x_\ell$, where $x_\ell$ is the root and $p(x_{i-1})=x_i$. The vertex $x_0$ has two further subtrees $B_{\ell-1}$ and $B_{\ell-2}$, so $h(x_0)=h_0$. For $i\ge1$, the vertex $x_i$ has, besides $x_{i-1}$, children that are roots of copies of $B_{h_i-2},B_{h_i-4},\dots$, down to height $0$ or $1$; we call them the *stairs* of $x_i$. By induction, $h(x_i)=h_i$, and $x_i$ is balanced, since it has children of heights $h_i-1$ and $h_i-2$. The weights increase along the path, $\mathit{wt}(x_0)<\mathit{wt}(x_1)<\dots<\mathit{wt}(x_{\ell-1})$, and all other weights are arbitrary and distinct. The second tree is a copy of $B_{2\ell+1}$ with root $b$. The forest has $n=O(\varphi^{2\ell})$ vertices, so $\ell=\Omega(\log n)$.

Insert an edge $(x_0,b,w)$ with $w<\mathit{wt}(x_0)$. In every call of [Algorithm 2](#alg-stitch), one endpoint is $b$, which is a root, and the other is some $x_i$, with pending weight $w$ or $\mathit{wt}(x_{i-1})$, which is smaller than $\mathit{wt}(x_i)$. So no endpoint climbs, and since $h(x_i)\le 2\ell<h(b)$, the vertex $x_i$ is attached below $b$, and the next call has the endpoints $x_{i+1}$ and $b$. This repeats until the root $x_\ell$ is attached below $b$. Now every $x_i$ with $i\ge1$ has lost its child $x_{i-1}$ and has only its stairs as children, whose heights differ by $2$. As long as its tallest remaining stair has height at least $2$, the vertex $x_i$ is unbalanced and $\mathrm{Fix}(x_i)$ promotes that stair, which removes it from the children of $x_i$ by either kind of promotion. Since the stairs are not affected by the processing of other vertices, $\mathrm{Fix}(x_i)$ performs one promotion for each stair of height at least $2$, that is, $\lfloor h_i/2\rfloor-1$ promotions. In total, the insertion performs $\sum_{i=1}^{\ell}(\lfloor h_i/2\rfloor-1)=\Omega(\ell^2)=\Omega(\log^2 n)$ promotions. $\square$

The two starting trees are reachable from singletons under the implementation's tie rule. Insert their edges by increasing depth of the new child, using strictly decreasing weights. Each partial tree contains every level above the deepest one entered, so the two final tallest children retain heights within one of the current maximum wherever the rival rule applies. The new vertex is a leaf, and with the new vertex first the stitch chooses it even on a height tie. The new weight is lighter than the old parent edges, so the new leaf attaches to its intended parent.

## References

<a name="ref-ding2025amtree"></a>

Xiangyun Ding, Yan Gu, and Yihan Sun (2025). *New Algorithms for Incremental Minimum Spanning Trees and Temporal Graph Applications*. arXiv `2504.04619`.

<a name="ref-knuth1997taocp1"></a>

Donald E. Knuth (1997). *The Art of Computer Programming, Volume 1: Fundamental Algorithms*. 3rd edition, Addison-Wesley, Section 1.2.10.

<a name="ref-notamir2024wdsu"></a>

not_amir (2024). *Maintaining MST with Online Edge Insertions — No LCT Needed*. Codeforces blog entry `130107`.

<a name="ref-song2024oec"></a>

Jingyi Song, Dong Wen, Lantian Xu, Lu Qin, Wenjie Zhang, and Xuemin Lin (2024). “On Querying Historical Connectivity in Temporal Graphs.” *Proceedings of the ACM on Management of Data*, 2(3). DOI `10.1145/3654960`.

<a name="ref-xu2026oecjournal"></a>

Lantian Xu, Dong Wen, Lu Qin, Wenjie Zhang, Xubo Wang, and Xuemin Lin (2026). “On Querying Historical Connectivity in Large-scale Temporal Graphs.” *The VLDB Journal*, 35(1), article 2. DOI `10.1007/s00778-025-00951-7`.

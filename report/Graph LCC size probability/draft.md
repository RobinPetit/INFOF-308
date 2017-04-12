# Conjecture

The conjecture states in a general sense that there exists a recurrence relation in order to determine the cardinality of $\Lambda_k^m(V, \cdot)$.
Furthermore, this recurrence relation looks something like:

$$\abs {\Lambda_k^m(V, \cdot)} = \binom {\abs V}k\sum_{\ell=1}^{\min(m, X(k))}\sum_{p=1}^k\abs {\chi_\ell(W)}\abs {\Lambda_p^{m-\ell}(V \setminus W, \cdot)}.$$

The idea is that in order to have a graph with largest connected component of given size $k \in \mathbb N$, one requires to be able to separate (unambiguously)
the graph's LCC and its complement.

The conjecture is *probably* (yeah...) true for $p \lneqq k$ because LCC is clearly uniquely defined (only one connected component has maximum size),
but this does not fit graphs having several connected components with same size (i.e. $k$) because several graphs are counted too many times. For instance, graph

> 1---2
>
> 4---3

is counted twice: once for $W = \{1, 2\} \subset V = \{1, 2, 3, 4\}$, and once for $W = \{3, 4\} \subset V$.

To remove this redundance, two options are possible:

+ find if a given proportion is redundant, thus divide the cardinality of $\chi_\ell(W, \cdot) \times \Lambda_k^{m-\ell}(V \setminus W, \cdot)$,
+ or change the expression in order to isolate the case where $p=k$, and find the right expression (would something like $\chi_\ell(W) \times \mathcal L_k^{m-\ell}(V, W)$, for:
$$\mathcal L_k^{m-\ell}(V, W) \coloneqq \Lambda_k^{m-\ell}(V \setminus W, \cdot) \setminus \mathfrak L_k^{m-\ell}(V, W),$$
for:
$$\mathfrak L_k^{m-\ell}(V, W) \coloneqq \left\{\Gamma(V, E) \in \Lambda_k^{m-\ell}(V \setminus W, \cdot) \text{ s.t. } \text{LCC}(\Gamma(V, E)) \subset \left\{v_1, \ldots, v_{\mu(W)}\right\}\right\}$$
work knowing that:
$$\mu(W) \coloneqq \max_{i = 1, \ldots, \abs V}i\mathbb I_{\left[v_i \in W\right]}$$
?)

## Function to prove cardinality equality

To prove that two sets have equal cardinality, a bijective function must be determined between these two. If $\mathfrak Q_k^m(V)$ is the set having the right cardinality,
the function will be:

$$\Omega : \Lambda_k^m(V) \to \mathfrak Q_k^m(V) : \Gamma(V, E) \mapsto \left(\Delta_{\text{LCC}(\Gamma(V, E))}(\Gamma(V, E)), \Delta_{V \setminus \text{LCC}(\Gamma(V, E))}(\Gamma(V, E))\right).$$

This $\Omega$ function is obviously injective. Now, the right set $\mathfrak Q_k^m(V)$ needs to be found in order to be surjective (the hard point is on graphs having more than one connected component of maximum size.)

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

# 04/13

Since the set $\Lambda_l^m(V, \cdot)$ has been split again into a disjoint union of $\Lambda_{k,\alpha}^m(V, \cdot)$, it has been proven that the conjecture stands for $\alpha=1$ and the $\beta$ coefficients equal to 1.

Yet, the formula has to be proven and arranged for $\alpha > 1$. Something like the follwing could work:
$$\abs {\Lambda_{k,\alpha}^m(V, \cdot)} = \abs {\mathfrak Q_{k,\alpha}^m(V)},$$
for:
$$\mathfrak Q_{k,\alpha}^m(V) \coloneqq \bigsqcup_{(W_1, \ldots, W_\alpha) \in \mathcal P_{k,\alpha}(V)}\bigsqcup_{\stackrel {i_1, \ldots, i_\alpha}{\sum_{j=1}^\alpha i_j} \leq m}
	\left[
		\left(\prod_{j=1}^\alpha \chi_{i_j}(W_j)\right)
			\times
		\left(\bigsqcup_{p=1}^{k-1}\Lambda_p^{m-\sum_{j=1}^\alpha i_j}(V, \cdot)\right)
	\right].$$
with:
$$\mathcal P_{k,\alpha}(V) \coloneqq \left\{(W_1, \ldots, W_\alpha) \in \mathcal P(V) \text{ s.t. }
	\begin{cases}
		&\forall i \in \{1, \ldots, \alpha\} : \abs {W_i} = k \\
		&\forall (i, j) \in \{1, \ldots, \alpha\}^2 : i \neq j \Leftrightarrow W_i \cap W_j = \emptyset
	\end{cases}
\right\}$$

**Remark:** Since tuples in $\mathcal P_{k,\alpha}(V)$ are sensitive to order (i.e. $(W_1, W_2) \neq (W_2, W_1)$), we have that:
$$\abs {\mathcal P_{k,\alpha}(V)} = \alpha!\prod_{i=1}^\alpha\binom {\abs V-i(k-1)}k = \alpha!\frac {\abs V!}{(k!)^\alpha(\abs V-k\alpha)!}.$$

# 04/21

Formula in report.pdf reduces problem of largest connected component size to problem of counting connected graphs having $n$ vertices and $k$ edges.
A result (recursive form) by Marko Riedel can be found on MSE (#689526) but only works for $n \leq 11$\ldots Let's then try to use the same type of reasonment as he did in order to find a correct formula:

$$\begin{aligned}
	\log\left(1+\sum_{m=1}^n(1+u)^{X(m)}\frac {z^m}{m!}\right) &= \sum_{q \geq 1}(-1)^{q+1}\frac 1q\left[\sum_{m=1}^n(1+u)^{X(m)}\frac {z^m}{m!}\right]^q \\
	&\simeq \sum_{q=1}^n(-1)^{q+1}\frac 1q\sum_{\abs \alpha = q}\binom {q}{\alpha}\prod_{m=1}^n\left((1+u)^{X(m)}\frac {z^m}{m!}\right)^{\alpha_m} \\
	&=: G(z, u).
\end{aligned}$$

Therefore, we find:
$$\begin{aligned}~
	[u^k]G(z, u) &= \sum_{q=1}^n(-1)^{q+1}\frac 1q\sum_{\abs \alpha=q}\binom q\alpha\left(\prod_{m=1}^n\frac {z^{m \cdot \alpha_m}}{(m!)^{\alpha_m}}\right)[u^k](1+u)^{\sum_{m=1}^nX(m)\alpha_m} \\
	&=\sum_{q=1}^n(-1)^{q+1}\frac 1q\sum_{\abs \alpha=q}\binom q\alpha\left(\prod_{m=1}^n\frac {z^{m \cdot \alpha_m}}{(m!)^{\alpha_m}}\right)\binom {\sum_{m=1}^nX(m)\alpha_m =: \beta_\alpha}{k},
\end{aligned}$$

and thus, for $\Theta_q$ lthe set of all vectors $\alpha$ in $\mathbb N^n$ such that $\sum_{m=1}^nm\alpha_m = q$ :

$$\begin{aligned}~
[u^k][z^n]G(z, u) &= \sum_{q=1}^n(-1)^{q+1}\frac 1q\sum_{\abs \alpha=q}\binom q\alpha[z^n]\left(\prod_{m=1}^n\frac {z^{m \cdot \alpha_m}}{(m!)^{\alpha_m}}\right)\binom {\beta_\alpha}{k} \\
&= \sum_{q=1}^n(-1)^{q+1}\frac 1q\sum_{\alpha \in \Theta_q}\binom q\alpha\binom {\beta_\alpha}k\prod_{m=1}^n\frac 1{(m!)^{\alpha_m}}
\end{aligned}$$

# 04/23

Riedel's formula:
$$q_{n,k} = \binom {X(n)}k - \sum_{m=0}^{n-2}\binom {n-1}m\sum_{p=0}^k\binom {X(n-m-1)}pq_{m+1,k-p}$$
is indeed correct, and thus complete formula works. Last step to do is to write about Riedel's solution in mathematical report.


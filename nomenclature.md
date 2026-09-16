# Nomenclature

**Project:** How does correlation relate to causation: seven insights from an analytic relationship between the two

**Last updated:** 2026-09-16

---

## 1. Typographic conventions

| Category | Rule | Examples |
|----------|------|---------|
| Matrices | Bold uppercase | $\mathbf{A}$, $\mathbf{I}$, $\mathbf{\Sigma}_Y$, $\mathbf{\Sigma}_W$, $\mathbf{K}$ |
| Vectors | Bold; lowercase for generic vectors, uppercase for state/process variables | $\mathbf{X}(t)$, $\boldsymbol{\xi}$, $\mathbf{B}_t$ |
| Scalars | Plain italic | $a_{ij}$, $\sigma_{ij}$, $\sigma_{W_i}$, $\rho_{ij}$, $n$, $t$ |
| Sets / spaces | Blackboard bold | $\mathbb{R}$, $\mathbb{R}^{n\times n}$ |
| Operators | Upright roman | $\operatorname{vec}$, $\operatorname{vech}$, $\operatorname{tr}$, $\operatorname{sgn}$, $\det$ |


---

## 2. Operators and decorations

| Symbol | Meaning |
|--------|---------|
| $Cov(\cdot)$ | Covariance |
| $\det(\cdot)$ | Determinant |
| $\delta(\cdot)$ | Dirac delta |
| $E[\cdot]$ | Expectation |
| $\otimes$ | Kronecker product |
| $\oplus$ | Kronecker sum, $\mathbf{A}\oplus\mathbf{A} = \mathbf{I}\otimes\mathbf{A}+\mathbf{A}\otimes\mathbf{I}$ |
| $\mathbb{R}$ | Real numbers |
| $\mathrm{Re}(\cdot)$ | Real part |
| $\operatorname{sgn}(\cdot)$ | Sign function |
| $\operatorname{tr}(\cdot)$ | Trace |
| ${}^T$ | Transpose |
| $\operatorname{vec}(\cdot)$ | Column-stacking vectorization |
| $\operatorname{vech}(\cdot)$ | Half-vectorization (lower triangle + diagonal) |

---

## 3. Symbols

### Scalars

| Symbol | Meaning | Notes |
|--------|---------|-------|
| $n$ | Number of nodes / variables in the network | Integer |
| $d$ | Dimension of the Brownian motion $\mathbf{B}_t$ | Appendix A; general case $d\neq n$, paper uses $d=n$ |
| $m$ | Dimension of the observation vector $\mathbf{Y}(t)$ | Appendix B; general case $m\neq n$, paper uses $m=n$ |
| $t$ | Time | Continuous |
| $a_{ij}$ | Causal influence of variable $j$ on variable $i$ | Off-diagonal entry of $\mathbf{A}$; note index order is $j\to i$, not $i\to j$ |
| $a_{ii}$ | Self-connection (e.g.,self-inhibition) | negative by default for stability |
| $a$ | Shared self-connection value | Used when $a_{11}=a_{22}=a_{33}=a$ in triadic topologies |
| $\lambda_i$ | Eigenvalue of $\mathbf{A}$ | Hurwitz stability: $\mathrm{Re}(\lambda_i) < 0\ \forall i$ |
| $\sigma_{W_i}$ | Variance (not std. dev.) of the noise entering node $i$ | Diagonal entry of $\mathbf{\Sigma}_W$ |
| $\sigma_{ij}$ | Entry $(i,j)$ of the steady-state covariance $\mathbf{\Sigma}_Y$ | $\sigma_{ii}$ = variance of node $i$ |
| $\rho_{ij}$ | Pairwise correlation coefficient | $\rho_{ij} = \sigma_{ij}/\sqrt{\sigma_{ii}\sigma_{jj}}$ |
| $N_{ij}$ | Numerator polynomial of $\rho_{ij}$ | $\rho_{ij} = N_{ij}/\sqrt{D_{ii}D_{jj}}$ |
| $D_{ii}$ | Denominator polynomial of $\rho_{ij}$ (per-node factor) | $\rho_{ij} = N_{ij}/\sqrt{D_{ii}D_{jj}}$ |
| $k_{ij}^{uv}$ | Entry of $\mathbf{K}$ at row $(i,j)$, column $(u,v)$ | Row/column indices flattened via the $n^2\times n^2$ Kronecker layout |
| $C_{ij}^{uv}$ | Cofactor of $\mathbf{K}$ corresponding to $k_{ij}^{uv}$ | Used to expand $\mathbf{K}^{-1}$ |
| $M_{\cdot,\cdot}$ | Minor of $\mathbf{K}$ | Determinant of a submatrix of $\mathbf{K}$ |
| $g_0$ | Polynomial function expressing $\vert\mathbf{K}\vert$ in the entries of $\mathbf{K}$ | |
| $g$ | Polynomial function expressing $\vert\mathbf{K}\vert$ after substituting entries of $\mathbf{A}$ | |
| $\tau$, $\tau_1$, $\tau_2$ | Dummy integration/time variables | Appendix A (stochastic integrals) |
| $\sigma$ | Dummy integration variable, $\sigma=t-\tau$ | Appendix A only; distinct from $\sigma_{ij}$ above |

### Vectors

| Symbol | Meaning | Notes |
|--------|---------|-------|
| $\mathbf{X}(t)$ | State vector of the linear stochastic dynamical system | $\mathbb{R}^n$ |
| $\mathbf{X}(0)$ | Initial condition of the state vector | $\mathbb{R}^n$ |
| $\mathbf{B}_t$ | Standard Brownian motion | $\mathbb{R}^d$ |
| $\boldsymbol{\xi}$ | Diffusion matrix scaling random fluctuations | $\mathbb{R}^{n\times d}$ |
| $\mathbf{Y}(t)$ | Observation vector | $\mathbb{R}^m$ |
| $\mathbf{V}(t)$ | Observation-noise vector | $\mathbb{R}^m$ |

### Matrices

| Symbol | Meaning | Notes |
|--------|---------|-------|
| $\mathbf{A}$ | Causal influence matrix (effective connectivity) | $\mathbb{R}^{n\times n}$; Hurwitz stable |
| $\mathbf{\Sigma}_Y$ | Steady-state covariance matrix | $\mathbb{R}^{n\times n}$; symmetric positive semi-definite |
| $\mathbf{\Sigma}_W$ | Noise covariance matrix | $\mathbb{R}^{n\times n}$; diagonal |
| $\mathbf{\Sigma}_0$ | Covariance of $\mathbf{X}(0)$ | $\mathbb{R}^{n\times n}$; drops out at steady state |
| $\mathbf{K}_{XX}(t_1,t_2)$ | Autocovariance of $\mathbf{X}(t)$ | $\mathbb{R}^{n\times n}$; Appendix A |
| $\mathbf{\Sigma}_{X(t)}$ | Covariance of $\mathbf{X}(t)$ | $\mathbb{R}^{n\times n}$; $\mathbf{K}_{XX}(t_1,t_2)$ at $t_1=t_2=t$ |
| $\mathbf{I}$ | Identity matrix | Size inferred from context |
| $\mathbf{K}$ | $\mathbf{I}\otimes\mathbf{A}+\mathbf{A}\otimes\mathbf{I}$ | $\mathbb{R}^{n^2\times n^2}$; $\vert\mathbf{K}\vert$ is the common denominator of every $\sigma_{ij}$ |
| $\mathbf{C}$ | Cofactor matrix of $\mathbf{K}$ | $\mathbb{R}^{n^2\times n^2}$; $\mathbf{K}^{-1} = \frac{1}{\vert\mathbf{K}\vert}\mathbf{C}^T$ |
| $\mathbf{L}_n$ | Elimination matrix | $\mathbb{R}^{\frac{n(n+1)}{2}\times n^2}$; $\operatorname{vech}(\cdot) = \mathbf{L}_n\operatorname{vec}(\cdot)$ |
| $\mathbf{D}_n$ | Duplication matrix | $\mathbb{R}^{n^2\times\frac{n(n+1)}{2}}$; $\operatorname{vec}(\cdot) = \mathbf{D}_n\operatorname{vech}(\cdot)$ |
| $\mathbf{H}$ | Observation/output matrix | $\mathbb{R}^{m\times n}$; $\mathbf{\Sigma}_Y = \mathbf{H}\mathbf{\Sigma}_\infty\mathbf{H}^T + \mathbf{\Sigma}_V$ |
| $\mathbf{\Sigma}_\infty$ | Steady-state covariance before observation | $\mathbb{R}^{n\times n}$ |
| $\mathbf{\Sigma}_V$ | Observation-noise covariance | $\mathbb{R}^{m\times m}$; $\mathbf{0}$ in the main paper |


---

## 4. Abbreviations and acronyms

| Abbreviation | Full term | Notes |
|-------------|----------|-------|
| LSDS | Linear stochastic dynamical system | The model class studied throughout |
| SDE | Stochastic differential equation | Itô form used in Appendix A |
| OU | Ornstein--Uhlenbeck (process) | $\mathbf{X}(t)$ is a multivariate OU process |
| VAR | Vector autoregressive (model) | Discrete-time analog of the OU process |
| PSD | Positive semi-definite | Property of $\mathbf{\Sigma}_W$ |
| SCM | Structural causal models | Appendix E, prior-work comparison |
| DCM | Dynamic causal modeling | Appendix E, prior-work comparison |
| GCLM | Graphical continuous Lyapunov model | Appendix E, prior-work comparison |

---


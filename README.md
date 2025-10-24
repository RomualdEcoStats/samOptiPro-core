# samOptiPro 0.1.0

**samOptiPro** — *Systematic Optimization of Bayesian Stock Assessment Models in NIMBLE*

A workflow-driven helper to **configure, assess, and optimise** MCMC sampling in [NIMBLE](https://r-nimble.org/), following a reproducible decision tree:

> general tools → assess → detect poor performance → identify bottlenecks →  
> (model surgery | custom samplers) → reassess → validate → iterate.

---

## ✳️ Overview

**samOptiPro**  provide an advanced, modular workflow for **diagnosing**, **benchmarking**, and **optimizing** hierarchical ecological models (SAM-like frameworks) built in NIMBLE.

It integrates:
- **Structural diagnostics**: detect non-differentiable nodes (truncations, Dirichlet, simplex constraints…)
- **Adaptive sampler configuration**: auto-assign Slice / AF_slice / RW / Block / HMC / NUTS
- **Performance analytics**: algorithmic (ESS, ESS/s) and computational (runtime) efficiency
- **Automatic visual reports**: bottlenecks, convergence, and Rhat distributions
- **Differentiability testing**: seamless handoff to gradient-based inference (nimbleHMC)
- **Hybrid family strategy**: use block samplers or adaptive HMC where correlations demand it

---

## ⚙️ Installation

```r
# From local development directory
devtools::load_all("samOptiPro")

# Or (once hosted)
# remotes::install_github("RomualdEcoStats/samOptiPro-core")


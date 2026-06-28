# When More Data Make the Wrong Model More Certain
## Network Dependence, Regularization, and Mean-Field Bias in Agent-Based SIR Models

---

## Abstract

Agent-based models are powerful because they can represent heterogeneous agents,
local interaction, network structure, and path dependence. They are controversial
for the same reason: the flexibility that makes them expressive can also make them
difficult to validate. A sufficiently flexible simulation may reproduce a pattern
without making clear which mechanisms are responsible, how the simulation connects
to data, or how its conclusions would change under alternative assumptions.

This paper proposes a statistical framing of agent-based models as generative
models. In this view, a simulator produces a latent social process, while observed
data arise through a separate observation process operating at a chosen level of
aggregation and time scale. The distinction is important: many social data sets do
not observe the process of interest directly. They record institutional traces,
aggregate counts, sampled events, or temporally coarsened measurements. A
simulation is therefore not yet a statistical model until it specifies how latent
simulated states become observed data.

---

## 1. The Problem

Agent-based modelers routinely simulate epidemic spread, opinion dynamics, and
market behavior at the level of individual agents. The simulation produces rich
internal trajectories: who interacted with whom, in what order, under what
conditions. But the data available for calibration are almost always coarser than
the simulation output. A modeler fitting an epidemic model to weekly case reports
is comparing a fine-grained stochastic process to aggregate counts measured at a
weekly clock.

This mismatch creates a problem that is easy to overlook. When the model is
calibrated to coarse data, the estimated parameters may be fitting the observation
process as much as the underlying mechanism. If the observation process discards
information that the mechanism depends on, the fitted parameters will not recover
the true microscopic rates---no matter how much data are collected.

---

## 2. Contribution

We reframe agent-based models as generative statistical models with two
components:

$$
\text{simulator} \to \text{latent social process} \to \text{observation process} \to \text{data}
$$

The simulator encodes a theory of how agents interact. The observation process
encodes how those interactions are measured---at what spatial and temporal
resolution, with what reporting rate, and with what aggregation structure. A
simulation is not yet a statistical model until both components are specified.

We demonstrate the consequences of fitting the wrong model to coarse data using
epidemic SIR models as a running example.

---

## 3. The SIR Example and the Social Assumption in the Product Term

The standard mean-field ordinary differential equation (ODE) SIR model is:

$$
\frac{dS}{dt} = -\beta_{ODE} \frac{SI}{N}, \quad
\frac{dI}{dt} = \beta_{ODE} \frac{SI}{N} - \gamma I, \quad
\frac{dR}{dt} = \gamma I.
$$

The term $SI$ is often treated as algebra. It is not. It encodes a strong social
assumption: that the rate at which susceptible individuals encounter infected
individuals is proportional to the product of their aggregate counts. This is the
uniform-mixing assumption. Every susceptible agent is equally likely to contact
any infected agent regardless of social structure, geography, or group membership.

The product term is a social assumption.

---

## 4. Network SIR and the Dependence on SI Edges

A network SIR model replaces the uniform-mixing assumption with an explicit
contact structure. Each agent $i$ interacts only with its neighbors in a contact
network $A$. The infection probability for susceptible agent $i$ over a short
interval $\Delta t$ is:

$$
P\bigl(i \text{ infected in } [t, t+\Delta t]\bigr) =
1 - \exp\!\Bigl(-\beta \,\Delta t \sum_j A_{ij} \mathbf{1}\{X_j(t) = I\}\Bigr),
$$

where $\beta$ is the per-contact, per-infected-neighbor transmission rate. This is
not the same parameter as $\beta_{ODE}$.

In the network model, aggregate incidence depends on $E_{SI}(t)$, the number of
edges connecting a susceptible node to an infected node:

$$
\text{incidence}_{network}(t) = \beta \, E_{SI}(t).
$$

The mean-field ODE approximates this as:

$$
E_{SI}(t) \approx c \, S(t) I(t),
$$

where $c$ is a constant determined by the contact structure. Under uniform mixing
with mean degree $\bar{d}$, this approximation holds with $c = \bar{d}/N$, giving
an effective mean-field rate:

$$
\beta_{mf,true} = \beta \, \bar{d}.
$$

When the contact network has community structure---dense within-group ties,
sparse between-group ties---the approximation $E_{SI} \approx c \, SI$ fails.
The number of active SI edges depends on where the epidemic front sits relative
to group boundaries, not just on aggregate $S$ and $I$ counts.

---

## 5. Regularization and the Uniform Graph as the Maximally Regularized Case

We formalize the relationship between network models and mean-field models through
regularization. Regularization penalizes model flexibility that is not supported
by the data. In the network setting, we generate graphs from a parametric family:

$$
\text{logit}\, P(A_{ij} = 1) = \alpha + \eta \, B_{ij},
$$

where $B_{ij} = 1$ if agents $i$ and $j$ belong to the same group and
$B_{ij} = -1/(K-1)$ otherwise. The parameter $\eta$ controls community structure
strength; $\alpha$ controls overall edge density. The regularization penalty is
$\lambda \eta^2$.

When $\eta = 0$, every dyad has equal probability $\sigma(\alpha)$ of an edge.
The graph is Erdős–Rényi and the epidemic mixes uniformly.

The uniform graph is not the neutral case. It is the maximally regularized case:
the network model that remains when all distinctions based on group membership
have been shrunk to zero.

The standard ODE is correspondingly the maximally regularized aggregate model. It
is not what happens when an ABM gets large. It is what happens when the ABM gets
large, coarse, and socially uniform.

---

## 6. Posterior Bias under Community Structure

We demonstrate the bias that results from fitting the mean-field ODE to data
generated by a community network SIR process. The simulation uses $N = 1000$
agents divided into $K = 4$ equal groups. The contact network is a stochastic
block model with $\alpha = -4.4$ and $\eta = 2.5$, yielding a mean degree of
approximately 37 within a highly clustered structure (within-group edge share
0.89). The true per-contact transmission rate is $\beta = 0.04$ per day and the
recovery rate is $\gamma = 0.10$ per day (mean infectious period 10 days).

The mean-field effective rate implied by this network under uniform mixing is:

$$
\beta_{mf,true} = \beta \cdot \bar{d} \approx 0.04 \times 36.7 \approx 1.47.
$$

An unbiased ODE calibration should recover $\beta_{ODE}$ near this value.

We generate a single 150-day epidemic, aggregate infections into 21 weekly counts,
and fit the ODE via a grid posterior over $\beta_{ODE}$ with $\gamma$ fixed at its
true value. The ODE likelihood uses a log-count Gaussian:

$$
\log(1 + Y_k) \sim \mathcal{N}\bigl(\log(1 + \hat{Y}_k),\, \sigma^2\bigr),
$$

where $Y_k$ is observed incidence in week $k$ and $\hat{Y}_k$ is the ODE
prediction.

The posterior concentrates near $\beta_{ODE}^* \approx 1.13$, well below the
reference $\beta_{mf,true} \approx 1.47$ (bias $\approx -0.33$, or 23%). For
comparison, the same procedure applied to data from a uniform-mixing network
recovers $\beta_{ODE}$ within 2% of its reference value.

Aggregation can make the ODE look right while making its parameters wrong. The
fitted epidemic curve tracks the observed weekly counts acceptably; the posterior
bias is invisible in the trajectory plot.

---

## 7. More Data Make the Wrong Model More Certain

The central result concerns what happens when we collect more data at the same
observational grain. We generate $R$ independent community-network outbreaks under
identical parameters and form a joint likelihood by multiplying across outbreaks:

$$
\log L(\beta_{ODE}) = \sum_{r=1}^{R} \log L_r(\beta_{ODE}).
$$

As $R$ increases from 1 to 20, the posterior 95% credible interval width shrinks
from 0.136 to 0.036---a fourfold reduction. The posterior MAP remains in the range
$[1.11, 1.22]$, always well below the reference $\beta_{mf,true} \approx 1.47$.

More data at the same observational grain do not recover the missing mechanism.
They estimate the wrong reduced model more precisely. Variance shrinks. Bias
remains.

---

## 8. The Regularization Path

The connection between community structure and ODE bias is monotone. We vary
$\eta$ from 0 to 3, pooling 20 independent outbreaks per level to obtain stable
pseudo-true estimates, ODE bias increases strictly monotonically: −0.044 (9% of
reference) at $\eta = 0$ and −0.48 (23% of reference) at $\eta = 3$. The
relationship is convex: bias is modest for small $\eta$ and accelerates as
community structure becomes strong. As community structure strengthens, the
epidemic increasingly spreads in group-wise waves; the ODE compensates by fitting
a lower effective transmission rate.

This is the regularization path made visible: as $\eta \to 0$, the community
network collapses toward the uniform graph, and ODE bias toward zero.

---

## 9. Robustness, Priors, and Perturbation Stability

This argument connects to PAC-Bayesian ideas, which formalize the tradeoff between
empirical fit and model flexibility. A model must pay for how far it moves from
prior plausible mechanisms in order to fit the data. In the network setting, prior
assumptions define a distribution over plausible contact structures. Fitted models
that require very specific, complex network configurations to explain the data are
penalized relative to simpler alternatives.

More generally, regularization can be understood as a requirement that conclusions
remain stable under modeling changes that the theory regards as scientifically
irrelevant. If a finding depends critically on exact network realization, specific
event ordering, or fine temporal resolution, it is not a robust finding. Stable
findings are those that survive perturbations within the set of admitted models.

This perspective motivates a practical workflow for model criticism:

- Compare richer simulations against their regularized mean-field limits.
- Vary the observation grain (daily, weekly, monthly) and examine how parameter
  estimates change.
- Vary the community structure parameter $\eta$ and track posterior bias.
- Run posterior predictive checks: does the fitted ODE reproduce features of the
  data that were not used in calibration?
- Test robustness across network realizations drawn from the same generator.

---

## 10. Conclusion

We have shown that fitting a mean-field ODE SIR model to data from a
community-structured network epidemic produces systematically biased parameter
estimates, and that collecting more data at the same coarse aggregation level
sharpens the bias rather than correcting it. The ODE posterior concentrates around
a pseudo-true parameter that absorbs the effects of community structure, temporal
aggregation, and observation grain into a single effective rate. The analyst
becomes more confident---but about the wrong model.

The broader goal is to make social simulation empirically disciplined: not merely
to generate plausible patterns, but to determine which mechanisms are needed, which
can be regularized away, and where a model's claims cease to be robust.

---

*Repository:* All simulation code, inference scripts, and figures are available
at [this repository]. Key outputs:

- `figures/fig_03_posterior_bias_single.png`: posterior bias for a single community-network outbreak.
- `figures/fig_04_more_data_shrinks_variance_not_bias.png`: posterior narrowing without bias reduction across $R = 1, 5, 10, 20$ outbreaks.
- `figures/fig_05_regularization_path_network_stats.png`: monotone increase in ODE bias with community structure strength $\eta$.

*Word count: approximately 1500 words (body text).*

# Slides Brief: "When More Data Make the Wrong Model More Certain"
## CSSSA 2025 — Reveal.js Presentation

---

## Context

This is a 20-minute academic conference presentation for CSSSA 2025. The paper demonstrates
that fitting a mean-field ODE SIR model to coarse aggregate data from a community-structured
network epidemic produces biased posteriors that **sharpen (not correct) with more data**.

The audience knows SIR models and agent-based modelling. They do not necessarily know the
distinction between network and ODE transmission rates.

**Target:** 15 slides. Clean, one idea per slide.

---

## Technical Setup

- **Framework:** Reveal.js 4.6.1 from CDN
- **Math:** KaTeX via the Reveal.js math plugin (CDN). Use `\(...\)` for inline, `\[...\]` for display math. Load the plugin properly:
  ```js
  plugins: [ RevealMath.KaTeX ],
  math: {
    katexScript: 'https://cdn.jsdelivr.net/npm/katex@0.16.9/dist/katex.min.js',
    katexStylesheet: 'https://cdn.jsdelivr.net/npm/katex@0.16.9/dist/katex.min.css',
  }
  ```
- **Slide size:** width 960, height 620, margin 0.06
- **Theme:** white
- **Transitions:** slide, fast
- Hash, slideNumber, controls, progress all on

### Critical CSS warning

Reveal.js computes its base font-size dynamically (roughly 40px for a 960×620 slide).
Any `em` units inside slides inherit from this large base. For structural layout elements
(boxes, diagrams, fixed-size containers), **always use `px` not `em`**, or the text will
overflow. For ordinary prose and list items, `em` is fine.

Display math equations (`\[...\]`) can easily be too wide for the slide. Keep display
equations short, or break them onto multiple lines using `\\` and `\quad` spacing rather
than putting all three ODE equations on one line.

---

## Available Figures

All figures are at `../figures/` relative to `slides/index.html`.

| File | What it shows |
|------|---------------|
| `fig_02_network_vs_ode_incidence.png` | Community-network weekly incidence (bars) vs. ODE fit (red line) vs. ODE at β_mf,true (dashed). Shows the two-phase shape the ODE cannot match. |
| `fig_03_posterior_bias_single.png` | ODE posterior under uniform data (blue, near reference) vs. community data (red, 23% below reference). Both tight. |
| `fig_04_slide.png` | Single-panel overlaid posteriors for R=1,5,10,20 outbreaks (faint→dark red). CI width annotated above each peak. Peak does not move toward reference. |
| `fig_05_regularization_path_network_stats.png` | Within-group edge share (green, left axis) and ODE bias (red, right axis) both increase monotonically with η. |
| `fig_beta_networks.png` | Two small network cartoons side by side. Left: Erdős–Rényi (circular layout, one colour, random edges). Right: SBM (4 colour-coded clusters, dense within-group edges dark, sparse between-group edges light). β annotations below each panel. |

---

## Key Parameters (use these exact values in annotations)

| Quantity | Value | What it is |
|----------|-------|-----------|
| N | 1000 | Agents |
| K | 4 | Equal groups |
| β_true | 0.04 day⁻¹ | Per-contact, per-infected-neighbor transmission rate (network ABM) |
| γ | 0.10 day⁻¹ | Recovery rate (mean infectious period 10 days) |
| d̄ (community network) | ≈ 37 | Mean degree |
| β_mf,true | ≈ 1.47 | β_true × d̄ — counterfactual reference (what unbiased ODE recovers from ER network of same density) |
| β_ODE\* | ≈ 1.13 | Where the posterior concentrates under community data (pseudo-true) |
| bias | −0.33 (−23%) | β_ODE\* − β_mf,true |
| η (community) | 2.5 | SBM community structure parameter |
| α | −4.4 | SBM baseline log-odds |

### The three β levels — critical for pedagogy

These are not the same quantity and cannot be compared directly:

1. **β_true = 0.04** — per-contact rate in the network ABM
2. **β_mf,true ≈ 1.47** — β_true × d̄; what an unbiased ODE would recover from an ER network of the same mean degree. This is a *counterfactual reference*, not a "true" ODE parameter.
3. **β_ODE\* ≈ 1.13** — where the ODE posterior actually concentrates when fitted to community-network data. This is 23% below the reference.

The jump from 0.04 to 1.47 is not a bias — it is a change of scale (per-contact → mass-action). The bias is the gap between 1.47 and 1.13.

---

## Slide Structure

### Slide 1 — Title

**When More Data Make the Wrong Model More Certain**
Network Dependence, Regularization, and Mean-Field Bias in Agent-Based SIR Models
CSSSA 2025 · [Author]

---

### Slide 2 — The SIR Model (ODE)

**Heading:** The SIR Model

Three compartments: Susceptible → Infected → Recovered.

Display equations (keep on separate lines, not side by side):
```
dS/dt = −β_ODE · SI/N
dI/dt =  β_ODE · SI/N − γ I
dR/dt =  γ I
```

Parameters:
- β_ODE — transmission rate (mass-action)
- γ — recovery rate; mean infectious period = 1/γ
- N = S + I + R (constant)

Brief note: standard, tractable, widely used — but contains a hidden assumption.

---

### Slide 3 — The Social Assumption

**Heading:** The Hidden Social Assumption

Focus on the infection term. Zoom in on `β_ODE · SI/N`.

The product `S · I` encodes **uniform mixing**:
- Every susceptible is equally likely to contact *any* infected agent
- Social structure, geography, group membership — all absent
- Contact opportunities scale with the product of aggregate counts

**Key statement (make visually prominent):**
> The product S · I is a social assumption, not algebra.
> It is wrong wherever people mostly interact within social groups.

---

### Slide 4 — Network SIR

**Heading:** Network SIR: Replacing the Product Term

Replace `SI/N` with an explicit contact network A. Agent i can only be infected by its neighbours N(i).

Infection probability for agent i over interval Δt:
```
P(i infected in [t, t+Δt]) = 1 − exp(−β · Δt · Σ_{j∈N(i)} 1{X_j = I})
```
(Derived from competing hazards of independent Poisson contact processes — Kiss, Miller & Simon 2017, Ch. 3.)

Key changes:
- β is now **per-contact, per-infected-neighbor** (different scale from β_ODE)
- Aggregate incidence depends on **E_SI(t)**: the number of susceptible–infected *edges*, not the product S·I
- Under community structure, E_SI(t) does not track S·I — it depends on where the epidemic front sits relative to group boundaries

---

### Slide 5 — The Two β Scales (figure)

**Heading:** Same Density, Different Structure — Different ODE Estimate

One line of text above the figure:
> Under uniform mixing with mean degree d̄: β_mf,true = β_true × d̄ ≈ 0.04 × 37 ≈ 1.47

Then show **fig_beta_networks.png** filling most of the slide.

The figure already contains the β annotations below each network panel. No additional annotation needed on this slide.

---

### Slide 6 — Stochastic Block Model

**Heading:** Generating Community Structure: The Stochastic Block Model

Edge probability formula:
```
logit P(A_ij = 1) = α + η · B_ij
```
where B_ij = +1 if i and j share a group, B_ij = −1/(K−1) otherwise.

- **η = 0**: all dyads equal → Erdős–Rényi, uniform mixing
- **η = 2.5**: strong community structure (within-group edge share ≈ 0.89)
- α = −4.4 controls overall density (mean degree ≈ 37)

**Key point:** η = 0 is not the neutral case. It is the *maximally regularised* case — the network that remains when all group distinctions have been shrunk to zero.

---

### Slide 7 — Two-Phase Dynamics (figure)

**Heading:** Community Structure Creates Two-Phase Epidemic Dynamics

Show **fig_02_network_vs_ode_incidence.png**.

Brief caption:
- Bars: observed weekly incidence from community-network simulation
- Red line: ODE fitted to this data — compromises between the two phases
- Dashed: ODE at β_mf,true — peaks too early

**Key point:** No single β_ODE can reproduce the two-phase shape. The fitted ODE lands on a pseudo-true β_ODE\* that minimises aggregate misfit — not the true transmission rate.

---

### Slide 8 — Calibration Setup

**Heading:** Calibration: Grid Posterior over β_ODE

How we fit the ODE to data:

Observation model (log-count Gaussian):
```
log(1 + Y_k) ~ N(log(1 + Ŷ_k), σ²)
```
where Y_k is observed incidence in week k, Ŷ_k is the ODE prediction.

Setup:
- β_ODE grid: [0.01, 4.0], 400 points, log-spaced
- γ fixed at true value (0.10 day⁻¹)
- σ = 0.3
- Posterior normalised via log-sum-exp

This is a minimal, transparent inference setup. The goal is to show the bias is structural — not an artefact of a particular inference method.

---

### Slide 9 — Result 1: Posterior Bias (figure)

**Heading:** Result 1: The ODE Posterior Is Systematically Biased

Show **fig_03_posterior_bias_single.png**.

Caption:
- Blue: posterior fitted to uniform-network data — correctly near the reference (dashed)
- Red: posterior fitted to community-network data — concentrated 23% below reference
- Both posteriors are tight. The bias is misspecification, not uncertainty.

**Key point:** The fitted epidemic curve looks acceptable. The bias is invisible in the trajectory plot — it only appears when you compare the posterior to the counterfactual reference.

---

### Slide 10 — Result 2: More Data (figure)

**Heading:** Result 2: More Data at the Same Grain

Show **fig_04_slide.png**.

Caption: R = 1, 5, 10, 20 independent community-network outbreaks (faint → dark red). 95% CI width narrows fourfold. The peak does not move toward the reference (dashed).

**Key statement (make visually prominent):**
> Variance shrinks. Bias remains.
> More data make the wrong model more certain.

---

### Slide 11 — Result 3: Regularisation Path (figure)

**Heading:** Result 3: Bias Grows Monotonically with Community Structure

Show **fig_05_regularization_path_network_stats.png**.

Caption: Varying η from 0 to 3. Within-group edge share (green, left axis) and ODE bias (red, right axis) both increase. At η = 0 (ER), bias ≈ 0. At η = 3, bias = −0.48.

---

### Slide 12 — The Regularisation Interpretation

**Heading:** A Regularisation Interpretation

The ODE is the maximally regularised limit of the network model:

- **η → 0**: community network collapses to ER; ODE bias → 0
- **η > 0**: group distinctions exist; ODE absorbs them into a lower effective rate
- **η large**: strong community structure; ODE bias is large and convex

The regularisation penalty is λη². At λ → ∞, the only surviving model is the uniform-mixing ODE.

This connects to PAC-Bayesian ideas (McAllester 1999): a model pays for moving far from prior plausible mechanisms to fit data. The ODE is what you get when community structure is penalised to zero.

---

### Slide 13 — What This Means for Practice

**Heading:** Implications for Calibration Practice

- **Specify the observation process explicitly** before fitting. What resolution? What aggregation? What reporting rate?
- **Compare richer simulations to their mean-field limits.** What can the ODE not represent?
- **Vary the observation grain** (daily/weekly/monthly) and track how estimates change.
- **Examine posterior bias**, not just posterior width. A tight posterior around the wrong value is worse than a wide honest one.
- **Run posterior predictive checks.** Does the fitted ODE reproduce features not used in calibration?

---

### Slide 14 — Conclusion

**Heading:** Conclusion

- The `S·I` product term encodes uniform mixing — a strong social assumption that fails under community structure
- Fitting the ODE to community-network data produces a biased posterior that absorbs group structure into a lower effective rate
- More data at the same coarse grain sharpens the bias — variance shrinks, bias remains
- ODE bias grows monotonically with η: stronger community structure → worse misspecification
- The uniform-mixing ODE is the maximally regularised limit of the network model — not a neutral baseline

**Closing statement:**
> The goal is not to generate plausible patterns. It is to determine which mechanisms are needed, which can be regularised away, and where a model's claims cease to be robust.

References: Kermack & McKendrick (1927) · Holland et al. (1983) · Newman (2002) · Kiss et al. (2017) · McAllester (1999)

---

## Style Notes

- **One idea per slide.** If a slide needs more than 4 bullet points, split it.
- **Equations:** never put three ODE fractions side by side on one line — they overflow. Stack them vertically or use a two-line display.
- **Figures:** `max-height: 400px`, centred, with a brief italic caption below. Do not over-annotate in the slide text if the figure already carries the information.
- **Callout boxes:** use sparingly — one per slide at most, for the key takeaway only.
- **Workflow diagram:** omitted per author's preference. The narrative carries the structure implicitly.
- **Colour conventions** (matching the figures):
  - Blue `rgb(38, 115, 191)`: uniform network, unbiased posterior
  - Red `rgb(204, 51, 26)`: community network, biased posterior
  - Green `rgb(51, 153, 77)`: network-side quantities (β_true, E_SI)
  - Black dashed line: β_mf,true reference

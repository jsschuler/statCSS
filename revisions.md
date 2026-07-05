# Revision Suggestions: Network SIR, Observation Grain, Clock, and PAC-Bayesian Regularization

## Core Reframing

The paper should be structured around a sharper inferential question:

> Networks matter for pandemic dynamics. Therefore, the mean-field SIR model is strictly misspecified as a mechanistic model. But the key empirical question is not merely whether the mean-field assumption is false. It is whether the available data, given their observational grain and clock, contain enough information to recover the relevant network structure.

This reframing avoids the stale “ODE bad, ABM good” argument. The stronger contribution is about **recoverability under observation constraints**.

The mature version of the argument is:

> More data make the wrong model more certain when the observation process is too coarse to identify the mechanism. Richer grain or finer clock data can recover network structure; coarse data require regularized uncertainty rather than false mechanistic confidence.

## Proposed Paper Spine

1. **Assume networks matter for pandemics.**

   The CSSSA audience is unlikely to challenge this. The paper does not need to spend much energy proving that pandemic transmission is network-dependent.

2. **Mean-field SIR is strictly incorrect as a mechanistic model.**

   The `SI/N` term encodes uniform mixing. It is not just algebra. It is a social assumption about how susceptible and infected people encounter each other.

3. **The real question is empirical recoverability.**

   Even if network structure matters in the latent process, the observation process may or may not allow us to discover it.

4. **Observation grain and clock determine what can be learned.**

   Fine-grained data may recover network structure. Coarse aggregate data may only support a reduced mean-field approximation.

5. **PAC-Bayesian regularization provides the inferential framework.**

   The ABM must allow heterogeneity, but regularization determines whether the data support that heterogeneity.

6. **Under ideal conditions, recover the correct network structure.**

   With sufficiently rich data, PAC-Bayesian inference should recover the correct network or the correct network-generating structure.

7. **Under non-ideal conditions, quantify uncertainty.**

   When data are too coarse, the method should avoid hallucinating structure and instead report uncertainty or support only reduced aggregate parameters.

## Key Conceptual Sentence

Use something like this early in the paper:

> Because network dependence is assumed to be mechanistically relevant, the central inferential problem is not whether the mean-field model is exactly true. It is not. The problem is whether the observation process is rich enough to distinguish network-dependent transmission from a lower-dimensional mean-field approximation.

This should appear near the introduction or contribution section. It prevents readers from mistaking the paper for a generic critique of ODE models.

## Define Grain and Clock Explicitly

The draft should distinguish observational grain from observational clock more sharply.

### Observational Grain

Grain refers to **what level of the latent process is observed**.

Examples:

- Individual infection histories
- Individual infection times
- Group-level incidence
- Network-edge exposure histories
- Aggregate population incidence

### Observational Clock

Clock refers to **when observations are recorded or reported**.

Examples:

- Continuous time
- Event time
- Daily reporting
- Weekly reporting
- Monthly reporting
- Irregular reporting intervals

Suggested paragraph:

> We distinguish observational grain from observational clock. Grain refers to the level at which the latent process is measured: individual infection events, group-level incidence, network-edge exposure histories, or aggregate population counts. Clock refers to the temporal resolution at which those measurements are reported: daily, weekly, event-time, or irregular reporting intervals. In the baseline experiment, the observation process maps a continuous-time network epidemic into weekly aggregate incidence counts. This discards both network position and within-week timing.

## Reframe the Simulation Design

The current weekly-aggregate simulation is useful, but it should become one condition in a broader simulation grid.

Suggested simulation grid:

| Grain | Clock | Expected Result |
|---|---|---|
| Individual infection histories | Event time / daily | Network structure most recoverable |
| Group-level incidence | Daily | Community structure recoverable, exact edges uncertain |
| Group-level incidence | Weekly | Community structure partially recoverable |
| Aggregate incidence | Daily | Some timing clues, weak structural recovery |
| Aggregate incidence | Weekly | Network structure mostly unrecoverable; ODE pseudo-true fit dominates |

The paper’s core empirical object should be a **recoverability map**, not only a bias demonstration.

## Reframe the Existing Result

The existing result is still valuable, but it should be described as the coarse-observation case:

> Under weekly aggregate incidence, additional outbreaks reduce posterior variance around a pseudo-true ODE parameter but do not recover the network mechanism. More data at the same grain and clock estimate the wrong reduced model more precisely.

This is stronger than simply saying “more data make the wrong model more certain.”

The fully qualified claim is:

> More data make the wrong model more certain when the additional data are collected under the same insufficient observation process.

## Clarify the Role of Regularization

The paper should not sound like regularization removes heterogeneity because heterogeneity is undesirable.

The better framing:

> The ABM must be capable of representing heterogeneous contact structure. Regularization does not remove that capacity; it governs whether the data support using it.

Suggested paragraph:

> The role of regularization is not to eliminate heterogeneity from the ABM. Heterogeneity is part of the model’s theoretical content: agents differ in contact neighborhoods, and epidemic transmission depends on the resulting susceptible-infected edges. Rather, regularization determines which forms of heterogeneity are supported by the data. In the stochastic block model, η parameterizes group-based network structure, while the penalty λη² shrinks unsupported group structure toward the uniform graph. Thus the uniform graph is not treated as substantively true; it is the limiting case obtained when the data do not justify retaining group-level contact heterogeneity.

## PAC-Bayesian Role

PAC-Bayes should not appear as a decorative citation near the end. It should become the inferential logic of the paper.

Suggested framing:

> In a PAC-Bayesian formulation, candidate network structures or network-generating parameters are evaluated by predictive performance penalized by departure from a prior over simpler or more regular mechanisms. Under rich observation, such as individual infection times or sufficiently fine group-level incidence, the posterior should concentrate near the correct network structure or community parameter. Under coarse observation, such as weekly aggregate incidence, the posterior should remain diffuse over network structures or concentrate only on reduced effective parameters.

The key point:

> PAC-Bayes lets the model recover structure when the data support it and quantify uncertainty when they do not.

## Strengthen the Contribution Section

Suggested replacement or revision:

> We study when network structure in an agent-based epidemic model is empirically recoverable from observed epidemic data. We assume that pandemic transmission is network-dependent, so the mean-field SIR model is strictly misspecified as a mechanistic account. However, misspecification alone does not determine whether the network mechanism can be learned. Recoverability depends on the observation process: the grain at which infections are measured and the clock at which they are reported.
>
> We use simulated network SIR epidemics to vary both observational grain and observational clock. Under rich observations, network structure can be recovered using regularized inference. Under coarse aggregate observations, the same latent process supports only a reduced mean-field approximation, and additional data at the same grain make that approximation more certain without recovering the missing network mechanism. We connect this to PAC-Bayesian regularization: the ABM must allow heterogeneity, but the data determine whether that heterogeneity is supported or should remain uncertain.

## Be Careful With the Reference Parameter

The paper currently uses:

> β_mf,true = β × d̄

This is useful, but it should be described as a benchmark, not a literal ground truth.

Suggested wording:

> β_mf,true is not a microscopic ground truth. It is a calibration benchmark: the ODE transmission rate recovered when the same per-contact process is run on a uniform network with the same mean degree and observed through the same reporting scheme.

The phrase “same reporting scheme” matters because otherwise readers may suspect that the benchmark mixes mechanism contrast with observation contrast.

## Separate the Sources of Failure

Avoid blending together:

1. Structural misspecification: the ODE assumes uniform mixing.
2. Grain loss: aggregate counts discard group and network position.
3. Clock loss: weekly reporting discards within-week timing.
4. Parameter translation: β in the network model is not β_ODE.

Suggested sentence:

> The bias arises from fitting a structurally misspecified mean-field model to data whose observation process has discarded the network information needed to distinguish community transmission from a lower effective aggregate rate.

## Improve the Figures

### Figure 1

Figure 1 should be made more self-explanatory.

Recommended changes:

- Label the uniform-network posterior directly.
- Label the community-network posterior directly.
- Color-match reference lines to the relevant posterior.
- Make clear why the red posterior is biased rather than merely different.
- Do not force the caption to rescue the figure from visual ambiguity.

### Figure 2

Figure 2 is strong. It directly supports the claim:

> Variance shrinks; bias remains.

Keep it, but make sure the caption emphasizes that all additional data are collected under the same observational grain and clock.

### Figure 3

Avoid saying “bias decreases monotonically” if the substantive claim is that bias worsens.

Use:

> Absolute bias increases.

or:

> Bias becomes increasingly negative.

Do not make readers solve a sign convention puzzle. They are already reading about stochastic block models. This is enough suffering.

## Stronger Posterior Predictive Point

The draft should emphasize that acceptable aggregate fit can coexist with mechanistic failure.

Suggested wording:

> Standard posterior concentration and acceptable aggregate fit can coexist with mechanistic failure. The fitted ODE may track weekly incidence while its parameter absorbs unobserved network structure into a biased effective rate.

This is one of the strongest points in the paper.

## Stronger Conclusion

Suggested language:

> Mean-field SIR is not wrong because it is mathematically crude; it is wrong because it encodes a social theory of uniform contact. Whether that wrongness matters empirically depends on the observation process. Fine-grained data can support network recovery. Coarse data may only support a reduced aggregate model. PAC-Bayesian regularization provides a way to distinguish these cases: retain heterogeneity when the data support it, penalize it when they do not, and quantify uncertainty when the observation process cannot decide.

Even sharper:

> A mechanism cannot be validated at a resolution where its identifying information has been discarded.

And:

> More data help only when they add information relevant to the mechanism, not merely more replications of the same reduced observation.

## Final Diagnosis

The paper should move from:

> More data make the wrong model more certain.

to:

> More data make the wrong model more certain when the observation process is too coarse to identify the mechanism. Richer grain or finer clock data can recover network structure; coarse data require regularized uncertainty rather than false mechanistic confidence.

That is the mature version of the argument: less slogan, more scalpel.
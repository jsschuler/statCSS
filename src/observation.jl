# observation.jl
# Coarsening of latent simulation output into observed aggregate data.
#
# The latent simulation tracks individual agent states at fine time resolution.
# Observed data are coarse aggregate incidence counts at weekly (or other) intervals.
# This coarsening hides event ordering and network structure.
#
# This file makes the observation process explicit and separate from the latent process.
# A simulation is not a statistical model until the observation process is specified.

module Observation

using DataFrames
using Distributions
using Random
using ..SIRTypes: ObservationParams

export observe_coarse_incidence, observe_group_incidence

"""
    observe_coarse_incidence(latent_df; obs_params, rng) -> DataFrame

Aggregate fine-grained latent simulation output into coarse observed counts.

latent_df must have columns: t, new_infections

obs_params controls:
  dt_obs        : observation interval (e.g. 7.0 for weekly)
  reporting_rate: fraction of true cases observed
  noise_sd      : log-normal noise std dev (0.0 = no noise)

Returns DataFrame with columns:
  period        : integer period index (1-indexed)
  t_start       : start time of observation window
  t_end         : end time of observation window
  true_incidence : total new infections in window (from latent simulation)
  reported_incidence : observed count (after reporting process)

When reporting_rate = 1.0 and noise_sd = 0.0 (perfect observation), then
  reported_incidence = true_incidence exactly.
This is the default for the bias demonstration, to isolate misspecification
from observation noise.
"""
function observe_coarse_incidence(
    latent_df::DataFrame;
    obs_params::ObservationParams,
    rng::AbstractRNG
)::DataFrame
    dt_obs = obs_params.dt_obs
    rho = obs_params.reporting_rate
    sigma = obs_params.noise_sd

    t_max = maximum(latent_df.t)
    n_periods = floor(Int, t_max / dt_obs)

    periods = Int[]
    t_starts = Float64[]
    t_ends = Float64[]
    true_incidence = Int[]
    reported_incidence = Float64[]

    for k in 1:n_periods
        t_start = (k - 1) * dt_obs
        t_end = k * dt_obs

        # Sum new infections in [t_start, t_end)
        # Note: latent_df.t is the time at the start of each step, and
        # new_infections is how many occurred during that step.
        mask = (latent_df.t .>= t_start) .& (latent_df.t .< t_end)
        Y_true = sum(latent_df.new_infections[mask])

        # Apply observation process
        Y_obs = if sigma == 0.0 && rho == 1.0
            # Perfect observation: no noise, no underreporting
            Float64(Y_true)
        elseif sigma == 0.0
            # Binomial reporting, no additional noise
            Float64(rand(rng, Binomial(Y_true, rho)))
        else
            # Reporting with log-normal noise
            # First apply reporting rate, then log-normal noise
            Y_reported = rand(rng, Binomial(Y_true, rho))
            if Y_reported == 0
                0.0
            else
                mean_log = log(Y_reported)
                exp(mean_log + sigma * randn(rng))
            end
        end

        push!(periods, k)
        push!(t_starts, t_start)
        push!(t_ends, t_end)
        push!(true_incidence, Y_true)
        push!(reported_incidence, Y_obs)
    end

    return DataFrame(
        period = periods,
        t_start = t_starts,
        t_end = t_ends,
        true_incidence = true_incidence,
        reported_incidence = reported_incidence,
    )
end

"""Aggregate group-specific latent incidence columns into a long observation table."""
function observe_group_incidence(
    latent_df::DataFrame;
    dt_obs::Float64
)::DataFrame
    group_cols = filter(n -> startswith(String(n), "new_infections_group_"),
                        propertynames(latent_df))
    isempty(group_cols) && throw(ArgumentError("latent data contain no group-incidence columns"))
    t_max = maximum(latent_df.t)
    n_periods = floor(Int, t_max / dt_obs)
    rows = NamedTuple[]
    for period in 1:n_periods
        t_start = (period - 1) * dt_obs
        t_end = period * dt_obs
        mask = (latent_df.t .>= t_start) .& (latent_df.t .< t_end)
        for (group, col) in enumerate(group_cols)
            push!(rows, (period=period, group=group, t_start=t_start, t_end=t_end,
                         incidence=Float64(sum(latent_df[mask, col]))))
        end
    end
    return DataFrame(rows)
end

end # module

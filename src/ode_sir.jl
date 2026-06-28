# ode_sir.jl
# Mean-field ODE SIR model and grid posterior inference.
#
# CRITICAL PARAMETERIZATION NOTE:
# The ODE uses beta_ODE as an effective mass-action transmission rate:
#   dS/dt = -beta_ODE * S * I / N
#   dI/dt =  beta_ODE * S * I / N - gamma * I
#   dR/dt =  gamma * I
#
# The network SIR uses beta_true as a per-contact, per-infected-neighbor rate:
#   P(infection | n_I infected neighbors) = 1 - exp(-beta_true * dt * n_I)
#
# These are NOT on the same scale. Under uniform mixing with mean degree d_bar:
#   beta_mf_true = beta_true * d_bar
#
# This is the reference value we plot against the ODE posterior.
# If the ODE posterior concentrates near beta_mf_true, there is no bias.
# If the ODE posterior concentrates away from beta_mf_true (which happens
# under strong community structure), that IS the bias we want to show.
#
# Equivalently: the ODE sees aggregate SI products, not SI edges. Community
# structure changes the relationship between aggregate SI and SI edges, so
# fitting the ODE to community-network data recovers a pseudo-true beta_ODE*
# that differs from beta_mf_true.

module ODESIR

using DifferentialEquations
using DataFrames
using Distributions
using Statistics

export solve_ode_sir, predict_weekly_incidence,
       loglik_ode_beta, grid_posterior_beta

"""
    solve_ode_sir(params, N, I0, R0, T; saveat) -> solution

Solve the mean-field ODE SIR system.

params : SIRParams with (beta_ODE, gamma)
N      : total population size
I0     : initial number infected
R0     : initial number recovered
T      : end time
saveat : times at which to save solution

ODE system:
  dS/dt = -beta * S * I / N
  dI/dt =  beta * S * I / N - gamma * I
  dR/dt =  gamma * I

Note: beta here is beta_ODE (mass-action), not the network per-contact rate.
"""
function solve_ode_sir(
    beta_ODE::Float64,
    gamma::Float64,
    N::Int,
    I0::Float64,
    R0::Float64,
    T::Float64;
    saveat::Union{Vector{Float64}, StepRangeLen, Nothing} = nothing
)
    S0 = N - I0 - R0
    u0 = [S0, I0, R0]
    tspan = (0.0, T)

    function sir_ode!(du, u, p, t)
        S, I, R = u
        beta, gamma, Nf = p
        du[1] = -beta * S * I / Nf
        du[2] =  beta * S * I / Nf - gamma * I
        du[3] =  gamma * I
    end

    prob = ODEProblem(sir_ode!, u0, tspan, [beta_ODE, gamma, Float64(N)])

    if saveat !== nothing
        return solve(prob, Tsit5(); saveat = saveat, abstol=1e-8, reltol=1e-8)
    else
        return solve(prob, Tsit5(); abstol=1e-8, reltol=1e-8)
    end
end

"""
    predict_weekly_incidence(beta_ODE, gamma, N, I0, R0, T, dt_obs) -> Vector{Float64}

Return the ODE-predicted weekly incidence for each observation period.

Predicted incidence in period k is the drop in S between t_k and t_{k+1}:
  Y_hat_k = S(t_k) - S(t_{k+1})

This is consistent with the observation model: we aggregate new infections
over each window.
"""
function predict_weekly_incidence(
    beta_ODE::Float64,
    gamma::Float64,
    N::Int,
    I0::Float64,
    R0::Float64,
    T::Float64,
    dt_obs::Float64
)::Vector{Float64}
    n_periods = floor(Int, T / dt_obs)
    save_times = [k * dt_obs for k in 0:n_periods]
    sol = solve_ode_sir(beta_ODE, gamma, N, I0, R0, T; saveat = save_times)

    # S is component 1
    S_vals = [sol(t)[1] for t in save_times]

    # Incidence in period k = drop in S from t_{k-1} to t_k
    incidence = [max(0.0, S_vals[k] - S_vals[k+1]) for k in 1:n_periods]
    return incidence
end

"""
    loglik_ode_beta(beta_ODE, gamma, obs_df, N, I0, R0, sigma) -> Float64

Log-likelihood of observed weekly incidence given ODE SIR parameters.

Observation model (log-count Gaussian):
  log(1 + Y_k) ~ Normal(log(1 + Y_hat_k), sigma^2)

where Y_k is observed incidence in period k, Y_hat_k is ODE prediction.

This is a reasonable likelihood for count data that avoids log(0) issues
and treats residuals symmetrically on a log scale.

obs_df must have column: reported_incidence
"""
function loglik_ode_beta(
    beta_ODE::Float64,
    gamma::Float64,
    obs_df::DataFrame,
    N::Int,
    I0::Float64,
    R0::Float64,
    sigma::Float64
)::Float64
    T = maximum(obs_df.t_end)
    dt_obs = obs_df.t_end[1] - obs_df.t_start[1]
    Y_hat = predict_weekly_incidence(beta_ODE, gamma, N, I0, R0, T, dt_obs)

    Y_obs = obs_df.reported_incidence
    n_periods = min(length(Y_obs), length(Y_hat))

    ll = 0.0
    for k in 1:n_periods
        log_obs = log(1.0 + Y_obs[k])
        log_hat = log(1.0 + Y_hat[k])
        ll += -0.5 * ((log_obs - log_hat) / sigma)^2 - log(sigma)
    end
    return ll
end

"""
    grid_posterior_beta(beta_grid, gamma, obs_df_list, N, I0, R0, sigma;
                        log_prior_fn) -> DataFrame

Compute grid posterior over beta_ODE given one or more observation DataFrames.

obs_df_list: vector of observation DataFrames (multiple independent outbreaks).
  For a single outbreak, pass [obs_df].
  For R repeated outbreaks, pass [obs_df_1, ..., obs_df_R].
  Log-likelihood is summed across outbreaks (independent).

log_prior_fn: function beta -> log prior (default: uniform, returns 0.0)

Returns DataFrame with columns:
  beta, loglik, logprior, logpost, posterior_weight

posterior_weight is normalized (sums to 1.0) using log-sum-exp.
"""
function grid_posterior_beta(
    beta_grid::Vector{Float64},
    gamma::Float64,
    obs_df_list::Vector{DataFrame},
    N::Int,
    I0::Float64,
    R0::Float64,
    sigma::Float64;
    log_prior_fn::Function = (b -> 0.0)
)::DataFrame
    n = length(beta_grid)
    logliks = zeros(Float64, n)
    logpriors = zeros(Float64, n)

    for (k, beta) in enumerate(beta_grid)
        # Sum log-likelihoods across independent outbreaks
        ll = 0.0
        for obs_df in obs_df_list
            ll += loglik_ode_beta(beta, gamma, obs_df, N, I0, R0, sigma)
        end
        logliks[k] = ll
        logpriors[k] = log_prior_fn(beta)
    end

    logposts = logliks .+ logpriors

    # Log-sum-exp normalization
    log_Z = log(sum(exp.(logposts .- maximum(logposts)))) + maximum(logposts)
    posterior_weights = exp.(logposts .- log_Z)

    return DataFrame(
        beta = beta_grid,
        loglik = logliks,
        logprior = logpriors,
        logpost = logposts,
        posterior_weight = posterior_weights,
    )
end

end # module

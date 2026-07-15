#!/usr/bin/env julia
# Validate likelihood-free recovery of eta under alternative observation processes.
# A diagonal Gaussian synthetic likelihood is estimated once per eta from an ABM
# library, then evaluated on independent held-out epidemics.

using Pkg
Pkg.activate(".")
include(joinpath(@__DIR__, "..", "src", "NetworkSIRRegularization.jl"))
using .NetworkSIRRegularization
using Random, Statistics, DataFrames, CSV, Printf

const N = 1000
const K = 4
const ALPHA = -4.4
const ETA_TRUE = 2.5
const BETA = 0.04
const GAMMA = 0.10
const DT = 0.25
const T = 150.0
const N_LIBRARY = 20
const N_HELDOUT = 30
const LAMBDA_ETA = 0.1

groups = [min(K, ((i-1) ÷ (N ÷ K))+1) for i in 1:N]
sizes = [count(==(k), groups) for k in 1:K]
pop = Population(N, groups)
target_density = logistic(ALPHA)
sir = SIRParams(BETA, GAMMA)
initial_infected = [1,2,3]

function simulate_at_eta(eta, seed)
    alpha = density_balanced_alpha(target_density, eta, sizes)
    g = generate_community_graph(pop, NetworkGeneratorParams(alpha, eta, K),
                                 MersenneTwister(seed))
    simulate_network_sir(g, sir; T=T, dt=DT, initial_infected=initial_infected,
                         rng=MersenneTwister(seed+100_000), groups=groups)
end

function aggregate_vector(latent, clock)
    obs = observe_coarse_incidence(latent;
        obs_params=ObservationParams(clock, 1.0, 0.0), rng=MersenneTwister(1))
    y = Float64.(obs.true_incidence)
    mids = obs.t_start .+ clock/2
    total = sum(y)
    cum = cumsum(y) ./ max(total, 1.0)
    quantile_times = [mids[something(findfirst(cum .>= q), length(mids))] / T
                      for q in (0.10, 0.25, 0.50, 0.75, 0.90)]
    positive = findall(y .> 0)
    duration = isempty(positive) ? 0.0 : (mids[last(positive)] - mids[first(positive)]) / T
    [log1p(total), log1p(maximum(y)), mids[argmax(y)]/T, duration, quantile_times...]
end

function group_vector(latent, clock)
    obs = observe_group_incidence(latent; dt_obs=clock)
    features = Float64[]
    for k in 1:K
        d = obs[obs.group .== k, :]
        y = d.incidence
        mids = d.t_start .+ clock/2
        total = sum(y)
        cum = cumsum(y) ./ max(total, 1.0)
        positive = findall(y .> 0)
        onset = isempty(positive) ? 1.0 : mids[first(positive)]/T
        median_time = mids[something(findfirst(cum .>= 0.5), length(mids))]/T
        append!(features, [log1p(total), log1p(maximum(y)), onset,
                           mids[argmax(y)]/T, median_time])
    end
    features
end

schemes = [
    (name="aggregate_daily", feature=x -> aggregate_vector(x, 1.0)),
    (name="aggregate_weekly", feature=x -> aggregate_vector(x, 7.0)),
    (name="aggregate_monthly", feature=x -> aggregate_vector(x, 30.0)),
    (name="group_weekly", feature=x -> group_vector(x, 7.0)),
]
eta_grid = collect(0.0:0.25:3.5)

function normalized(w)
    w = copy(w); w ./= sum(w); w
end

prior_specs = NamedTuple[]
for (strength, scale) in (("weak",3.0), ("moderate",1.5), ("strong",0.75))
    push!(prior_specs, (family="half_normal", strength=strength,
                        hyperparameter=scale,
                        weight=normalized(exp.(-0.5 .* (eta_grid ./ scale).^2))))
end
for (strength, rate) in (("weak",0.33), ("moderate",1.0), ("strong",2.0))
    push!(prior_specs, (family="exponential", strength=strength,
                        hyperparameter=rate,
                        weight=normalized(exp.(-rate .* eta_grid))))
end
for (strength, spike) in (("weak",0.25), ("moderate",0.50), ("strong",0.75))
    slab = normalized(exp.(-0.5 .* (eta_grid[2:end] ./ 2.0).^2))
    weight = vcat(spike, (1-spike) .* slab)
    push!(prior_specs, (family="spike_slab", strength=strength,
                        hyperparameter=spike, weight=weight))
end

# Estimate a reusable diagonal synthetic likelihood at each eta.
library_stats = Dict{Tuple{String,Float64},NamedTuple}()
for eta in eta_grid
    library = [simulate_at_eta(eta, 20_000 + round(Int,100eta) + r)
               for r in 1:N_LIBRARY]
    for scheme in schemes
        X = reduce(hcat, (scheme.feature(x) for x in library))'
        mu = vec(mean(X, dims=1))
        variance = vec(var(X, dims=1; corrected=true))
        variance .= max.(variance, 0.02^2)
        library_stats[(scheme.name, eta)] = (mu=mu, variance=variance)
    end
end

posterior_records = NamedTuple[]
summary_records = NamedTuple[]
heldouts = [(true_eta=true_eta, heldout=h,
             latent=simulate_at_eta(true_eta, 900_000 + round(Int,100true_eta) + h))
            for true_eta in (0.0, ETA_TRUE) for h in 1:N_HELDOUT]

for item in heldouts, scheme in schemes
    true_eta, h, heldout = item.true_eta, item.heldout, item.latent
    y = scheme.feature(heldout)
    loglik = Float64[]
    for eta in eta_grid
        stat = library_stats[(scheme.name, eta)]
        ll = -0.5sum(log.(2π .* stat.variance) .+
                     (y .- stat.mu).^2 ./ stat.variance)
        push!(loglik, ll)
    end
    for prior in prior_specs
        logpost = loglik .+ log.(prior.weight)
        weights = exp.(logpost .- maximum(logpost)); weights ./= sum(weights)
        map_eta = eta_grid[argmax(weights)]
        cum = cumsum(weights)
        lo = eta_grid[findfirst(cum .>= 0.025)]
        hi = eta_grid[findfirst(cum .>= 0.975)]
        entropy = -sum(w > 0 ? w*log(w) : 0.0 for w in weights)
        eta_mean = sum(eta_grid .* weights)
        eta_sd = sqrt(sum((eta_grid .- eta_mean).^2 .* weights))
        kl_from_prior = sum(w > 0 ? w * log(w / prior.weight[i]) : 0.0
                            for (i, w) in enumerate(weights))
        expected_loss = -sum(weights .* loglik)
        mass_near_uniform = sum(weights[eta_grid .<= 0.25])
        attack_rate = heldout.R[end] / N
        push!(summary_records, (true_eta=true_eta, heldout=h, scheme=scheme.name,
            prior_family=prior.family, prior_strength=prior.strength,
            prior_hyperparameter=prior.hyperparameter, attack_rate=attack_rate,
            major_epidemic=attack_rate >= 0.5, eta_map=map_eta, eta_mean=eta_mean,
            ci_lo=lo, ci_hi=hi, ci_width=hi-lo, covered=(lo <= true_eta <= hi),
            entropy=entropy, eta_sd=eta_sd, kl_from_prior=kl_from_prior,
            expected_loss=expected_loss, mass_near_uniform=mass_near_uniform))
        for (eta, weight) in zip(eta_grid, weights)
            push!(posterior_records, (true_eta=true_eta, heldout=h, scheme=scheme.name,
                prior_family=prior.family, prior_strength=prior.strength,
                eta=eta, posterior_weight=weight))
        end
    end
end

posteriors = DataFrame(posterior_records)
heldout_summary = DataFrame(summary_records)
function summarize_validation(df)
    combine(groupby(df, [:true_eta, :scheme, :prior_family, :prior_strength]),
    [:eta_map, :true_eta] => ((x,t) -> mean(x .- t)) => :map_bias,
    [:eta_map, :true_eta] => ((x,t) -> sqrt(mean((x .- t).^2))) => :map_rmse,
    :eta_map => mean => :mean_map,
    :eta_mean => mean => :mean_posterior_mean,
    :ci_width => mean => :mean_ci_width,
    :covered => mean => :coverage,
    :entropy => mean => :mean_entropy,
    :eta_sd => mean => :mean_eta_sd,
    :kl_from_prior => mean => :mean_kl_from_prior,
    :expected_loss => mean => :mean_expected_loss,
    :mass_near_uniform => mean => :mean_mass_near_uniform,
    :heldout => (x -> length(unique(x))) => :n)
end
validation = summarize_validation(heldout_summary)

# Mean posterior across held-outs is used only as a compact visualization.
mean_posterior = combine(groupby(posteriors,
    [:true_eta, :scheme, :prior_family, :prior_strength, :eta]),
                         :posterior_weight => mean => :posterior_weight)

mkpath("output/posterior")
CSV.write("output/posterior/posterior_eta_validation_long.csv", posteriors)
CSV.write("output/posterior/posterior_eta_validation_heldouts.csv", heldout_summary)
CSV.write("output/posterior/posterior_eta_validation_summary.csv", validation)
CSV.write("output/posterior/posterior_eta_by_observation.csv", mean_posterior)

println(validation)

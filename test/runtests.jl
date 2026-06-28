#!/usr/bin/env julia
# test/runtests.jl
#
# Acceptance tests for NetworkSIRRegularization.
# Covers the six tests specified in the project plan.
#
# Run: julia --project=. test/runtests.jl

using Test
using Random
using DataFrames
using Statistics
using CSV

push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
include(joinpath(@__DIR__, "..", "src", "NetworkSIRRegularization.jl"))
using .NetworkSIRRegularization

# Shared test parameters — kept small so the suite runs fast
const N_TEST    = 200
const K_TEST    = 4
const ALPHA_TEST = -3.5   # mean degree ~5-6
const ETA_TEST  = 2.0
const BETA_TEST  = 0.05
const GAMMA_TEST = 0.10
const DT_SIM    = 0.5
const T_SIM     = 60.0
const DT_OBS    = 7.0
const SEED      = 1234

function make_groups(N, K)
    groups = Vector{Int}(undef, N)
    for i in 1:N
        groups[i] = ((i - 1) ÷ (N ÷ K)) + 1
    end
    groups[groups .> K] .= K
    return groups
end

groups = make_groups(N_TEST, K_TEST)
pop    = Population(N_TEST, groups)

# ============================================================
# Test 1: Uniform Graph Limit
# When η=0, community generator should produce approximately
# equal within- and between-group edge probabilities.
# We check that within-group edge share is close to 1/K.
# ============================================================
@testset "Test 1: Uniform graph limit (η=0)" begin
    params_eta0 = NetworkGeneratorParams(ALPHA_TEST, 0.0, K_TEST)
    rng = Random.MersenneTwister(SEED)

    # Generate several graphs and average within-group share to reduce noise
    n_graphs = 20
    shares = Float64[]
    for _ in 1:n_graphs
        g = generate_community_graph(pop, params_eta0, rng)
        push!(shares, within_group_edge_share(g, groups))
    end
    mean_share = mean(shares)
    expected   = 1.0 / K_TEST   # 0.25 for K=4

    # Allow ±5 percentage points — stochastic with N=200 so some variance expected
    @test abs(mean_share - expected) < 0.05

    # Also verify the uniform graph generator itself gives ~1/K share
    p0 = logistic(ALPHA_TEST)
    rng2 = Random.MersenneTwister(SEED + 1)
    shares_u = Float64[]
    for _ in 1:n_graphs
        g = generate_uniform_graph(N_TEST, p0, rng2)
        push!(shares_u, within_group_edge_share(g, groups))
    end
    @test abs(mean(shares_u) - expected) < 0.05
end

# ============================================================
# Test 2: Regularization Penalty
# network_regularization(η=0) == 0
# network_regularization(η>0) > 0 and strictly increasing in η
# ============================================================
@testset "Test 2: Regularization penalty" begin
    lambda = 1.0
    params0  = NetworkGeneratorParams(ALPHA_TEST, 0.0, K_TEST)
    params1  = NetworkGeneratorParams(ALPHA_TEST, 1.0, K_TEST)
    params2  = NetworkGeneratorParams(ALPHA_TEST, 2.0, K_TEST)

    pen0 = network_regularization(params0, lambda)
    pen1 = network_regularization(params1, lambda)
    pen2 = network_regularization(params2, lambda)

    @test pen0 == 0.0
    @test pen1 > pen0
    @test pen2 > pen1

    # Penalty scales as lambda * eta^2
    @test pen1 ≈ lambda * 1.0^2
    @test pen2 ≈ lambda * 2.0^2
end

# ============================================================
# Test 3: SIR Conservation
# S(t) + I(t) + R(t) = N at every time step.
# ============================================================
@testset "Test 3: SIR conservation (S+I+R=N)" begin
    rng = Random.MersenneTwister(SEED)
    g   = generate_uniform_graph(N_TEST, logistic(ALPHA_TEST), rng)

    sir_params = SIRParams(BETA_TEST, GAMMA_TEST)
    rng2 = Random.MersenneTwister(SEED + 10)
    df = simulate_network_sir(
        g, sir_params;
        T = T_SIM, dt = DT_SIM,
        initial_infected = [1, 2, 3],
        rng = rng2
    )

    totals = df.S .+ df.I .+ df.R
    @test all(totals .== N_TEST)

    # Also test on community graph
    params_c = NetworkGeneratorParams(ALPHA_TEST, ETA_TEST, K_TEST)
    rng3 = Random.MersenneTwister(SEED + 11)
    g_c  = generate_community_graph(pop, params_c, rng3)
    rng4 = Random.MersenneTwister(SEED + 12)
    df_c = simulate_network_sir(
        g_c, sir_params;
        T = T_SIM, dt = DT_SIM,
        initial_infected = [1, 2, 3],
        rng = rng4
    )

    totals_c = df_c.S .+ df_c.I .+ df_c.R
    @test all(totals_c .== N_TEST)
end

# ============================================================
# Test 4: Observation Aggregation
# Sum of reported_incidence (with perfect observation) equals
# total new infections from the latent simulation.
# ============================================================
@testset "Test 4: Observation aggregation (perfect reporting)" begin
    rng = Random.MersenneTwister(SEED)
    g   = generate_uniform_graph(N_TEST, logistic(ALPHA_TEST), rng)

    sir_params = SIRParams(BETA_TEST, GAMMA_TEST)
    rng2 = Random.MersenneTwister(SEED + 20)
    df = simulate_network_sir(
        g, sir_params;
        T = T_SIM, dt = DT_SIM,
        initial_infected = [1, 2],
        rng = rng2
    )

    obs_params = ObservationParams(DT_OBS, 1.0, 0.0)   # perfect observation
    rng3 = Random.MersenneTwister(SEED + 21)
    obs = observe_coarse_incidence(df; obs_params = obs_params, rng = rng3)

    latent_total  = sum(df.new_infections)
    observed_total = sum(obs.true_incidence)

    @test observed_total == latent_total

    # Each period's true_incidence should equal the latent count in that window
    for row in eachrow(obs)
        mask = (df.t .>= row.t_start) .& (df.t .< row.t_end)
        expected = sum(df.new_infections[mask])
        @test row.true_incidence == expected
    end
end

# ============================================================
# Test 5: Posterior Normalization
# Grid posterior weights must sum to 1.0 (within floating-point tolerance).
# ============================================================
@testset "Test 5: Posterior normalization" begin
    rng = Random.MersenneTwister(SEED)
    g   = generate_uniform_graph(N_TEST, logistic(ALPHA_TEST), rng)

    sir_params = SIRParams(BETA_TEST, GAMMA_TEST)
    rng2 = Random.MersenneTwister(SEED + 30)
    df = simulate_network_sir(
        g, sir_params;
        T = T_SIM, dt = DT_SIM,
        initial_infected = [1],
        rng = rng2
    )

    obs_params = ObservationParams(DT_OBS, 1.0, 0.0)
    rng3 = Random.MersenneTwister(SEED + 31)
    obs = observe_coarse_incidence(df; obs_params = obs_params, rng = rng3)

    beta_grid = exp.(range(log(0.01), log(3.0), length = 100))
    sigma     = 0.3
    post = grid_posterior_beta(
        beta_grid, GAMMA_TEST, [obs], N_TEST, 1.0, 0.0, sigma
    )

    @test abs(sum(post.posterior_weight) - 1.0) < 1e-10

    # Also test with multiple outbreaks
    obs2 = observe_coarse_incidence(df; obs_params = obs_params,
                                    rng = Random.MersenneTwister(SEED + 32))
    post2 = grid_posterior_beta(
        beta_grid, GAMMA_TEST, [obs, obs2], N_TEST, 1.0, 0.0, sigma
    )
    @test abs(sum(post2.posterior_weight) - 1.0) < 1e-10

    # All weights should be non-negative
    @test all(post.posterior_weight .>= 0.0)
    @test all(post2.posterior_weight .>= 0.0)
end

# ============================================================
# Test 6: Bias demonstration smoke test
# Runs the full pipeline on a small community network and checks
# that a posterior and summary CSV can be produced.
# Warns (but does not fail) if bias is not clearly visible —
# since this is a small N test, bias may be noisy.
# ============================================================
@testset "Test 6: Bias demonstration smoke test" begin
    rng = Random.MersenneTwister(SEED)
    params_c = NetworkGeneratorParams(ALPHA_TEST, ETA_TEST, K_TEST)
    g_c = generate_community_graph(pop, params_c, rng)

    stats  = compute_all_stats(g_c, groups)
    # PARAMETERIZATION: beta_mf_true = beta_true * mean_degree
    # This is the reference against which ODE posterior bias is measured.
    beta_mf_true = BETA_TEST * stats.mean_degree

    sir_params = SIRParams(BETA_TEST, GAMMA_TEST)
    rng2 = Random.MersenneTwister(SEED + 40)
    df_c = simulate_network_sir(
        g_c, sir_params;
        T = T_SIM, dt = DT_SIM,
        initial_infected = [1, 2],
        rng = rng2
    )

    obs_params = ObservationParams(DT_OBS, 1.0, 0.0)
    rng3 = Random.MersenneTwister(SEED + 41)
    obs_c = observe_coarse_incidence(df_c; obs_params = obs_params, rng = rng3)

    beta_grid = exp.(range(log(0.01), log(5.0), length = 200))
    post = grid_posterior_beta(
        beta_grid, GAMMA_TEST, [obs_c], N_TEST, 2.0, 0.0, 0.3
    )

    beta_map  = post.beta[argmax(post.posterior_weight)]
    beta_mean = sum(post.beta .* post.posterior_weight)
    bias_map  = beta_map - beta_mf_true

    # Must produce a valid DataFrame with expected columns
    @test post isa DataFrame
    @test "beta" in names(post)
    @test "posterior_weight" in names(post)
    @test "loglik" in names(post)

    # Must produce a non-degenerate posterior (not all weight on one point)
    @test sum(post.posterior_weight .> 1e-6) > 5

    # Save smoke-test summary CSV
    mkpath("output/posterior")
    summary_df = DataFrame(
        beta_mf_true = [beta_mf_true],
        beta_map     = [beta_map],
        beta_mean    = [beta_mean],
        bias_map     = [bias_map],
    )
    CSV.write("output/posterior/smoke_test_summary.csv", summary_df)
    @test isfile("output/posterior/smoke_test_summary.csv")

    # Warn if bias is not visible — small N so don't hard-fail
    bias_pct = abs(bias_map) / beta_mf_true * 100
    if bias_pct < 5.0
        @warn "Bias is small ($(round(bias_pct, digits=1))% of reference). " *
              "Consider increasing ETA_TEST or N_TEST for a clearer demonstration."
    else
        println("    Bias detected: MAP=$( round(beta_map, digits=4)), " *
                "reference=$(round(beta_mf_true, digits=4)), " *
                "bias=$(round(bias_map, sigdigits=3)) " *
                "($(round(bias_pct, digits=1))%)")
    end
end

println("\nAll tests passed.")

#!/usr/bin/env julia
# 01_generate_data.jl
#
# Generate simulated network SIR data for all experiments.
# Saves latent simulation outputs and coarse observations to output/simulations/.
#
# Run: julia --project=. scripts/01_generate_data.jl

using Pkg
Pkg.activate(".")

push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
include(joinpath(@__DIR__, "..", "src", "NetworkSIRRegularization.jl"))
using .NetworkSIRRegularization
using Random, DataFrames, CSV, Statistics, Printf, JSON3

mkpath("output/simulations")
mkpath("output/tables")

# ============================================================
# Canonical parameters (must match 00_calibrate.jl)
# ============================================================
const SEED       = 42
const N          = 1000
const K          = 4
const ALPHA      = -4.4    # gives mean degree ~12 for uniform graph
const ETA_U      = 0.0     # uniform mixing
const ETA_C      = 2.5     # community structure
const BETA_TRUE  = 0.04    # per-contact, per-infected-neighbor per day
const GAMMA_TRUE = 1.0/10  # recovery rate
const DT_SIM     = 0.25    # simulation time step (days)
const DT_OBS     = 7.0     # observation interval (days, = 1 week)
const T_SIM      = 150.0   # total simulation time (days)
const N_REPEATS  = 20      # number of independent outbreaks for variance demo
const N_ETA_REPS = 20     # reps per eta on regularization path (pooled for stable pseudo-true estimate)

# ============================================================
# Helper: build equal-sized groups
# ============================================================
function make_equal_groups(N, K)
    groups = Vector{Int}(undef, N)
    for i in 1:N
        groups[i] = ((i - 1) ÷ (N ÷ K)) + 1
    end
    groups[groups .> K] .= K
    return groups
end

groups = make_equal_groups(N, K)
pop    = Population(N, groups)
obs_params = ObservationParams(DT_OBS, 1.0, 0.0)  # perfect observation
sir_params = SIRParams(BETA_TRUE, GAMMA_TRUE)

net_params_u = NetworkGeneratorParams(ALPHA, ETA_U, K)
net_params_c = NetworkGeneratorParams(ALPHA, ETA_C, K)

println("Generating data. Parameters:")
@printf("  N=%d, K=%d, ALPHA=%.2f, ETA_C=%.1f\n", N, K, ALPHA, ETA_C)
@printf("  beta_true=%.4f, gamma_true=%.4f\n", BETA_TRUE, GAMMA_TRUE)
@printf("  DT_SIM=%.2f, DT_OBS=%.1f, T_SIM=%.0f\n", DT_SIM, DT_OBS, T_SIM)

# ============================================================
# Experiment 1: Canonical single outbreaks (uniform and community)
# Used for: fig_02 (data vs ODE fit), fig_03 (single posterior)
# ============================================================
println("\n[1] Canonical single outbreaks")

rng_u_net = Random.MersenneTwister(SEED)
g_uniform = generate_uniform_graph(N, logistic(ALPHA), rng_u_net)

rng_c_net = Random.MersenneTwister(SEED + 1)
g_community = generate_community_graph(pop, net_params_c, rng_c_net)

# Network statistics for reference
stats_u = compute_all_stats(g_uniform, groups)
stats_c = compute_all_stats(g_community, groups)

# Mean-field reference values — the scale on which ODE bias is measured.
# Under uniform mixing with mean degree d_bar:
#   beta_mf_true = beta_true * d_bar
# This is what an unbiased ODE should recover when fitting to network data
# from the corresponding network (under perfect observation).
beta_mf_uniform   = BETA_TRUE * stats_u.mean_degree
beta_mf_community = BETA_TRUE * stats_c.mean_degree

@printf("  Uniform:   mean_degree=%.2f, beta_mf_true=%.4f\n",
        stats_u.mean_degree, beta_mf_uniform)
@printf("  Community: mean_degree=%.2f, beta_mf_true=%.4f\n",
        stats_c.mean_degree, beta_mf_community)
@printf("  NOTE: Community is denser due to approximate density balancing.\n")
@printf("        beta_mf values are computed per-network for apples-to-apples comparison.\n")

# Seed infection only in group 1 (one seed per group-1 agent index)
# Seeding in one group forces epidemic to cross community boundaries — this
# is what makes the epidemic shape visibly non-ODE.
initial_infected = [1, 2, 3]   # 3 agents in group 1

rng_sim_u = Random.MersenneTwister(SEED + 10)
df_latent_u = simulate_network_sir(
    g_uniform, sir_params;
    T = T_SIM, dt = DT_SIM,
    initial_infected = initial_infected,
    rng = rng_sim_u
)

rng_sim_c = Random.MersenneTwister(SEED + 11)
df_latent_c = simulate_network_sir(
    g_community, sir_params;
    T = T_SIM, dt = DT_SIM,
    initial_infected = initial_infected,
    rng = rng_sim_c
)

rng_obs = Random.MersenneTwister(SEED + 20)
obs_u = observe_coarse_incidence(df_latent_u; obs_params = obs_params, rng = rng_obs)
obs_c = observe_coarse_incidence(df_latent_c; obs_params = obs_params, rng = rng_obs)

CSV.write("output/simulations/latent_uniform.csv",   df_latent_u)
CSV.write("output/simulations/latent_community.csv",  df_latent_c)
CSV.write("output/simulations/observed_uniform.csv",  obs_u)
CSV.write("output/simulations/observed_community.csv", obs_c)

@printf("  Uniform:   %d infections, %d weeks observed\n",
        df_latent_u.R[end], nrow(obs_u))
@printf("  Community: %d infections, %d weeks observed\n",
        df_latent_c.R[end], nrow(obs_c))

# ============================================================
# Experiment 2: Repeated community outbreaks for variance demo
# Used for: fig_04 (more data shrinks variance, not bias)
# R = N_REPEATS independent outbreaks, each on a new community graph
# ============================================================
println("\n[2] Repeated community outbreaks (N=$N_REPEATS)")

repeated_obs = DataFrame[]

for r in 1:N_REPEATS
    rng_net = Random.MersenneTwister(SEED + 1000 + r)
    g_r = generate_community_graph(pop, net_params_c, rng_net)

    rng_sim = Random.MersenneTwister(SEED + 2000 + r)
    df_r = simulate_network_sir(
        g_r, sir_params;
        T = T_SIM, dt = DT_SIM,
        initial_infected = initial_infected,
        rng = rng_sim
    )

    rng_obs_r = Random.MersenneTwister(SEED + 3000 + r)
    obs_r = observe_coarse_incidence(df_r; obs_params = obs_params, rng = rng_obs_r)
    obs_r[!, :outbreak_id] .= r
    push!(repeated_obs, obs_r)
end

df_repeated = vcat(repeated_obs...)
CSV.write("output/simulations/observed_community_repeated.csv", df_repeated)
println("  Saved $N_REPEATS outbreaks × $(nrow(repeated_obs[1])) weeks each")

# ============================================================
# Experiment 3: Regularization path — vary eta
# Used for: fig_05 (regularization path)
# ============================================================
println("\n[3] Regularization path (varying eta)")

eta_values = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0]

eta_records = []

for eta in eta_values
    net_params_eta = NetworkGeneratorParams(ALPHA, eta, K)
    total_within_shares = Float64[]
    total_mean_degrees  = Float64[]
    obs_this_eta = DataFrame[]

    for r in 1:N_ETA_REPS
        rng_net = Random.MersenneTwister(SEED + 4000 + round(Int, eta*10)*100 + r)
        g_eta = generate_community_graph(pop, net_params_eta, rng_net)
        stats_eta = compute_all_stats(g_eta, groups)

        push!(total_within_shares, stats_eta.within_group_edge_share)
        push!(total_mean_degrees,  stats_eta.mean_degree)

        rng_sim = Random.MersenneTwister(SEED + 5000 + round(Int, eta*10)*100 + r)
        df_eta = simulate_network_sir(
            g_eta, sir_params;
            T = T_SIM, dt = DT_SIM,
            initial_infected = initial_infected,
            rng = rng_sim
        )

        rng_obs_e = Random.MersenneTwister(SEED + 6000 + round(Int, eta*10)*100 + r)
        obs_eta = observe_coarse_incidence(df_eta; obs_params = obs_params, rng = rng_obs_e)
        obs_eta[!, :eta] .= eta
        obs_eta[!, :rep]  .= r
        push!(obs_this_eta, obs_eta)
    end

    push!(eta_records, (
        eta = eta,
        mean_within_share = mean(total_within_shares),
        mean_degree       = mean(total_mean_degrees),
        beta_mf_true      = BETA_TRUE * mean(total_mean_degrees),
        n_reps            = N_ETA_REPS,
    ))

    df_eta_all = vcat(obs_this_eta...)
    fname = @sprintf("output/simulations/observed_eta_%.1f.csv", eta)
    CSV.write(replace(fname, " " => ""), df_eta_all)
    @printf("  eta=%.1f: mean_degree=%.1f, within_share=%.3f, beta_mf=%.4f\n",
            eta,
            mean(total_mean_degrees),
            mean(total_within_shares),
            BETA_TRUE * mean(total_mean_degrees))
end

df_eta_summary = DataFrame(eta_records)
CSV.write("output/simulations/regularization_path_summary.csv", df_eta_summary)

# ============================================================
# Experiment 4: Grain/clock coarsening
# Used for: fig_06
# ============================================================
println("\n[4] Grain/clock coarsening (same latent, different dt_obs)")

obs_daily   = observe_coarse_incidence(df_latent_c;
    obs_params = ObservationParams(1.0, 1.0, 0.0), rng = Random.MersenneTwister(SEED + 90))
obs_weekly  = obs_c  # already computed
obs_monthly = observe_coarse_incidence(df_latent_c;
    obs_params = ObservationParams(30.0, 1.0, 0.0), rng = Random.MersenneTwister(SEED + 91))

CSV.write("output/simulations/observed_community_daily.csv",   obs_daily)
CSV.write("output/simulations/observed_community_monthly.csv", obs_monthly)
println("  Daily: $(nrow(obs_daily)) periods, weekly: $(nrow(obs_weekly)), monthly: $(nrow(obs_monthly))")

# ============================================================
# Save canonical metadata for use by downstream scripts
# ============================================================
metadata = Dict(
    "N"               => N,
    "K"               => K,
    "ALPHA"           => ALPHA,
    "ETA_U"           => ETA_U,
    "ETA_C"           => ETA_C,
    "BETA_TRUE"       => BETA_TRUE,
    "GAMMA_TRUE"      => GAMMA_TRUE,
    "DT_SIM"          => DT_SIM,
    "DT_OBS"          => DT_OBS,
    "T_SIM"           => T_SIM,
    "N_REPEATS"       => N_REPEATS,
    "initial_infected" => initial_infected,
    "beta_mf_uniform"  => beta_mf_uniform,
    "beta_mf_community" => beta_mf_community,
    "mean_degree_uniform"   => stats_u.mean_degree,
    "mean_degree_community" => stats_c.mean_degree,
)
open("output/simulations/metadata.json", "w") do f
    JSON3.write(f, metadata)
end

println("\nAll data saved to output/simulations/")
println("Metadata: output/simulations/metadata.json")
println("\nProceed to: julia --project=. scripts/02_fit_ode_grid_posterior.jl")

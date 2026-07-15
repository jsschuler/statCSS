#!/usr/bin/env julia
# 00_calibrate.jl
#
# CALIBRATION SCRIPT — run this before anything else.
#
# Purpose: verify the models are correct and find parameters that produce
# clearly visible ODE posterior bias under community network structure.
#
# What this checks:
# 1. Conservation: S+I+R = N at all times.
# 2. Network structure: eta > 0 produces visible community structure.
# 3. Mean-field reference: beta_mf_true = beta_true * mean_degree.
# 4. ODE posterior bias: posterior concentrates away from beta_mf_true
#    when fitting to community-network data.
# 5. More data shrinks variance: running R outbreaks narrows the posterior
#    without centering it on beta_mf_true.
#
# Run: julia --project=. scripts/00_calibrate.jl

using Pkg
Pkg.activate(".")

push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
include(joinpath(@__DIR__, "..", "src", "NetworkSIRRegularization.jl"))
using .NetworkSIRRegularization
using Random, DataFrames, Statistics, Printf

println("="^60)
println("CALIBRATION RUN")
println("="^60)

# ============================================================
# Parameters — tune these if bias is not visible
# ============================================================
const SEED = 42
const N = 1000
const K = 4

# Network parameters
# ALPHA controls overall edge density:
#   E[degree] ≈ logistic(ALPHA) * (N - 1)
# For mean degree ~12: logistic(ALPHA) = 12/999 => ALPHA ≈ log(12/987) ≈ -4.41
const ALPHA = -4.4  # gives mean degree ~12-13 (epidemiologically realistic)

const ETA_UNIFORM   = 0.0   # uniform mixing (maximally regularized)
const ETA_COMMUNITY = 2.5   # strong community structure

# SIR parameters (network model)
# beta_true is per-infected-neighbor per-time-unit transmission rate.
# With mean degree ~12 and gamma=1/10:
#   beta_mf_true = beta_true * d_bar ≈ 0.04 * 12 = 0.48
#   R0_network ≈ beta_mf_true / gamma = 0.48 / 0.1 = 4.8 (produces clear epidemics)
const BETA_TRUE  = 0.04   # per contact per day
const GAMMA_TRUE = 1.0/10  # recovery rate (mean infectious period = 10 days)

# Observation
const DT_SIM = 0.25      # simulation time step (days)
const DT_OBS = 7.0       # weekly observations
const T_SIM  = 150.0     # total simulation duration (days)

# ODE posterior grid
# beta_mf_true ≈ 0.04 * 12 = 0.48 for uniform network.
# Grid must comfortably contain this and any community-biased pseudo-true.
const BETA_GRID_MIN = 0.01
const BETA_GRID_MAX = 3.0
const BETA_GRID_N   = 300
const SIGMA_OBS     = 0.3   # log-count observation noise in likelihood

# ============================================================
# Helper
# ============================================================
function make_equal_groups(N, K)
    groups = Vector{Int}(undef, N)
    for i in 1:N
        groups[i] = ((i - 1) ÷ (N ÷ K)) + 1
    end
    groups[groups .> K] .= K  # handle remainder
    return groups
end

function print_divider()
    println("-"^60)
end

# ============================================================
# STEP 1: Check network generators
# ============================================================
println("\n[1] Network generator checks")
print_divider()

rng = Random.MersenneTwister(SEED)
groups = make_equal_groups(N, K)
pop = Population(N, groups)
group_sizes = [count(==(k), groups) for k in 1:K]

net_params_uniform = NetworkGeneratorParams(ALPHA, ETA_UNIFORM, K)
alpha_community = density_balanced_alpha(logistic(ALPHA), ETA_COMMUNITY, group_sizes)
net_params_community = NetworkGeneratorParams(alpha_community, ETA_COMMUNITY, K)

g_uniform   = generate_uniform_graph(N, logistic(ALPHA), rng)
g_community = generate_community_graph(pop, net_params_community, rng)

stats_u = compute_all_stats(g_uniform, groups)
stats_c = compute_all_stats(g_community, groups)

println("Uniform graph:")
@printf("  Mean degree:           %.2f\n", stats_u.mean_degree)
@printf("  Edge density:          %.4f\n", stats_u.edge_density)
@printf("  Within-group edge share: %.3f (expected ~1/K = %.3f)\n",
        stats_u.within_group_edge_share, 1/K)
@printf("  Clustering coefficient: %.4f\n", stats_u.clustering)

println("Community graph (eta=$(ETA_COMMUNITY)):")
@printf("  Mean degree:           %.2f\n", stats_c.mean_degree)
@printf("  Edge density:          %.4f\n", stats_c.edge_density)
@printf("  Within-group edge share: %.3f (should be >> 1/K = %.3f)\n",
        stats_c.within_group_edge_share, 1/K)
@printf("  Clustering coefficient: %.4f\n", stats_c.clustering)

# Network regularization penalty
reg0 = network_regularization(net_params_uniform, 1.0)
reg2 = network_regularization(net_params_community, 1.0)
@printf("\nRegularization penalty (lambda=1): eta=0 -> %.2f, eta=%.1f -> %.2f\n",
        reg0, ETA_COMMUNITY, reg2)
@assert reg0 == 0.0 "TEST FAILED: regularization at eta=0 should be 0"
@assert reg2 > reg0 "TEST FAILED: regularization should increase with eta"
println("  [PASS] regularization penalty tests")

# ============================================================
# STEP 2: Network SIR simulation + conservation check
# ============================================================
println("\n[2] Network SIR simulation")
print_divider()

sir_params = SIRParams(BETA_TRUE, GAMMA_TRUE)
initial_infected = [1, N÷K + 1]  # seed two agents in group 1

rng2 = Random.MersenneTwister(SEED + 1)
df_uniform = simulate_network_sir(
    g_uniform, sir_params;
    T = T_SIM, dt = DT_SIM,
    initial_infected = initial_infected,
    rng = rng2
)

rng3 = Random.MersenneTwister(SEED + 2)
df_community = simulate_network_sir(
    g_community, sir_params;
    T = T_SIM, dt = DT_SIM,
    initial_infected = initial_infected,
    rng = rng3
)

# Conservation check: S + I + R = N at all times
max_err_uniform   = maximum(abs.(df_uniform.S .+ df_uniform.I .+ df_uniform.R .- N))
max_err_community = maximum(abs.(df_community.S .+ df_community.I .+ df_community.R .- N))
@assert max_err_uniform   == 0 "TEST FAILED: S+I+R ≠ N in uniform simulation"
@assert max_err_community == 0 "TEST FAILED: S+I+R ≠ N in community simulation"
println("  [PASS] S+I+R = N conservation (both networks)")

total_inf_uniform   = df_uniform.R[end]
total_inf_community = df_community.R[end]
@printf("  Uniform network: total infected = %d (%.1f%%)\n",
        total_inf_uniform, 100*total_inf_uniform/N)
@printf("  Community network: total infected = %d (%.1f%%)\n",
        total_inf_community, 100*total_inf_community/N)

if total_inf_uniform < 5 || total_inf_community < 5
    println("  WARNING: Very few infections. Consider increasing beta_true,")
    println("  mean degree, or number of seeds.")
end

# ============================================================
# STEP 3: Mean-field reference value
# ============================================================
println("\n[3] Mean-field reference (parameterization alignment)")
print_divider()

d_bar_uniform   = stats_u.mean_degree
d_bar_community = stats_c.mean_degree

beta_mf_uniform   = BETA_TRUE * d_bar_uniform
beta_mf_community = BETA_TRUE * d_bar_community

@printf("  beta_true (per-contact):           %.4f\n", BETA_TRUE)
@printf("  Mean degree (uniform):             %.2f\n", d_bar_uniform)
@printf("  Mean degree (community):           %.2f\n", d_bar_community)
@printf("  beta_mf_true (uniform network):    %.4f  [= beta_true * d_bar_uniform]\n",
        beta_mf_uniform)
@printf("  beta_mf_true (community network):  %.4f  [= beta_true * d_bar_community]\n",
        beta_mf_community)
println()
println("  These are the reference values for ODE bias comparison.")
println("  If the ODE posterior under community data concentrates away from")
println("  beta_mf_community, that demonstrates the bias.")

# ============================================================
# STEP 4: Observation model
# ============================================================
println("\n[4] Observation model")
print_divider()

using ..NetworkSIRRegularization: observe_coarse_incidence
obs_params = ObservationParams(DT_OBS, 1.0, 0.0)  # perfect observation

rng4 = Random.MersenneTwister(SEED + 3)
obs_uniform   = observe_coarse_incidence(df_uniform;   obs_params = obs_params, rng = rng4)
obs_community = observe_coarse_incidence(df_community; obs_params = obs_params, rng = rng4)

# Check aggregation: sum of observed should equal total new infections
sum_obs_uniform   = sum(obs_uniform.true_incidence)
sum_obs_community = sum(obs_community.true_incidence)
total_latent_uniform   = sum(df_uniform.new_infections)
total_latent_community = sum(df_community.new_infections)

@assert sum_obs_uniform == total_latent_uniform "TEST FAILED: observation aggregation mismatch (uniform)"
@assert sum_obs_community == total_latent_community "TEST FAILED: observation aggregation mismatch (community)"
println("  [PASS] Weekly incidence sums match total latent infections")

@printf("  Uniform observation periods: %d, total observed: %d\n",
        nrow(obs_uniform), sum_obs_uniform)
@printf("  Community observation periods: %d, total observed: %d\n",
        nrow(obs_community), sum_obs_community)

# ============================================================
# STEP 5: ODE grid posterior — check normalization and bias
# ============================================================
println("\n[5] ODE grid posterior (single outbreak)")
print_divider()

beta_grid = exp.(range(log(BETA_GRID_MIN), log(BETA_GRID_MAX), length = BETA_GRID_N))
I0 = Float64(length(initial_infected))
R0 = 0.0

# Fit ODE to uniform-network data (expect posterior near beta_mf_uniform)
post_uniform = grid_posterior_beta(
    beta_grid, GAMMA_TRUE, [obs_uniform], N, I0, R0, SIGMA_OBS
)
@assert abs(sum(post_uniform.posterior_weight) - 1.0) < 1e-6 "TEST FAILED: posterior not normalized (uniform)"

# Fit ODE to community-network data (expect posterior near beta_mf_community, or biased)
post_community = grid_posterior_beta(
    beta_grid, GAMMA_TRUE, [obs_community], N, I0, R0, SIGMA_OBS
)
@assert abs(sum(post_community.posterior_weight) - 1.0) < 1e-6 "TEST FAILED: posterior not normalized (community)"
println("  [PASS] Posterior weights sum to 1.0")

# Compute posterior summaries
beta_map_uniform   = post_uniform.beta[argmax(post_uniform.posterior_weight)]
beta_map_community = post_community.beta[argmax(post_community.posterior_weight)]
beta_mean_uniform   = sum(post_uniform.beta .* post_uniform.posterior_weight)
beta_mean_community = sum(post_community.beta .* post_community.posterior_weight)

println()
println("  Fitting ODE to UNIFORM network data:")
@printf("    beta_mf_true (reference):   %.4f\n", beta_mf_uniform)
@printf("    ODE posterior MAP:          %.4f\n", beta_map_uniform)
@printf("    ODE posterior mean:         %.4f\n", beta_mean_uniform)
@printf("    Bias (MAP - reference):     %.4f\n", beta_map_uniform - beta_mf_uniform)

println()
println("  Fitting ODE to COMMUNITY network data (key result):")
@printf("    beta_mf_true (reference):   %.4f\n", beta_mf_community)
@printf("    ODE posterior MAP:          %.4f\n", beta_map_community)
@printf("    ODE posterior mean:         %.4f\n", beta_mean_community)
@printf("    Bias (MAP - reference):     %.4f  <-- should be nonzero\n",
        beta_map_community - beta_mf_community)

bias = abs(beta_map_community - beta_mf_community)
if bias < 0.01 * beta_mf_community
    println()
    println("  WARNING: Bias is very small relative to reference value.")
    println("  Suggestions to increase bias:")
    println("    - Increase ETA_COMMUNITY (stronger community structure)")
    println("    - Seed infection in one group only (initial_infected from group 1)")
    println("    - Increase DT_OBS (coarser observation grain)")
    println("    - Decrease ALPHA (sparser within-group network)")
else
    println()
    println("  [GOOD] Visible bias detected. Parameters look usable.")
end

# ============================================================
# STEP 6: More data → variance shrinks, bias remains
# ============================================================
println("\n[6] More data at same grain")
print_divider()

# Generate R = 20 independent community-network outbreaks
R_values = [1, 5, 20]
beta_maps = Float64[]
posterior_widths = Float64[]

all_obs = DataFrame[]
for r in 1:20
    rng_r = Random.MersenneTwister(SEED + 100 + r)
    g_r = generate_community_graph(pop, net_params_community, rng_r)
    rng_sim = Random.MersenneTwister(SEED + 200 + r)
    df_r = simulate_network_sir(
        g_r, sir_params;
        T = T_SIM, dt = DT_SIM,
        initial_infected = [1],
        rng = rng_sim
    )
    rng_obs = Random.MersenneTwister(SEED + 300 + r)
    obs_r = observe_coarse_incidence(df_r; obs_params = obs_params, rng = rng_obs)
    push!(all_obs, obs_r)
end

for R in R_values
    post_r = grid_posterior_beta(
        beta_grid, GAMMA_TRUE, all_obs[1:R], N, I0, R0, SIGMA_OBS
    )
    map_r = post_r.beta[argmax(post_r.posterior_weight)]

    # Posterior width: 95% credible interval
    cum_weight = cumsum(post_r.posterior_weight)
    lo_idx = findfirst(cum_weight .>= 0.025)
    hi_idx = findfirst(cum_weight .>= 0.975)
    width_r = post_r.beta[hi_idx] - post_r.beta[lo_idx]

    push!(beta_maps, map_r)
    push!(posterior_widths, width_r)

    @printf("  R=%2d outbreaks: MAP = %.4f, 95%% CI width = %.4f, bias = %.4f\n",
            R, map_r, width_r, map_r - beta_mf_community)
end

# Check variance is shrinking
if posterior_widths[end] < posterior_widths[1]
    println("  [GOOD] Posterior variance shrinks with more data.")
else
    println("  WARNING: Posterior variance did not shrink as expected.")
end

# Check bias is not going away
bias_R1  = abs(beta_maps[1] - beta_mf_community)
bias_R20 = abs(beta_maps[end] - beta_mf_community)
if bias_R20 > 0.2 * bias_R1  # still has substantial bias
    println("  [GOOD] Bias remains (does not vanish with more data).")
else
    println("  WARNING: Bias appears to have vanished. Check parameterization.")
end

# ============================================================
# Summary
# ============================================================
println("\n" * "="^60)
println("CALIBRATION SUMMARY")
println("="^60)
@printf("beta_true (per-contact):            %.4f\n", BETA_TRUE)
@printf("gamma_true:                         %.4f\n", GAMMA_TRUE)
@printf("Mean degree (uniform):              %.2f\n", d_bar_uniform)
@printf("Mean degree (community):            %.2f\n", d_bar_community)
@printf("beta_mf_true (reference):           %.4f\n", beta_mf_community)
@printf("ODE posterior MAP (1 outbreak):     %.4f\n", beta_maps[1])
@printf("Bias (MAP - reference):             %.4f\n", beta_maps[1] - beta_mf_community)
@printf("95%% CI width (1 outbreak):          %.4f\n", posterior_widths[1])
@printf("95%% CI width (20 outbreaks):        %.4f\n", posterior_widths[end])
println()
println("If bias is visible and variance shrinks: proceed to 01_generate_data.jl")
println("If not: tune ALPHA, ETA_COMMUNITY, or DT_OBS above and rerun.")

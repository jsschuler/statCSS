#!/usr/bin/env julia
# 02_fit_ode_grid_posterior.jl
#
# Fit mean-field ODE SIR to all coarse observation datasets via grid posterior.
# Saves posterior summaries and per-beta posterior weights to output/posterior/.
#
# PARAMETERIZATION REMINDER:
# ODE beta_ODE is a mass-action rate. Network beta_true is per-contact.
# The unbiased reference for comparison is:
#   beta_mf_true = beta_true * mean_degree
# We compare the ODE posterior MAP and mean to this reference.
# Visible bias means: posterior concentrates near some beta_ODE* ≠ beta_mf_true.
#
# Run: julia --project=. scripts/02_fit_ode_grid_posterior.jl

using Pkg
Pkg.activate(".")

push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
include(joinpath(@__DIR__, "..", "src", "NetworkSIRRegularization.jl"))
using .NetworkSIRRegularization
using DataFrames, CSV, JSON3, Printf, Statistics

mkpath("output/posterior")

# ============================================================
# Load metadata
# ============================================================
meta = JSON3.read("output/simulations/metadata.json")

N            = Int(meta["N"])
K            = Int(meta["K"])
BETA_TRUE    = Float64(meta["BETA_TRUE"])
GAMMA_TRUE   = Float64(meta["GAMMA_TRUE"])
DT_OBS       = Float64(meta["DT_OBS"])
T_SIM        = Float64(meta["T_SIM"])
N_REPEATS    = Int(meta["N_REPEATS"])
beta_mf_u    = Float64(meta["beta_mf_uniform"])
beta_mf_c    = Float64(meta["beta_mf_community"])

# Initial conditions: 3 infected, 0 recovered
I0 = 3.0
R0 = 0.0

# ODE posterior grid over beta_ODE
# Range must contain both beta_mf_uniform and beta_mf_community.
# Using log spacing for better resolution near smaller values.
const BETA_GRID_MIN = 0.01
const BETA_GRID_MAX = 4.0
const BETA_GRID_N   = 2000
const SIGMA_OBS     = 0.3

beta_grid = exp.(range(log(BETA_GRID_MIN), log(BETA_GRID_MAX), length = BETA_GRID_N))

println("ODE Grid Posterior Inference")
@printf("  beta_mf_uniform:   %.4f (unbiased reference for uniform data)\n", beta_mf_u)
@printf("  beta_mf_community: %.4f (unbiased reference for community data)\n", beta_mf_c)
@printf("  beta grid: [%.3f, %.3f], %d points, log-spaced\n",
        BETA_GRID_MIN, BETA_GRID_MAX, BETA_GRID_N)

# ============================================================
# Posterior summary helper
# ============================================================
function posterior_summary(post_df, beta_mf_ref, label)
    w = post_df.posterior_weight
    b = post_df.beta

    beta_map  = b[argmax(w)]
    beta_mean = sum(b .* w)
    beta_var  = sum((b .- beta_mean).^2 .* w)
    beta_sd   = sqrt(beta_var)

    cum = cumsum(w)
    lo = b[findfirst(cum .>= 0.025)]
    hi = b[findfirst(cum .>= 0.975)]

    bias_map  = beta_map  - beta_mf_ref
    bias_mean = beta_mean - beta_mf_ref

    @printf("  %s\n", label)
    @printf("    beta_mf_true (reference): %.4f\n", beta_mf_ref)
    @printf("    MAP:    %.4f  (bias = %+.4f)\n", beta_map, bias_map)
    @printf("    Mean:   %.4f  (bias = %+.4f)\n", beta_mean, bias_mean)
    @printf("    SD:     %.4f\n", beta_sd)
    @printf("    95%% CI: [%.4f, %.4f]  (width = %.4f)\n", lo, hi, hi - lo)

    return (
        label        = label,
        beta_mf_true = beta_mf_ref,
        beta_map     = beta_map,
        beta_mean    = beta_mean,
        beta_sd      = beta_sd,
        ci_lo        = lo,
        ci_hi        = hi,
        ci_width     = hi - lo,
        bias_map     = bias_map,
        bias_mean    = bias_mean,
    )
end

all_summaries = []

# ============================================================
# Fit 1: Uniform network (single outbreak)
# Expected: low bias — ODE should recover approximately beta_mf_uniform
# ============================================================
println("\n[1] Uniform network (single outbreak)")
obs_u = CSV.read("output/simulations/observed_uniform.csv", DataFrame)
post_u = grid_posterior_beta(beta_grid, GAMMA_TRUE, [obs_u], N, I0, R0, SIGMA_OBS)
CSV.write("output/posterior/posterior_uniform.csv", post_u)
push!(all_summaries, posterior_summary(post_u, beta_mf_u, "uniform_single"))

# ============================================================
# Fit 2: Community network (single outbreak)
# Key result: ODE posterior should be biased relative to beta_mf_community
# ============================================================
println("\n[2] Community network (single outbreak) — KEY BIAS RESULT")
obs_c = CSV.read("output/simulations/observed_community.csv", DataFrame)
post_c = grid_posterior_beta(beta_grid, GAMMA_TRUE, [obs_c], N, I0, R0, SIGMA_OBS)
CSV.write("output/posterior/posterior_community_1.csv", post_c)
push!(all_summaries, posterior_summary(post_c, beta_mf_c, "community_1outbreak"))

# ============================================================
# Fit 3: Community network (increasing outbreaks)
# KEY DEMO: more data shrinks posterior variance but not bias
# ============================================================
println("\n[3] Community network — more data, same grain")

df_repeated = CSV.read("output/simulations/observed_community_repeated.csv", DataFrame)
outbreak_ids = sort(unique(df_repeated.outbreak_id))

R_values = [1, 5, 10, 20]

for R in R_values
    ids = outbreak_ids[1:R]
    obs_list = [filter(row -> row.outbreak_id == id, df_repeated) for id in ids]
    post_r = grid_posterior_beta(beta_grid, GAMMA_TRUE, obs_list, N, I0, R0, SIGMA_OBS)
    fname = "output/posterior/posterior_community_$(R)outbreaks.csv"
    CSV.write(fname, post_r)
    push!(all_summaries, posterior_summary(post_r, beta_mf_c, "community_$(R)outbreaks"))
end

# ============================================================
# Fit 4: Regularization path — one posterior per eta value
# ============================================================
println("\n[4] Regularization path posteriors")

eta_values = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0]
eta_summary_df = CSV.read("output/simulations/regularization_path_summary.csv", DataFrame)

for eta in eta_values
    fname = replace(@sprintf("output/simulations/observed_eta_%.1f.csv", eta), " " => "")
    df_eta = CSV.read(fname, DataFrame)

    rep_ids = sort(unique(df_eta.rep))
    obs_list = [filter(row -> row.rep == r, df_eta) for r in rep_ids]

    # Pool all reps to get a stable estimate of the pseudo-true beta* for this eta.
    # A single outbreak is too noisy to reliably locate the pseudo-true parameter.
    post_eta = grid_posterior_beta(beta_grid, GAMMA_TRUE, obs_list, N, I0, R0, SIGMA_OBS)

    out_fname = replace(@sprintf("output/posterior/posterior_eta_%.1f.csv", eta), " " => "")
    CSV.write(out_fname, post_eta)

    # Look up beta_mf_true for this eta
    row = filter(r -> r.eta == eta, eta_summary_df)
    ref = isempty(row) ? NaN : row.beta_mf_true[1]
    push!(all_summaries, posterior_summary(post_eta, ref, "eta_$(eta)"))
end

# ============================================================
# Save combined summary table
# ============================================================
df_summary = DataFrame(all_summaries)
CSV.write("output/posterior/posterior_summary.csv", df_summary)

println("\n" * "="^60)
println("POSTERIOR SUMMARY TABLE")
println("="^60)
println(df_summary[:, [:label, :beta_mf_true, :beta_map, :bias_map, :ci_width]])

println("\nAll posteriors saved to output/posterior/")
println("Summary: output/posterior/posterior_summary.csv")
println("\nProceed to: julia --project=. scripts/03_make_figures.jl")

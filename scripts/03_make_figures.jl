#!/usr/bin/env julia
# 03_make_figures.jl
#
# Generate all static figures using CairoMakie.
# Reads from output/simulations/ and output/posterior/.
# Writes to figures/.
#
# NOTE on posterior plots: beta_grid is log-spaced, so bin widths vary.
# We plot posterior_weight directly (y-axis: "posterior weight per grid point"),
# not a continuous density. The peak location and relative widths are what matter.
#
# Run: julia --project=. scripts/03_make_figures.jl

using Pkg
Pkg.activate(".")

push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
include(joinpath(@__DIR__, "..", "src", "NetworkSIRRegularization.jl"))
using .NetworkSIRRegularization
using CairoMakie
using DataFrames, CSV, JSON3, Printf, Statistics, Random

mkpath("figures")

# ============================================================
# Load metadata
# ============================================================
meta         = JSON3.read("output/simulations/metadata.json")
N            = Int(meta["N"])
BETA_TRUE    = Float64(meta["BETA_TRUE"])
GAMMA_TRUE   = Float64(meta["GAMMA_TRUE"])
DT_OBS       = Float64(meta["DT_OBS"])
T_SIM        = Float64(meta["T_SIM"])
beta_mf_u    = Float64(meta["beta_mf_uniform"])
beta_mf_c    = Float64(meta["beta_mf_community"])
I0           = 3.0
R0           = 0.0

# Shared style
const FONT_SIZE  = 14
const TITLE_SIZE = 15
const LABEL_SIZE = 13
const W          = 900
const H          = 540

const COL_TRUE     = :black
const COL_BIASED   = RGBf(0.80, 0.20, 0.10)
const COL_UNBIASED = RGBf(0.15, 0.45, 0.75)
const COL_NETWORK  = RGBf(0.20, 0.60, 0.30)
const COL_OBS      = RGBf(0.50, 0.30, 0.70)

println("Making figures...")

# ============================================================
# Figure 1: Model workflow diagram
# ============================================================
println("  fig_01: model workflow")

fig1 = Figure(size = (W, 300), fontsize = FONT_SIZE)
ax1 = Axis(fig1[1, 1],
    limits = (0, 10, 0, 1),
    xticksvisible = false, yticksvisible = false,
    xticklabelsvisible = false, yticklabelsvisible = false,
    leftspinevisible = false, rightspinevisible = false,
    topspinevisible = false, bottomspinevisible = false,
    title = "Generative Statistical Model for Network Epidemics",
    titlesize = TITLE_SIZE,
)

boxes = [
    (0.3,  0.55, "Mechanism\n(theory)"),
    (1.7,  0.55, "Network\nABM"),
    (3.1,  0.55, "Latent\nepidemic"),
    (4.5,  0.55, "Observation\nprocess"),
    (5.9,  0.55, "Coarse\ndata"),
    (7.3,  0.55, "ODE\ncalibration"),
    (8.7,  0.55, "Posterior\n+ criticism"),
]
box_w, box_h = 0.95, 0.52

for (cx, cy, label) in boxes
    poly!(ax1,
        Point2f[(cx - box_w/2, cy - box_h/2), (cx + box_w/2, cy - box_h/2),
                (cx + box_w/2, cy + box_h/2), (cx - box_w/2, cy + box_h/2)],
        color = RGBAf(0.92, 0.94, 0.98, 1.0),
        strokecolor = :gray40, strokewidth = 1.5)
    text!(ax1, cx, cy; text = label, align = (:center, :center),
          fontsize = 10.5, color = :black)
end

# Arrows between boxes (use text "→" to avoid arrows API version issues)
for b in boxes[1:end-1]
    x_mid = b[1] + box_w/2 + 0.07
    text!(ax1, x_mid, 0.55; text = "→", align = (:center, :center),
          fontsize = 14, color = :gray40)
end

text!(ax1, 4.5, 0.06;
    text = "← coarsening hides network structure →",
    align = (:center, :bottom), fontsize = 9.5, color = :gray35, font = :italic)

save("figures/fig_01_model_workflow.png", fig1, px_per_unit = 2)

# ============================================================
# Figure 2: Network ABM data vs ODE fit
# ============================================================
println("  fig_02: network data vs ODE fit")

df_latent = CSV.read("output/simulations/latent_community.csv",    DataFrame)
obs_c     = CSV.read("output/simulations/observed_community.csv",  DataFrame)
post_c1   = CSV.read("output/posterior/posterior_community_1.csv", DataFrame)

beta_map_c = post_c1.beta[argmax(post_c1.posterior_weight)]
obs_mids   = obs_c.t_start .+ DT_OBS / 2

ode_fit = predict_weekly_incidence(beta_map_c, GAMMA_TRUE, N, I0, R0, T_SIM, DT_OBS)
ode_ref = predict_weekly_incidence(beta_mf_c,  GAMMA_TRUE, N, I0, R0, T_SIM, DT_OBS)
n_fit   = min(length(ode_fit), length(obs_mids))

fig2 = Figure(size = (W, H), fontsize = FONT_SIZE)
ax2  = Axis(fig2[1, 1],
    xlabel = "Time (days)",
    ylabel = "New infections per week",
    title  = "Network SIR Data vs. ODE Fit (Community Network)",
    titlesize = TITLE_SIZE)

barplot!(ax2, obs_mids, Float64.(obs_c.true_incidence);
    width = DT_OBS * 0.85, color = (COL_NETWORK, 0.35),
    label = "True weekly incidence (latent)")

scatter!(ax2, obs_mids, Float64.(obs_c.reported_incidence);
    color = COL_OBS, markersize = 9, marker = :diamond,
    label = "Observed counts")

lines!(ax2, obs_mids[1:n_fit], ode_fit[1:n_fit];
    color = COL_BIASED, linewidth = 2.5,
    label = @sprintf("ODE fit  (β̂_ODE = %.3f)", beta_map_c))

lines!(ax2, obs_mids[1:n_fit], ode_ref[1:n_fit];
    color = COL_UNBIASED, linewidth = 2.0, linestyle = :dash,
    label = @sprintf("ODE at β_mf,true = %.3f", beta_mf_c))

axislegend(ax2, position = :rt, labelsize = LABEL_SIZE)

save("figures/fig_02_network_vs_ode_incidence.png", fig2, px_per_unit = 2)

# ============================================================
# Figure 3: Posterior bias — single outbreak
# ============================================================
println("  fig_03: posterior bias (single outbreak)")

post_u1 = CSV.read("output/posterior/posterior_uniform.csv", DataFrame)

fig3 = Figure(size = (W, H), fontsize = FONT_SIZE)
# NOTE: x-axis expanded to [0.0, 2.2] to include β_true = 0.04 (the microscopic
# per-contact rate). This makes the three-level parameterization visible:
#   β_true (microscopic) << β_mf,true (uniform-mixing reference) > β_ODE* (biased)
# The wide gap between β_true and the ODE quantities is intentional — it shows
# that β_true and β_ODE are not on the same scale and should not be compared directly.
ax3  = Axis(fig3[1, 1],
    xlabel = "Transmission rate",
    ylabel = "Posterior weight",
    title  = "ODE Posterior Bias under Community Structure",
    titlesize = TITLE_SIZE,
    limits = (0.0, 2.2, nothing, nothing))

lines!(ax3, post_u1.beta, post_u1.posterior_weight;
    color = COL_UNBIASED, linewidth = 2.5,
    label = "Posterior: uniform network data")

lines!(ax3, post_c1.beta, post_c1.posterior_weight;
    color = COL_BIASED, linewidth = 2.5,
    label = "Posterior: community network data")

# β_true: the microscopic per-contact rate. NOT directly comparable to β_ODE —
# shown here to make the unit difference explicit.
vlines!(ax3, [BETA_TRUE];
    color = :gray40, linestyle = :dot, linewidth = 2.0,
    label = @sprintf("β_true (per-contact) = %.3f", BETA_TRUE))

vlines!(ax3, [beta_mf_u];
    color = COL_UNBIASED, linestyle = :dash, linewidth = 1.8,
    label = @sprintf("β·d̄ ref. (uniform) = %.3f", beta_mf_u))

vlines!(ax3, [beta_mf_c];
    color = COL_TRUE, linestyle = :dash, linewidth = 1.8,
    label = @sprintf("β·d̄ ref. (community) = %.3f", beta_mf_c))

axislegend(ax3, position = :rt, labelsize = LABEL_SIZE - 1)

# Annotate β_true to explain why it is far left
text!(ax3, BETA_TRUE + 0.02, maximum(post_c1.posterior_weight) * 0.85;
    text = "β_true = $(BETA_TRUE)\n(per-contact rate;\nnot ODE-scale)",
    fontsize = 10, color = :gray40, font = :italic, align = (:left, :center))

# Bias annotation
beta_map_annot = post_c1.beta[argmax(post_c1.posterior_weight)]
text!(ax3, beta_map_annot - 0.02, maximum(post_c1.posterior_weight) * 0.55;
    text = "ODE posterior\n(biased)",
    fontsize = 11, color = COL_BIASED, font = :italic, align = (:right, :center))
text!(ax3, beta_mf_c + 0.03, maximum(post_c1.posterior_weight) * 0.20;
    text = "Counterfactual\nreference (β·d̄)",
    fontsize = 11, color = :black, font = :italic, align = (:left, :center))

save("figures/fig_03_posterior_bias_single.png", fig3, px_per_unit = 2)

# ============================================================
# Figure 4: More data shrinks variance, not bias (KEY FIGURE)
# ============================================================
println("  fig_04: more data shrinks variance not bias")

R_values = [1, 5, 10, 20]
r_labels = ["R = 1 outbreak", "R = 5 outbreaks", "R = 10 outbreaks", "R = 20 outbreaks"]

# NOTE: vertically stacked panels, one per R value, shared x-axis.
# Each panel uses the same x limits so narrowing is directly visible across rows.
# This is cleaner than overlaid curves: the reader can see each posterior
# individually rather than trying to disentangle four curves on one axis.
fig4 = Figure(size = (W, 680), fontsize = FONT_SIZE)
Label(fig4[0, :],
    text = "More Data at the Same Grain: Variance Shrinks, Bias Remains",
    fontsize = TITLE_SIZE, font = :bold)

# Compute shared y-max from the single-outbreak posterior (widest, lowest peak)
# so all panels use the same y scale — making the peak growth visible.
post_R1 = CSV.read("output/posterior/posterior_community_1outbreaks.csv", DataFrame)
post_R20 = CSV.read("output/posterior/posterior_community_20outbreaks.csv", DataFrame)
y_max = maximum(post_R20.posterior_weight) * 1.12

axes4 = Axis[]
for (i, (R, lbl)) in enumerate(zip(R_values, r_labels))
    post_r = CSV.read("output/posterior/posterior_community_$(R)outbreaks.csv", DataFrame)

    is_bottom = (i == length(R_values))
    ax = Axis(fig4[i, 1],
        ylabel = "Post. weight",
        xlabel = is_bottom ? "Transmission rate" : "",
        title  = lbl,
        titlesize = LABEL_SIZE,
        limits = (0.0, 2.2, 0.0, y_max),
        xticklabelsvisible = is_bottom,
        xticksvisible      = is_bottom,
    )
    push!(axes4, ax)

    # Posterior curve
    lines!(ax, post_r.beta, post_r.posterior_weight;
        color = COL_BIASED, linewidth = 2.2)

    # Shared reference lines on every panel
    vlines!(ax, [BETA_TRUE];
        color = :gray40, linestyle = :dot, linewidth = 1.8)
    vlines!(ax, [beta_mf_c];
        color = COL_TRUE, linestyle = :dash, linewidth = 1.8)

    # Label reference lines only on top panel to avoid repetition
    if i == 1
        text!(ax, BETA_TRUE + 0.02, y_max * 0.85;
            text = "β_true = $(BETA_TRUE)\n(per-contact; not ODE-scale)",
            fontsize = 9, color = :gray40, font = :italic, align = (:left, :top))
        text!(ax, beta_mf_c + 0.03, y_max * 0.85;
            text = "β·d̄ ref. = $(round(beta_mf_c, digits=3))\n(counterfactual)",
            fontsize = 9, color = :black, font = :italic, align = (:left, :top))
    end

    # 95% CI width annotation on right side of each panel
    cum  = cumsum(post_r.posterior_weight)
    lo   = post_r.beta[findfirst(cum .>= 0.025)]
    hi   = post_r.beta[findfirst(cum .>= 0.975)]
    text!(ax, 1.95, y_max * 0.75;
        text = @sprintf("95%% CI\nwidth = %.3f", hi - lo),
        fontsize = 9, color = :gray30, align = (:right, :center))
end

# Tighten row gaps
rowgap!(fig4.layout, 4)

save("figures/fig_04_more_data_shrinks_variance_not_bias.png", fig4, px_per_unit = 2)

# ============================================================
# Figure 5: Regularization path
# ============================================================
println("  fig_05: regularization path")

df_eta_sum  = CSV.read("output/simulations/regularization_path_summary.csv", DataFrame)
df_post_sum = CSV.read("output/posterior/posterior_summary.csv", DataFrame)
eta_rows    = filter(r -> startswith(r.label, "eta_"), df_post_sum)
eta_rows[!, :eta] = parse.(Float64, replace.(eta_rows.label, "eta_" => ""))
sort!(eta_rows, :eta)

fig5 = Figure(size = (W, 640), fontsize = FONT_SIZE)
Label(fig5[0, :],
    text = "Regularization Path: Community Structure and ODE Bias",
    fontsize = TITLE_SIZE, font = :bold)

ax5a = Axis(fig5[1, 1],
    xlabel = "",
    ylabel = "Within-group edge share",
    ylabelcolor      = COL_NETWORK,
    yticklabelcolor  = COL_NETWORK,
    xticklabelsvisible = false,
    xticksvisible      = false)

ax5b = Axis(fig5[2, 1],
    xlabel = "Community structure (η)",
    ylabel = "Posterior bias  (MAP − β_mf,true)",
    ylabelcolor      = COL_BIASED,
    yticklabelcolor  = COL_BIASED)

lines!(ax5a, df_eta_sum.eta, df_eta_sum.mean_within_share;
    color = COL_NETWORK, linewidth = 2.2)
scatter!(ax5a, df_eta_sum.eta, df_eta_sum.mean_within_share;
    color = COL_NETWORK, markersize = 10)
hlines!(ax5a, [1/4];
    color = (COL_NETWORK, 0.4), linestyle = :dash, linewidth = 1.5)
text!(ax5a, 0.08, 1/4 + 0.015;
    text = "uniform mixing (1/K = 0.25)", fontsize = 10, color = (COL_NETWORK, 0.75))

lines!(ax5b, eta_rows.eta, eta_rows.bias_map;
    color = COL_BIASED, linewidth = 2.2)
scatter!(ax5b, eta_rows.eta, eta_rows.bias_map;
    color = COL_BIASED, markersize = 10, marker = :rect)
hlines!(ax5b, [0.0];
    color = (COL_BIASED, 0.35), linestyle = :dash, linewidth = 1.5)

rowgap!(fig5.layout, 6)
linkxaxes!(ax5a, ax5b)

save("figures/fig_05_regularization_path_network_stats.png", fig5, px_per_unit = 2)

# ============================================================
# Figure 6: Grain/clock coarsening
# ============================================================
println("  fig_06: grain/clock coarsening")

obs_daily   = CSV.read("output/simulations/observed_community_daily.csv",   DataFrame)
obs_weekly  = CSV.read("output/simulations/observed_community.csv",         DataFrame)
obs_monthly = CSV.read("output/simulations/observed_community_monthly.csv", DataFrame)

fig6  = Figure(size = (W + 100, H), fontsize = FONT_SIZE)
ttls  = ["Daily (Δt = 1 day)", "Weekly (Δt = 7 days)", "Monthly (Δt = 30 days)"]
dsets = [obs_daily, obs_weekly, obs_monthly]
cols  = [COL_UNBIASED, COL_BIASED, COL_TRUE]

ymax = maximum(obs_daily.true_incidence) * 1.1

for (i, (obs, ttl, col)) in enumerate(zip(dsets, ttls, cols))
    ax = Axis(fig6[1, i],
        xlabel = "Time (days)",
        ylabel = i == 1 ? "Incidence (new infections per period)" : "",
        title  = ttl, titlesize = 12,
        ylabelsize = 11)
    mids = obs.t_start .+ (obs.t_end .- obs.t_start) ./ 2
    w    = obs.t_end[1] - obs.t_start[1]
    barplot!(ax, mids, Float64.(obs.true_incidence);
        width = w * 0.85,
        color = (col, 0.65), strokecolor = col, strokewidth = 0.5)
    ylims!(ax, 0, ymax)
end

Label(fig6[0, :],
    text = "Same Latent Epidemic — Different Observation Grains",
    fontsize = TITLE_SIZE, font = :bold)
Label(fig6[2, :],
    text = "Coarser grains erase event ordering and mask network structure",
    fontsize = 11, color = :gray30, font = :italic)

save("figures/fig_06_grain_clock_coarsening.png", fig6, px_per_unit = 2)

# ============================================================
# Figure 7: Perturbation smoothing
# ============================================================
println("  fig_07: perturbation smoothing")

df_rep = CSV.read("output/simulations/observed_community_repeated.csv", DataFrame)
peak_incidences = Float64[]
for id in sort(unique(df_rep.outbreak_id))
    sub = filter(r -> r.outbreak_id == id, df_rep)
    push!(peak_incidences, Float64(maximum(sub.true_incidence)))
end

rng_perturb = Random.MersenneTwister(99)
n_bins = 8

fig7 = Figure(size = (W, H), fontsize = FONT_SIZE)
ax7  = Axis(fig7[1, 1],
    xlabel = "Weekly peak incidence (new infections)",
    ylabel = "Density (normalized)",
    title  = "Perturbation Stability: Histogram Averaging",
    titlesize = TITLE_SIZE)

# Perturbed histograms in the background
for _ in 1:40
    noise     = randn(rng_perturb, length(peak_incidences)) .* 4.0
    perturbed = peak_incidences .+ noise
    hist!(ax7, perturbed;
        bins = n_bins, normalization = :pdf,
        color = (COL_BIASED, 0.05), strokecolor = :transparent)
end

# Raw histogram
hist!(ax7, peak_incidences;
    bins = n_bins, normalization = :pdf,
    color = (:transparent, 0.0),
    strokecolor = COL_TRUE, strokewidth = 2.5,
    label = "Raw histogram (20 outbreaks)")

# Legend proxies
lines!(ax7, Float64[], Float64[];
    color = (COL_BIASED, 0.5), linewidth = 8,
    label = "Perturbed histograms (40 draws)")
lines!(ax7, Float64[], Float64[];
    color = COL_TRUE, linewidth = 2.5,
    label = "Raw histogram (observed data)")

axislegend(ax7, position = :rt, labelsize = LABEL_SIZE)

text!(ax7, minimum(peak_incidences) + 5,
    0.010;
    text = "Features that survive small\nperturbations are stable summaries",
    fontsize = 11, color = :gray30, font = :italic)

save("figures/fig_07_histogram_perturbation_smoothing.png", fig7, px_per_unit = 2)

println("\nAll figures saved to figures/")
println("  fig_01_model_workflow.png")
println("  fig_02_network_vs_ode_incidence.png")
println("  fig_03_posterior_bias_single.png")
println("  fig_04_more_data_shrinks_variance_not_bias.png")
println("  fig_05_regularization_path_network_stats.png")
println("  fig_06_grain_clock_coarsening.png")
println("  fig_07_histogram_perturbation_smoothing.png")

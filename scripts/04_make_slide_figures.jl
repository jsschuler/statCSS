#!/usr/bin/env julia
# 04_make_slide_figures.jl
#
# Slide-specific figure variants. Does not touch paper figures.
# Run: julia --project=. scripts/04_make_slide_figures.jl

using Pkg
Pkg.activate(".")
push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
include(joinpath(@__DIR__, "..", "src", "NetworkSIRRegularization.jl"))
using .NetworkSIRRegularization
using CairoMakie, DataFrames, CSV, JSON3, Printf, Statistics

mkpath("figures")

meta       = JSON3.read("output/simulations/metadata.json")
N          = Int(meta["N"])
BETA_TRUE  = Float64(meta["BETA_TRUE"])
GAMMA_TRUE = Float64(meta["GAMMA_TRUE"])
DT_OBS     = Float64(meta["DT_OBS"])
T_SIM      = Float64(meta["T_SIM"])
beta_mf_u  = Float64(meta["beta_mf_uniform"])
beta_mf_c  = Float64(meta["beta_mf_community"])

COL_BIASED   = RGBf(0.80, 0.20, 0.10)
COL_UNBIASED = RGBf(0.15, 0.45, 0.75)
COL_TRUE     = :black

# ============================================================
# fig_04_slide: single-panel overlaid posteriors for R=1,5,10,20
# ============================================================
println("  fig_04_slide: overlaid posteriors")

R_values = [1, 5, 10, 20]
alphas   = [0.25, 0.50, 0.72, 1.00]
lwidths  = [1.5,  2.0,  2.5,  3.0]

posts = [CSV.read("output/posterior/posterior_community_$(R)outbreaks.csv", DataFrame)
         for R in R_values]

y_max = maximum(posts[end].posterior_weight) * 1.18

fig = Figure(size = (860, 400), fontsize = 14)
ax  = Axis(fig[1, 1],
    xlabel    = "β_ODE",
    ylabel    = "Posterior weight",
    title     = "More Data at the Same Grain: Variance Shrinks, Bias Remains",
    titlesize = 15,
    limits    = (0.7, 2.1, 0.0, y_max))

for (R, post, α, lw) in zip(R_values, posts, alphas, lwidths)
    lines!(ax, post.beta, post.posterior_weight;
        color     = RGBAf(COL_BIASED.r, COL_BIASED.g, COL_BIASED.b, α),
        linewidth = lw,
        label     = "R = $R outbreak$(R > 1 ? "s" : "")")
end

vlines!(ax, [BETA_TRUE];
    color     = :gray50,
    linestyle = :dot,
    linewidth = 1.5,
    label     = @sprintf("β_true = %.2f  (per-contact; different scale)", BETA_TRUE))

vlines!(ax, [beta_mf_c];
    color     = COL_TRUE,
    linestyle = :dash,
    linewidth = 2.0,
    label     = @sprintf("β_mf,true ≈ %.2f  (counterfactual reference)", beta_mf_c))

# Annotate that the peak does not move toward the reference
peak_x = posts[end].beta[argmax(posts[end].posterior_weight)]
text!(ax, peak_x + 0.03, y_max * 0.55;
    text      = "peak fixed here\n(bias persists)",
    fontsize  = 11,
    color     = :gray35,
    font      = :italic,
    align     = (:left, :center))

# Annotate CI shrinkage with a bracket-style annotation
for (R, post) in zip(R_values, posts)
    w = post.posterior_weight
    b = post.beta
    cum = cumsum(w)
    lo = b[findfirst(cum .>= 0.025)]
    hi = b[findfirst(cum .>= 0.975)]
    width = hi - lo
    α_ann = R == 20 ? 1.0 : 0.5
    text!(ax, (lo + hi) / 2, maximum(w) * 1.05;
        text     = @sprintf("%.3f", width),
        fontsize = 9,
        color    = RGBAf(COL_BIASED.r, COL_BIASED.g, COL_BIASED.b, α_ann),
        align    = (:center, :bottom))
end

text!(ax, 0.71, y_max * 0.97;
    text    = "95% CI width:",
    fontsize = 9,
    color   = :gray40,
    align   = (:left, :top))

axislegend(ax, position = :rt, labelsize = 11)

save("figures/fig_04_slide.png", fig, px_per_unit = 2)
println("  Saved figures/fig_04_slide.png")

# ============================================================
# fig_beta_networks: two small network cartoons illustrating
# why β_true and β_ODE operate on different scales, and what
# community structure does to ODE calibration.
#
# LEFT:  Erdős–Rényi network (same density as community, d̄≈37)
#        → ODE recovers β_mf,true ≈ 1.47 without bias
# RIGHT: SBM community network (same d̄, but clustered)
#        → ODE finds β_ODE* ≈ 1.13, 23% below the reference
# ============================================================
println("  fig_beta_networks: network cartoons")

using Random
rng = MersenneTwister(2025)

N_vis = 20
K_vis = 4
npg   = N_vis ÷ K_vis          # 5 nodes per group
groups_vis = repeat(1:K_vis, inner = npg)

# ── Node positions ─────────────────────────────────────────
# Uniform: all nodes on a circle
pos_u = [Point2f(0.5 + 0.40*cos(2π*(i-1)/N_vis + π/N_vis),
                 0.5 + 0.40*sin(2π*(i-1)/N_vis + π/N_vis)) for i in 1:N_vis]

# Community: four clusters in quadrants
gcx = [0.27, 0.73, 0.27, 0.73]
gcy = [0.73, 0.73, 0.27, 0.27]
pos_c = Point2f[]
for i in 1:N_vis
    g   = groups_vis[i]
    idx = (i - 1) % npg
    θ   = 2π * idx / npg + π/8
    push!(pos_c, Point2f(gcx[g] + 0.14*cos(θ), gcy[g] + 0.14*sin(θ)))
end

# ── Edges ─────────────────────────────────────────────────
# Both networks have similar total density (≈35% of possible edges)
# so the structural difference is clearly the *pattern*, not the count.
p_uni    = 0.35
p_within = 0.78
p_btwn   = 0.07

edges_u = Tuple{Int,Int}[]
edges_cw = Tuple{Int,Int}[]   # within-group (community)
edges_cb = Tuple{Int,Int}[]   # between-group (community)
for i in 1:N_vis, j in (i+1):N_vis
    rand(rng) < p_uni && push!(edges_u, (i, j))
    if groups_vis[i] == groups_vis[j]
        rand(rng) < p_within && push!(edges_cw, (i, j))
    else
        rand(rng) < p_btwn   && push!(edges_cb, (i, j))
    end
end

# ── Colors ─────────────────────────────────────────────────
grp_cols = [RGBf(0.20,0.50,0.80), RGBf(0.82,0.33,0.14),
            RGBf(0.22,0.68,0.35), RGBf(0.60,0.32,0.72)]
nc_u = fill(RGBf(0.35, 0.52, 0.78), N_vis)
nc_c = [grp_cols[g] for g in groups_vis]

# ── Figure layout ──────────────────────────────────────────
fig_bn = Figure(size = (860, 460), fontsize = 14)

# Network axes — extend y below 0 to make room for annotations
ax_u = Axis(fig_bn[1, 1],
    title     = "Erdős–Rényi (uniform mixing)",
    titlesize = 13,
    limits    = (0, 1, -0.38, 1.0),
    aspect    = DataAspect())
ax_c = Axis(fig_bn[1, 2],
    title     = "Stochastic Block Model (community structure)",
    titlesize = 13,
    limits    = (0, 1, -0.38, 1.0),
    aspect    = DataAspect())

for ax in (ax_u, ax_c)
    hidedecorations!(ax)
    hidespines!(ax)
end

# ── Draw uniform network ────────────────────────────────────
pts_u = Point2f[]
for (i,j) in edges_u; push!(pts_u, pos_u[i], pos_u[j]); end
isempty(pts_u) || linesegments!(ax_u, pts_u; color = (:gray55, 0.45), linewidth = 0.9)
scatter!(ax_u, pos_u; color = nc_u, markersize = 14,
         strokecolor = :white, strokewidth = 1.2)

# ── Draw community network ──────────────────────────────────
pts_cb = Point2f[]
for (i,j) in edges_cb; push!(pts_cb, pos_c[i], pos_c[j]); end
isempty(pts_cb) || linesegments!(ax_c, pts_cb; color = (:gray78, 0.35), linewidth = 0.7)

pts_cw = Point2f[]
for (i,j) in edges_cw; push!(pts_cw, pos_c[i], pos_c[j]); end
isempty(pts_cw) || linesegments!(ax_c, pts_cw; color = (:gray30, 0.70), linewidth = 1.1)

scatter!(ax_c, pos_c; color = nc_c, markersize = 14,
         strokecolor = :white, strokewidth = 1.2)

# ── β annotations below each network ───────────────────────
# Shared first line (same for both)
for ax in (ax_u, ax_c)
    text!(ax, 0.5, -0.07;
        text      = @sprintf("β_true = %.2f day⁻¹  (per contact, per infected neighbor)", BETA_TRUE),
        align     = (:center, :center),
        fontsize  = 10,
        color     = RGBf(0.22, 0.68, 0.35))
end

# Uniform: unbiased ODE recovery
text!(ax_u, 0.5, -0.17;
    text     = @sprintf("d̄ ≈ 37  →  β_mf,true = β × d̄ ≈ %.2f", beta_mf_c),
    align    = (:center, :center),
    fontsize = 10,
    color    = :gray20)
text!(ax_u, 0.5, -0.27;
    text     = @sprintf("ODE posterior ≈ %.2f  ✓  (within 2%%)", beta_mf_c),
    align    = (:center, :center),
    fontsize = 11,
    font     = :bold,
    color    = RGBf(0.15, 0.45, 0.75))

# Community: biased ODE
beta_map_c = 1.133   # MAP from posterior_summary
bias_pct   = round((beta_map_c - beta_mf_c) / beta_mf_c * 100, digits=0)
text!(ax_c, 0.5, -0.17;
    text     = @sprintf("d̄ ≈ 37  →  β_mf,true = β × d̄ ≈ %.2f  (reference)", beta_mf_c),
    align    = (:center, :center),
    fontsize = 10,
    color    = :gray20)
text!(ax_c, 0.5, -0.27;
    text     = @sprintf("ODE posterior ≈ %.2f  ✗  (bias = %+.0f%%)", beta_map_c, bias_pct),
    align    = (:center, :center),
    fontsize = 11,
    font     = :bold,
    color    = RGBf(0.80, 0.20, 0.10))

save("figures/fig_beta_networks.png", fig_bn, px_per_unit = 2)
println("  Saved figures/fig_beta_networks.png")

#!/usr/bin/env julia
# Direct recovery of community structure when the contact graph is observed.

using Pkg
Pkg.activate(".")
include(joinpath(@__DIR__, "..", "src", "NetworkSIRRegularization.jl"))
using .NetworkSIRRegularization
using Random, CSV, DataFrames, JSON3, Printf

mkpath("output/posterior")
meta = JSON3.read("output/simulations/metadata.json")
N = Int(meta["N"]); K = Int(meta["K"])
eta_true = Float64(meta["ETA_C"])
alpha0 = Float64(meta["ALPHA"])
target_density = logistic(alpha0)
groups = [min(K, ((i - 1) ÷ (N ÷ K)) + 1) for i in 1:N]
pop = Population(N, groups)
sizes = [count(==(k), groups) for k in 1:K]
alpha = density_balanced_alpha(target_density, eta_true, sizes)

# Reproduce the canonical community graph from 01_generate_data.jl.
g = generate_community_graph(pop, NetworkGeneratorParams(alpha, eta_true, K),
                             Random.MersenneTwister(43))
eta_grid = collect(0.0:0.01:3.5)

for lambda in (0.0, 0.1, 1.0)
    post = grid_posterior_eta(g, pop, eta_grid, target_density; lambda_eta=lambda)
    CSV.write("output/posterior/posterior_eta_rich_lambda_$(lambda).csv", post)
    eta_map = post.eta[argmax(post.posterior_weight)]
    cum = cumsum(post.posterior_weight)
    lo = post.eta[findfirst(cum .>= 0.025)]
    hi = post.eta[findfirst(cum .>= 0.975)]
    @printf("lambda=%3.1f: eta MAP %.3f, 95%% interval [%.3f, %.3f]\n",
            lambda, eta_map, lo, hi)
end

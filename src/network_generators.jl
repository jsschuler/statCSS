# network_generators.jl
# Uniform and community network generators for the SIR regularization demo.

module NetworkGenerators

using Graphs
using Random
using ..SIRTypes: Population, NetworkGeneratorParams

export logistic, generate_uniform_graph, generate_community_graph,
       network_regularization, density_balanced_alpha, expected_edge_density

"""
    logistic(x)

Standard logistic function: 1 / (1 + exp(-x)).
"""
logistic(x::Float64) = 1.0 / (1.0 + exp(-x))

"""Exact expected SBM edge density for finite group sizes."""
function expected_edge_density(
    alpha::Float64,
    eta::Float64,
    group_sizes::Vector{Int}
)::Float64
    K = length(group_sizes)
    N = sum(group_sizes)
    total_dyads = N * (N - 1) / 2
    within_dyads = sum(n * (n - 1) / 2 for n in group_sizes)
    between_dyads = total_dyads - within_dyads
    p_within = logistic(alpha + eta)
    p_between = logistic(alpha - eta / (K - 1))
    return (within_dyads * p_within + between_dyads * p_between) / total_dyads
end

"""Solve for the SBM intercept that preserves a requested expected density."""
function density_balanced_alpha(
    target_density::Float64,
    eta::Float64,
    group_sizes::Vector{Int};
    tol::Float64 = 1e-12
)::Float64
    0.0 < target_density < 1.0 || throw(ArgumentError("target_density must lie in (0, 1)"))
    lo, hi = -30.0, 30.0
    while hi - lo > tol
        mid = (lo + hi) / 2
        if expected_edge_density(mid, eta, group_sizes) < target_density
            lo = mid
        else
            hi = mid
        end
    end
    return (lo + hi) / 2
end

"""
    generate_uniform_graph(N, p, rng) -> SimpleGraph

Erdős–Rényi G(N, p) graph. This is the maximally regularized (eta=0) case.
Every dyad has the same edge probability p, encoding uniform mixing.
"""
function generate_uniform_graph(N::Int, p::Float64, rng::AbstractRNG)::SimpleGraph
    g = SimpleGraph(N)
    for i in 1:N
        for j in (i+1):N
            if rand(rng) < p
                add_edge!(g, i, j)
            end
        end
    end
    return g
end

"""
    generate_community_graph(pop, params, rng) -> SimpleGraph

Stochastic block model with community regularization parameter eta.

logit P(A_ij = 1) = alpha + eta * B_ij

B_ij = 1           if same group
B_ij = -1/(K-1)    if different group

The normalization b = 1/(K-1) is chosen so that for equal group sizes the
average B_ij across all dyads is zero:
  (within-group dyads) * 1 + (between-group dyads) * (-1/(K-1)) ≈ 0

This means eta entering only changes within/between edge ratio without
systematically raising or lowering overall density. Check with network_stats.

When eta = 0: all dyads have probability logistic(alpha), uniform mixing.
When eta > 0: within-group edges more probable, between-group less.
"""
function generate_community_graph(
    pop::Population,
    params::NetworkGeneratorParams,
    rng::AbstractRNG
)::SimpleGraph
    N = pop.N
    groups = pop.groups
    K = params.n_groups
    alpha = params.alpha
    eta = params.eta
    b = 1.0 / (K - 1)  # normalization to stabilize density across eta values

    g = SimpleGraph(N)
    for i in 1:N
        for j in (i+1):N
            B_ij = (groups[i] == groups[j]) ? 1.0 : -b
            p_ij = logistic(alpha + eta * B_ij)
            if rand(rng) < p_ij
                add_edge!(g, i, j)
            end
        end
    end
    return g
end

"""
    network_regularization(params, lambda_eta) -> Float64

Penalty term for network community structure.
Returns lambda_eta * eta^2.

When eta = 0 (uniform graph): penalty is 0.
Larger eta (stronger community structure) incurs larger penalty.
This encodes the prior belief that uniform mixing is the simpler model,
and departures from it must be justified by the data.
"""
function network_regularization(params::NetworkGeneratorParams, lambda_eta::Float64)::Float64
    return lambda_eta * params.eta^2
end

end # module

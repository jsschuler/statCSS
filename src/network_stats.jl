# network_stats.jl
# Intensive (per-agent or per-dyad) network statistics for comparing graphs
# across regularization levels.
#
# All statistics are normalized so they are comparable as N changes.
# Do not use raw edge counts or raw triangle counts in comparison figures.

module NetworkStats

using Graphs
using Statistics

export mean_degree, degree_variance, edge_density, clustering_coefficient,
       within_group_edge_share, compute_all_stats

"""
    mean_degree(g) -> Float64

Average degree per node. Intensive: comparable across N.
"""
mean_degree(g::SimpleGraph)::Float64 = mean(degree(g))

"""
    degree_variance(g) -> Float64

Variance of degree distribution.
"""
degree_variance(g::SimpleGraph)::Float64 = var(degree(g))

"""
    edge_density(g) -> Float64

Fraction of possible edges that exist. In [0, 1].
Intensive: comparable across N.
"""
function edge_density(g::SimpleGraph)::Float64
    N = nv(g)
    max_edges = N * (N - 1) / 2
    return max_edges > 0 ? ne(g) / max_edges : 0.0
end

"""
    clustering_coefficient(g) -> Float64

Global clustering coefficient: fraction of paths of length 2 that close
into a triangle. Returns 0.0 if no paths of length 2 exist.
Intensive: comparable across N.
"""
function clustering_coefficient(g::SimpleGraph)::Float64
    N = nv(g)
    triangles = 0
    paths2 = 0
    for v in vertices(g)
        nbrs = neighbors(g, v)
        k = length(nbrs)
        if k < 2
            continue
        end
        paths2 += k * (k - 1) ÷ 2
        for i in 1:length(nbrs)
            for j in (i+1):length(nbrs)
                if has_edge(g, nbrs[i], nbrs[j])
                    triangles += 1
                end
            end
        end
    end
    return paths2 > 0 ? triangles / paths2 : 0.0
end

"""
    within_group_edge_share(g, groups) -> Float64

Fraction of edges that connect two nodes in the same group.
Intensive: comparable across N (it is a share, not a count).

Values near 1/K (where K = number of groups) indicate uniform mixing.
Values near 1.0 indicate strong community structure.
"""
function within_group_edge_share(g::SimpleGraph, groups::Vector{Int})::Float64
    total_edges = ne(g)
    if total_edges == 0
        return 0.0
    end
    within = 0
    for e in edges(g)
        if groups[src(e)] == groups[dst(e)]
            within += 1
        end
    end
    return within / total_edges
end

"""
    compute_all_stats(g, groups) -> NamedTuple

Convenience function returning all network statistics as a named tuple.
"""
function compute_all_stats(g::SimpleGraph, groups::Vector{Int})
    return (
        mean_degree = mean_degree(g),
        degree_variance = degree_variance(g),
        edge_density = edge_density(g),
        clustering = clustering_coefficient(g),
        within_group_edge_share = within_group_edge_share(g, groups),
        n_nodes = nv(g),
        n_edges = ne(g),
    )
end

end # module

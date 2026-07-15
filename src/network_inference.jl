module NetworkInference

using DataFrames
using Graphs
using ..SIRTypes: Population
using ..NetworkGenerators: density_balanced_alpha, logistic

export grid_posterior_eta

"""
    grid_posterior_eta(g, pop, eta_grid, target_density; lambda_eta=0)

Exact grid posterior for the SBM community parameter when the contact graph and
group memberships are observed. The intercept is recalibrated at every eta so
all candidates have the same expected density. `lambda_eta * eta^2` is used as
a log-prior penalty, making the regularization operational rather than merely
descriptive.
"""
function grid_posterior_eta(
    g::SimpleGraph,
    pop::Population,
    eta_grid::Vector{Float64},
    target_density::Float64;
    lambda_eta::Float64 = 0.0
)::DataFrame
    K = maximum(pop.groups)
    sizes = [count(==(k), pop.groups) for k in 1:K]
    loglik = zeros(length(eta_grid))
    alphas = zeros(length(eta_grid))
    possible_within = sum(n * (n - 1) ÷ 2 for n in sizes)
    possible_total = pop.N * (pop.N - 1) ÷ 2
    possible_between = possible_total - possible_within
    observed_within = count(e -> pop.groups[src(e)] == pop.groups[dst(e)], edges(g))
    observed_between = ne(g) - observed_within

    for (q, eta) in enumerate(eta_grid)
        alpha = density_balanced_alpha(target_density, eta, sizes)
        alphas[q] = alpha
        p_within = logistic(alpha + eta)
        p_between = logistic(alpha - eta / (K - 1))
        loglik[q] = observed_within * log(p_within) +
                    (possible_within - observed_within) * log1p(-p_within) +
                    observed_between * log(p_between) +
                    (possible_between - observed_between) * log1p(-p_between)
    end

    logprior = -lambda_eta .* eta_grid.^2
    logpost = loglik .+ logprior
    shifted = logpost .- maximum(logpost)
    weights = exp.(shifted)
    weights ./= sum(weights)
    return DataFrame(eta=eta_grid, alpha=alphas, loglik=loglik,
                     logprior=logprior, logpost=logpost,
                     posterior_weight=weights)
end

end

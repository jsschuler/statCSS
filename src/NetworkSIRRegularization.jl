# NetworkSIRRegularization.jl
# Main module. Includes all submodules in dependency order.

module NetworkSIRRegularization

include("sir_types.jl")
include("network_generators.jl")
include("network_stats.jl")
include("sir_simulation.jl")
include("observation.jl")
include("ode_sir.jl")

using .SIRTypes
using .NetworkGenerators
using .NetworkStats
using .SIRSimulation
using .Observation
using .ODESIR

export SIRParams, Population, NetworkGeneratorParams, ObservationParams
export logistic, generate_uniform_graph, generate_community_graph, network_regularization
export mean_degree, degree_variance, edge_density, clustering_coefficient,
       within_group_edge_share, compute_all_stats
export simulate_network_sir, simulate_network_sir_frames
export observe_coarse_incidence
export solve_ode_sir, predict_weekly_incidence, loglik_ode_beta, grid_posterior_beta

end # module

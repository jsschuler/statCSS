# sir_types.jl
# Core types for network SIR regularization project.

module SIRTypes

export SIRParams, Population, NetworkGeneratorParams, ObservationParams

"""
    SIRParams(beta, gamma)

Transmission parameters. IMPORTANT: interpretation depends on context.

- In network SIR: beta is the per-contact, per-infected-neighbor transmission
  rate (units: probability per contact per time step).
- In ODE SIR: beta_ODE is an effective mass-action rate (units: 1/(N * time)).

These are NOT automatically on the same scale. The mean-field effective rate
implied by network beta under uniform mixing is:

    beta_mf_true = beta_true * mean_degree

where mean_degree is the average number of contacts per agent. This is the
reference value against which ODE posterior bias is measured.
"""
struct SIRParams
    beta::Float64
    gamma::Float64
end

"""
    Population(N, groups)

N agents divided into groups (group membership by index).
groups[i] ∈ 1:K gives the group of agent i.
"""
struct Population
    N::Int
    groups::Vector{Int}
end

"""
    NetworkGeneratorParams(alpha, eta, n_groups)

Parameters for community network generator.

alpha: baseline log-odds of an edge (controls overall density).
eta: community structure strength.
  - eta = 0: uniform random graph (Erdős–Rényi with p = logistic(alpha)).
  - eta > 0: within-group edges are more likely than between-group.

The logit probability of edge (i,j) is:
  alpha + eta * B_ij

where B_ij = 1 if same group, B_ij = -1/(K-1) if different group.
The choice b = 1/(K-1) keeps expected edge density approximately stable as eta
varies, because the average B_ij across all dyads is zero when groups are equal
size.
"""
struct NetworkGeneratorParams
    alpha::Float64
    eta::Float64
    n_groups::Int
end

"""
    ObservationParams(dt_obs, reporting_rate, noise_sd)

Controls the coarsening of latent simulation data into observed data.
dt_obs: observation interval in simulation time units.
reporting_rate: fraction of true cases observed (1.0 = perfect reporting).
noise_sd: standard deviation of log-normal observation noise (0.0 = no noise).
"""
struct ObservationParams
    dt_obs::Float64
    reporting_rate::Float64
    noise_sd::Float64
end

end # module

# sir_simulation.jl
# Discrete-time network SIR simulation.
#
# PARAMETERIZATION NOTE (critical):
# beta here is the per-contact, per-infected-neighbor transmission rate.
# P(susceptible i gets infected | t infected neighbors) = 1 - exp(-beta * dt * n_I)
#
# This is NOT the same as ODE beta_ODE, which is a mass-action aggregate rate.
# The mean-field effective rate under uniform mixing with mean degree d_bar is:
#   beta_mf_true = beta_true * d_bar
# That is the reference for comparing ODE posterior bias.

module SIRSimulation

using Graphs
using DataFrames
using Random
using ..SIRTypes: SIRParams, Population

export simulate_network_sir, simulate_network_sir_frames

"""
    simulate_network_sir(g, params; T, dt, initial_infected, rng) -> DataFrame

Discrete-time network SIR simulation with synchronous update.

Returns a DataFrame with columns:
  t, S, I, R, new_infections, new_recoveries

Arguments:
  g                 : contact network (SimpleGraph)
  params            : SIRParams with (beta, gamma)
  T                 : total simulation time
  dt                : time step size
  initial_infected  : vector of node indices initially infected
  rng               : random number generator (pass explicit seed for reproducibility)

Beta interpretation: probability of transmission per infected neighbor per unit
time is 1 - exp(-beta * dt). So for small beta*dt, it is approximately beta*dt.
"""
function simulate_network_sir(
    g::SimpleGraph,
    params::SIRParams;
    T::Float64,
    dt::Float64,
    initial_infected::Vector{Int},
    rng::AbstractRNG
)::DataFrame
    N = nv(g)
    # Status codes: 1=S, 2=I, 3=R
    status = ones(Int, N)  # all susceptible
    for i in initial_infected
        status[i] = 2
    end

    beta = params.beta
    gamma = params.gamma
    # Precompute per-step recovery probability
    p_recover = 1.0 - exp(-gamma * dt)

    times = Float64[]
    S_counts = Int[]
    I_counts = Int[]
    R_counts = Int[]
    new_inf = Int[]
    new_rec = Int[]

    t = 0.0
    n_steps = round(Int, T / dt)

    for step in 0:n_steps
        S = count(==(1), status)
        I = count(==(2), status)
        R = count(==(3), status)
        push!(times, t)
        push!(S_counts, S)
        push!(I_counts, I)
        push!(R_counts, R)

        if step == n_steps
            push!(new_inf, 0)
            push!(new_rec, 0)
            break
        end

        # Compute transitions (synchronous: use current status for all decisions)
        new_status = copy(status)
        infections_this_step = 0
        recoveries_this_step = 0

        for i in 1:N
            if status[i] == 1  # susceptible
                # Count infected neighbors using current status
                n_infected_neighbors = 0
                for j in neighbors(g, i)
                    if status[j] == 2
                        n_infected_neighbors += 1
                    end
                end
                if n_infected_neighbors > 0
                    p_infect = 1.0 - exp(-beta * dt * n_infected_neighbors)
                    if rand(rng) < p_infect
                        new_status[i] = 2
                        infections_this_step += 1
                    end
                end
            elseif status[i] == 2  # infected
                if rand(rng) < p_recover
                    new_status[i] = 3
                    recoveries_this_step += 1
                end
            end
        end

        push!(new_inf, infections_this_step)
        push!(new_rec, recoveries_this_step)
        status = new_status
        t += dt
    end

    return DataFrame(
        t = times,
        S = S_counts,
        I = I_counts,
        R = R_counts,
        new_infections = new_inf,
        new_recoveries = new_rec,
    )
end

"""
    simulate_network_sir_frames(g, pop, params; T, dt, initial_infected, rng, positions)
    -> DataFrame

Like simulate_network_sir but also returns per-agent state at each time step,
for animation.

Returns DataFrame with columns:
  t, agent_id, x, y, group, status

status is :S, :I, or :R.
positions: optional (N x 2) matrix of (x, y) coordinates for each agent.
  If not provided, uses group-clustered random layout.
"""
function simulate_network_sir_frames(
    g::SimpleGraph,
    pop::Population,
    params::SIRParams;
    T::Float64,
    dt::Float64,
    initial_infected::Vector{Int},
    rng::AbstractRNG,
    positions::Union{Matrix{Float64}, Nothing} = nothing
)::DataFrame
    N = nv(g)
    groups = pop.groups
    K = maximum(groups)

    # Generate group-clustered positions if not provided
    if positions === nothing
        positions = Matrix{Float64}(undef, N, 2)
        # Place group centers around a circle, agents scattered within
        for i in 1:N
            g_idx = groups[i]
            angle = 2π * (g_idx - 1) / K
            cx = 0.5 + 0.3 * cos(angle)
            cy = 0.5 + 0.3 * sin(angle)
            positions[i, 1] = cx + 0.1 * randn(rng)
            positions[i, 2] = cy + 0.1 * randn(rng)
        end
    end

    status = ones(Int, N)
    for i in initial_infected
        status[i] = 2
    end

    beta = params.beta
    gamma = params.gamma
    p_recover = 1.0 - exp(-gamma * dt)

    status_sym = [:S, :I, :R]

    rows = []
    t = 0.0
    n_steps = round(Int, T / dt)

    for step in 0:n_steps
        for i in 1:N
            push!(rows, (
                t = t,
                agent_id = i,
                x = positions[i, 1],
                y = positions[i, 2],
                group = groups[i],
                status = status_sym[status[i]],
            ))
        end

        if step == n_steps
            break
        end

        new_status = copy(status)
        for i in 1:N
            if status[i] == 1
                n_inf_nbrs = count(j -> status[j] == 2, neighbors(g, i))
                if n_inf_nbrs > 0
                    p_infect = 1.0 - exp(-beta * dt * n_inf_nbrs)
                    if rand(rng) < p_infect
                        new_status[i] = 2
                    end
                end
            elseif status[i] == 2
                if rand(rng) < p_recover
                    new_status[i] = 3
                end
            end
        end
        status = new_status
        t += dt
    end

    return DataFrame(rows)
end

end # module

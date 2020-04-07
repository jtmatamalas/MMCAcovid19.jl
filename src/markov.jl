
"""
    update_prob!(Pᵢᵍ::Array{Float64, 2},
                 Pᵢᴬᵍ::Array{Float64, 2},
                 Pᵢᴵᵍ::Array{Float64, 2},
                 Sᵢᵍ::Array{Float64, 2},
                 τᵢᵍ::Array{Float64, 2},
                 epi_params::Epidemic_Params,
                 population::Population_Params,
                 M::Int64,
                 G::Int64,
                 edgelist::Array{Int64, 2},
                 Rᵢⱼ::Array{Float64, 1},
                 C::Array{Float64, 2},
                 κ₀::Float64,
                 ϕ::Float64,
                 pᵍ::Array{Float64, 1},
                 t::Int64,
                 tᶜ::Int64)

Updates the probabilities of the model using the equations described in the paper.
"""
function update_prob!(Pᵢᵍ::Array{Float64, 2},
                      Pᵢᴬᵍ::Array{Float64, 2},
                      Pᵢᴵᵍ::Array{Float64, 2},
                      Sᵢᵍ::Array{Float64, 2},
                      τᵢᵍ::Array{Float64, 2},
                      epi_params::Epidemic_Params,
                      population::Population_Params,
                      M::Int64,
                      G::Int64,
                      edgelist::Array{Int64, 2},
                      Rᵢⱼ::Array{Float64, 1},
                      C::Array{Float64, 2},
                      κ₀::Float64,
                      ϕ::Float64,
                      pᵍ::Array{Float64, 1},
                      t::Int64,
                      tᶜ::Int64)

    # Shortcuts to parameters
    ρˢʸ = epi_params.ρˢᵍ
    ρᴱʸ = epi_params.ρᴱᵍ
    ρᴬʸ = epi_params.ρᴬᵍ
    ρᴵʸ = epi_params.ρᴵᵍ
    ρᴿʸ = epi_params.ρᴿᵍ
    ηᵍ = epi_params.ηᵍ
    αᵍ = epi_params.αᵍ
    μᵍ = epi_params.μᵍ
    γᵍ = epi_params.γᵍ
    χᵍ = epi_params.χᵍ
    ψᵍ = epi_params.ψᵍ
    ωᵍ = epi_params.ωᵍ
    ρˢᵍ = epi_params.ρˢᵍ
    ρᴱᵍ = epi_params.ρᴱᵍ
    ρᴬᵍ = epi_params.ρᴬᵍ
    ρᴵᵍ = epi_params.ρᴵᵍ
    ρᴴᵍ = epi_params.ρᴴᵍ
    ρᴰᵍ = epi_params.ρᴰᵍ
    ρᴿᵍ = epi_params.ρᴿᵍ
    CHᵢᵍ = epi_params.CHᵢᵍ

    # Intervention at time tᶜ
    if tᶜ == t
        pᵍ[:] .= (1 - κ₀) .* pᵍ
        population.kᵍ[:] .= [population.σ - 1, (1 - κ₀) * population.kᵍ[2] + κ₀ * (population.σ - 1), population.σ - 1]
        update_population_param!(population, pᵍ, edgelist, Rᵢⱼ, M, G)
    end

    # Get P and compute Q
    compute_P!(Pᵢᵍ, Pᵢᴬᵍ, Pᵢᴵᵍ, Sᵢᵍ, pᵍ, ρˢᵍ, ρᴬᵍ, ρᴵᵍ,
               epi_params.Qᵢᵍ, population.nᵢᵍ_eff, population.mobilityᵍ,
               population.normᵍ, epi_params.βᴬ, epi_params.βᴵ, epi_params.νᵍ,
               edgelist, Rᵢⱼ, C, M, G, t)

    # Comptue τᵢᵍ
    for indx_e in 1:length(Rᵢⱼ)
        i = edgelist[indx_e, 1]
        j = edgelist[indx_e, 2]

        @simd for g in 1:G
            τᵢᵍ[g, i] += Rᵢⱼ[indx_e] * Pᵢᵍ[g, j]
        end
    end

    # Update probabilities
    for i in 1:M

        # Compute secure households
        CHᵢ = 0.0
        if tᶜ == t
            @simd for g in 1:G
                CHᵢ += (ρˢᵍ[g, i, t] + ρᴿᵍ[g, i, t]) * population.nᵢᵍ[g, i]
            end
            CHᵢ = κ₀ * (1 - ϕ) * (CHᵢ / population.nᵢ[i]) ^ population.σ
            CHᵢᵍ[:, i] .= ρˢᵍ[:, i, t] .* CHᵢ
        end

        # Update compartmental probabilities
        @simd for g in 1:G
            Πᵢᵍ = (1 - pᵍ[g]) * Pᵢᵍ[g, i] + pᵍ[g] * τᵢᵍ[g, i]

            # Spreading
            ρˢᵍ[g, i, t + 1] = (1 - Πᵢᵍ) * (1 - CHᵢ) * ρˢᵍ[g, i, t]
            ρᴱᵍ[g, i, t + 1] = (1 - ηᵍ[g]) * ρᴱᵍ[g, i, t] + Πᵢᵍ * (1 - CHᵢ) * ρˢᵍ[g, i, t]
            ρᴬᵍ[g, i, t + 1] = (1 - αᵍ[g]) * ρᴬᵍ[g, i, t] + ηᵍ[g] * ρᴱᵍ[g, i, t]
            ρᴵᵍ[g, i, t + 1] = (1 - μᵍ[g]) * ρᴵᵍ[g, i, t] + αᵍ[g] * ρᴬᵍ[g, i, t]
            ρᴴᵍ[g, i, t + 1] = μᵍ[g] * γᵍ[g] * ρᴵᵍ[g, i, t] + ωᵍ[g] * (1 - ψᵍ[g]) * ρᴴᵍ[g, i, t] + (1 - ωᵍ[g]) * (1 - χᵍ[g]) * ρᴴᵍ[g, i, t]
            ρᴰᵍ[g, i, t + 1] = ωᵍ[g] * ψᵍ[g] * ρᴴᵍ[g, i, t] + ρᴰᵍ[g, i, t]
            ρᴿᵍ[g, i, t + 1] = μᵍ[g] * (1 - γᵍ[g]) * ρᴵᵍ[g, i, t] + (1 - ωᵍ[g]) * χᵍ[g] * ρᴴᵍ[g, i, t] + ρᴿᵍ[g, i, t]

            # Reset values
            τᵢᵍ[g, i] = 0.
            Pᵢᵍ[g, i] = 0.
        end
    end
end

"""
    compute_P!(Pᵢᵍ::Array{Float64, 2},
                    Pᵢᴬᵍ::Array{Float64, 2},
                    Pᵢᴵᵍ::Array{Float64, 2},
                    Sᵢᵍ::Array{Float64, 2},
                    pᵍ::Array{Float64, 1},
                    ρˢᵍ::Array{Float64, 3},
                    ρᴬᵍ::Array{Float64, 3},
                    ρᴵᵍ::Array{Float64, 3},
                    Qᵢᵍ::Array{Float64, 3},
                    nᵢᵍ_eff::Array{Float64, 2},
                    mobilityᵍ::Array{Float64, 2},
                    normᵍ::Array{Float64, 2},
                    βᴬ::Float64,
                    βᴵ::Float64,
                    νᵍ::Array{Float64, 1},
                    edgelist::Array{Int64, 2},
                    Rᵢⱼ::Array{Float64, 1},
                    C::Array{Float64, 2},
                    M::Int64,
                    G::Int64,
                    t::Int64)

Function to compute ``P_i^g(t)`` and ``Q_i^g(t)`` as described in the referenced paper.
"""
function compute_P!(Pᵢᵍ::Array{Float64, 2},
                    Pᵢᴬᵍ::Array{Float64, 2},
                    Pᵢᴵᵍ::Array{Float64, 2},
                    Sᵢᵍ::Array{Float64, 2},
                    pᵍ::Array{Float64, 1},
                    ρˢᵍ::Array{Float64, 3},
                    ρᴬᵍ::Array{Float64, 3},
                    ρᴵᵍ::Array{Float64, 3},
                    Qᵢᵍ::Array{Float64, 3},
                    nᵢᵍ_eff::Array{Float64, 2},
                    mobilityᵍ::Array{Float64, 2},
                    normᵍ::Array{Float64, 2},
                    βᴬ::Float64,
                    βᴵ::Float64,
                    νᵍ::Array{Float64, 1},
                    edgelist::Array{Int64, 2},
                    Rᵢⱼ::Array{Float64, 1},
                    C::Array{Float64, 2},
                    M::Int64,
                    G::Int64,
                    t::Int64)

    # Init. aux variables
    Sᵢᵍ .= 0.
    Pᵢᴬᵍ .= 0.
    Pᵢᴵᵍ .= 0.

    for indx_e in 1:length(Rᵢⱼ)
        i = edgelist[indx_e, 1]
        j = edgelist[indx_e, 2]

        # Get effective S, A and I
        for g in 1:G
            nˢᵍ_ij = ρˢᵍ[g, i, t] * mobilityᵍ[g, indx_e]
            Sᵢᵍ[g, j] +=  nˢᵍ_ij / nᵢᵍ_eff[g, j]
            @simd for h in 1:G
                nᴬᵍ_ij = ρᴬᵍ[h, i, t] * mobilityᵍ[h, indx_e]
                nᴵᵍ_ij = ρᴵᵍ[h, i, t] * mobilityᵍ[h, indx_e]
                Pᵢᴬᵍ[g, j] += C[g, h] * nᴬᵍ_ij / nᵢᵍ_eff[h, j]
                Pᵢᴵᵍ[g, j] += C[g, h] * νᵍ[h] * nᴵᵍ_ij / nᵢᵍ_eff[h, j]
            end
        end
    end

    # Get P and effective ρ
    for i in 1:M
        @simd for g in 1:G
            Pᵢᵍ[g, i] = 1 - (1 - βᴬ)^(normᵍ[g, i] * Pᵢᴬᵍ[g, i]) * (1 - βᴵ)^(normᵍ[g, i] * Pᵢᴵᵍ[g, i])
            Sᵢᵍ[g, i] = normᵍ[g, i] * Sᵢᵍ[g, i]
        end
    end

    # Compute Q to get the effective R
    for indx_e in 1:length(Rᵢⱼ)
        i = edgelist[indx_e, 1]
        j = edgelist[indx_e, 2]

        for g in 1:G
            for h in 1:G
                Qᵢᵍ[g, i, t] += C[g, h] * Sᵢᵍ[h, j] * (pᵍ[g] * Rᵢⱼ[indx_e] + (1 - pᵍ[g]) * (i == j ? 1. : 0.))
            end
        end
    end
end

"""
    run_epidemic_spreading_mmca!(epi_param::Epidemic_Params,
                                 population::Population_Params,
                                 M::Int64,
                                 G::Int64,
                                 edgelist::Array{Int64,2},
                                 Rᵢⱼ::Array{Float64,1},
                                 C::Array{Float64,2},
                                 pᵍ::Array{Float64,1},
                                 timesteps::Int64;
                                 κ₀::Float64 = 0.0,
                                 ϕ::Float64 = 1.0,
                                 tᶜ::Int64 = -1,
                                 t₀::Int64 = 1,
                                 verbose::Bool = false)

This functions computes the evolution of the epidemic parameters over time,
updating the variables stored in epi_param. It also provides, through optional
arguments, the application of a containmnet on specific date.

# Arguments

- `epi_param::Epidemic_Params`: Structure that contains all related epidemic
  parameters (see `Epidemic_Params` documentation for more infomation), including
  the ones that will be updated.
- `population::Population_Params`: Structure that contains all the parameters
  related with the population, (see `Population_Params` documentation for more
  information).
- `M::Int64`: Number of patches.
- `G::Int64`: Number of stratifications of the population.
- `edgelist::Array{Float64, 2}`: ``M \\times 2`` matrix containing the directed
  edgelist between patches. The IDs of the patches have to go from 1 to M
- `Rᵢⱼ::Array{Float64, 1}`: Vector containing the transition probabilities for
  each edge in the edgelist.
- `C::Array{Float64, 2}`: ``G \\times G`` contact matrix containing the
  probabilities of contact between different population strata.
- `pᵍ::Array{Float64, 1}`: Degree of mobility of each strata ranged between 0 and 1.
- `timesteps::Int64`: Number of timesteps to evolve the system.

## Optional

- `κ⁰::Float64 = 0.0`: Confinement factor.
- `ϕ::Float64 = 1.0`: Permeability factor.
- `tᶜ::Int64 = -1`: Timestep where the confinment measures will be applied.
- `t₀::Int64 = 1`: Initial timestep (useful to define multiple containment
  strategies).
- `verbose::Bool = false`: if it is set to `true` it prints useful information
  about the evolution of the epidemic process.
"""
function run_epidemic_spreading_mmca!(epi_param::Epidemic_Params,
                                      population::Population_Params,
                                      M::Int64,
                                      G::Int64,
                                      edgelist::Array{Int64, 2},
                                      Rᵢⱼ::Array{Float64, 1},
                                      C::Array{Float64, 2},
                                      pᵍ::Array{Float64, 1},
                                      timesteps::Int64;
                                      κ₀::Float64 = 0.,
                                      ϕ::Float64 = 1.,
                                      tᶜ::Int64 = -1,
                                      t₀::Int64 = 1,
                                      verbose::Bool = false
                                      )

    # Initialize τᵢ (Π = (1 - p) P + pτ) and Pᵢ for markov chain
    τᵢᵍ = zeros(Float64, G, M)
    Pᵢᵍ = zeros(Float64, G, M)

    # Auxiliar arrays to compute P (avoid the allocation of additional memory)
    Pᵢᴬᵍ = zeros(Float64, G, M)
    Pᵢᴵᵍ = zeros(Float64, G, M)
    Sᵢᵍ = zeros(Float64, G, M)

    # Start loop for time evoluiton
    @inbounds for t in t₀:(timesteps - 1)
        update_prob!(Pᵢᵍ, Pᵢᴬᵍ, Pᵢᴵᵍ, Sᵢᵍ, τᵢᵍ, epi_param, population, M, G, edgelist, Rᵢⱼ, C, κ₀, ϕ, pᵍ, t, tᶜ)

        if verbose
            players = sum((epi_param.ρˢᵍ[:, :, t] .+ epi_param.ρᴱᵍ[:, :, t] .+ epi_param.ρᴬᵍ[:, :, t] .+
                           epi_param.ρᴵᵍ[:, :, t] .+ epi_param.ρᴿᵍ[:, :, t] .+ epi_param.ρᴴᵍ[:, :, t] .+
                           epi_param.ρᴰᵍ[:, :, t]) .* population.nᵢᵍ[:, :])
            cases = sum((epi_param.ρᴿᵍ[:, :, t] .+ epi_param.ρᴴᵍ[:, :, t] .+ epi_param.ρᴰᵍ[:, :, t]) .* population.nᵢᵍ[:, :])
            hospital = sum(epi_param.ρᴴᵍ[:, :, t] .* population.nᵢᵍ[:, :])
            deaths = sum(epi_param.ρᴰᵍ[:, :, t] .* population.nᵢᵍ[:, :])
            infected = sum(epi_param.ρᴵᵍ[:, :, t] .* population.nᵢᵍ[:, :] + epi_param.ρᴬᵍ[:, :, t] .* population.nᵢᵍ[:, :])

            @printf("Time: %d, players: %.2f, infected: %.2f, cases: %.2f, hospital: %.2f, deaths: %.2f\n", t, players, infected, cases, hospital, deaths)
        end
    end
end

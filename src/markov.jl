"""
    update_prob!(Pᵢᵍ::Array{Float64, 2},
                 Pᵢᴬᵍ::Array{Float64, 2},
                 Pᵢᴵᵍ::Array{Float64, 2},
                 Sᵢᵍ::Array{Float64, 2},
                 τᵢᵍ::Array{Float64, 2},
                 epi_params::Epidemic_Params,
                 population::Population_Params,
                 κ₀::Float64,
                 ϕ::Float64,
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
                      κ₀::Float64,
                      ϕ::Float64,
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
    G = population.G
    M = population.M
    pᵍ = population.pᵍ
    C = population.C
    edgelist = population.edgelist
    Rᵢⱼ = population.Rᵢⱼ

    # Intervention at time tᶜ
    if tᶜ == t
        pᵍ[:] .= (1 - κ₀) .* pᵍ
        population.kᵍ[:] .= [population.σ - 1, (1 - κ₀) * population.kᵍ[2] + κ₀ * (population.σ - 1), population.σ - 1]
        update_population_params!(population)
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

Compute ``P_i^g(t)`` and ``Q_i^g(t)`` as described in the referenced paper.
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
    print_status(epi_params::Epidemic_Params,
                 population::Population_Params,
                 t::Int64)

Print the status of the epidemic spreading.
"""
function print_status(epi_params::Epidemic_Params,
                      population::Population_Params,
                      t::Int64)
    players = sum((epi_params.ρˢᵍ[:, :, t] .+ epi_params.ρᴱᵍ[:, :, t] .+ epi_params.ρᴬᵍ[:, :, t] .+
                   epi_params.ρᴵᵍ[:, :, t] .+ epi_params.ρᴿᵍ[:, :, t] .+ epi_params.ρᴴᵍ[:, :, t] .+
                   epi_params.ρᴰᵍ[:, :, t]) .* population.nᵢᵍ[:, :])
    infected = sum(epi_params.ρᴵᵍ[:, :, t] .* population.nᵢᵍ[:, :] + epi_params.ρᴬᵍ[:, :, t] .* population.nᵢᵍ[:, :])
    cases = sum((epi_params.ρᴿᵍ[:, :, t] .+ epi_params.ρᴴᵍ[:, :, t] .+ epi_params.ρᴰᵍ[:, :, t]) .* population.nᵢᵍ[:, :])
    icus = sum(epi_params.ρᴴᵍ[:, :, t] .* population.nᵢᵍ[:, :])
    deaths = sum(epi_params.ρᴰᵍ[:, :, t] .* population.nᵢᵍ[:, :])

    @printf("Time: %d, players: %.2f, infected: %.2f, cases: %.2f, icus: %.2f, deaths: %.2f\n", t, players, infected, cases, icus, deaths)
end


"""
    run_epidemic_spreading_mmca!(epi_params::Epidemic_Params,
                                 population::Population_Params;
                                 t₀::Int64 = 1,
                                 tᶜ::Int64 = -1,
                                 κ₀::Float64 = 0.0,
                                 ϕ::Float64 = 1.0,
                                 verbose::Bool = false)

Computes the evolution of the epidemic spreading over time, updating the variables stored in epi_params. It also provides, through optional arguments, the application of a containmnet on a specific date.

# Arguments

- `epi_params::Epidemic_Params`: Structure that contains all epidemic parameters and the epidemic spreading information.
- `population::Population_Params`: Structure that contains all the parameters related with the population.

## Optional

- `t₀::Int64 = 1`: Initial timestep (useful to define multiple containment strategies).
- `tᶜ::Int64 = -1`: Timestep where the confinment measures is applied, -1 for no containment.
- `κ⁰::Float64 = 0.0`: Confinement factor.
- `ϕ::Float64 = 1.0`: Permeability factor.
- `verbose::Bool = false`: If `true`, prints useful information about the evolution of the epidemic process.
"""
function run_epidemic_spreading_mmca!(epi_params::Epidemic_Params,
                                      population::Population_Params;
                                      t₀::Int64 = 1,
                                      tᶜ::Int64 = -1,
                                      κ₀::Float64 = 0.,
                                      ϕ::Float64 = 1.,
                                      verbose::Bool = false)

    G = population.G
    M = population.M
    T = epi_params.T

    # Initialize τᵢ (Π = (1 - p) P + pτ) and Pᵢ for markov chain
    τᵢᵍ = zeros(Float64, G, M)
    Pᵢᵍ = zeros(Float64, G, M)

    # Auxiliar arrays to compute P (avoid the allocation of additional memory)
    Pᵢᴬᵍ = zeros(Float64, G, M)
    Pᵢᴵᵍ = zeros(Float64, G, M)
    Sᵢᵍ = zeros(Float64, G, M)

    # Initial state
    if verbose
        print_status(epi_params, population, t₀)
    end

    # Start loop for time evoluiton
    for t in t₀:(T - 1)
        update_prob!(Pᵢᵍ, Pᵢᴬᵍ, Pᵢᴵᵍ, Sᵢᵍ, τᵢᵍ, epi_params, population, κ₀, ϕ, t, tᶜ)

        if verbose
            print_status(epi_params, population, t + 1)
        end
    end
end

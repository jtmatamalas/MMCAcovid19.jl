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
                 δ::Float64,
                 t::Int64,
                 tᶜ::Int64)

Updates the probabilities of the model using the equations described in the
paper.
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
                      δ::Float64,
                      t::Int64,
                      tᶜ::Int64)

    # Shortcuts to parameters
    ηᵍ = epi_params.ηᵍ
    αᵍ = epi_params.αᵍ
    μᵍ = epi_params.μᵍ
    θᵍ = epi_params.θᵍ
    γᵍ = epi_params.γᵍ
    ζᵍ = epi_params.ζᵍ
    λᵍ = epi_params.λᵍ
    ωᵍ = epi_params.ωᵍ
    ψᵍ = epi_params.ψᵍ
    χᵍ = epi_params.χᵍ
    ρˢᵍ = epi_params.ρˢᵍ
    ρᴱᵍ = epi_params.ρᴱᵍ
    ρᴬᵍ = epi_params.ρᴬᵍ
    ρᴵᵍ = epi_params.ρᴵᵍ
    ρᴾᴴᵍ = epi_params.ρᴾᴴᵍ
    ρᴾᴰᵍ = epi_params.ρᴾᴰᵍ
    ρᴴᴿᵍ = epi_params.ρᴴᴿᵍ
    ρᴴᴰᵍ = epi_params.ρᴴᴰᵍ
    ρᴰᵍ = epi_params.ρᴰᵍ
    ρᴿᵍ = epi_params.ρᴿᵍ
    CHᵢᵍ = epi_params.CHᵢᵍ
    G = population.G
    M = population.M
    pᵍ_eff = population.pᵍ_eff
    C = population.C
    edgelist = population.edgelist
    Rᵢⱼ = population.Rᵢⱼ
    kᵍ_h = population.kᵍ_h
    kᵍ_hw = population.kᵍ_h .+ population.kᵍ_w

    # Intervention at time tᶜ
    if tᶜ == t
        pᵍ_eff[:] .= (1 - κ₀) .* population.pᵍ

        population.kᵍ_eff .= kᵍ_h * κ₀ .+ kᵍ_hw * (1 - δ) * (1 - κ₀)
        # elder keep home contacts during confinement
        population.kᵍ_eff[G] = kᵍ_h[G]

        update_population_params!(population)
    end

    # Get P and compute Q
    compute_P!(Pᵢᵍ, Pᵢᴬᵍ, Pᵢᴵᵍ, Sᵢᵍ, pᵍ_eff, ρˢᵍ, ρᴬᵍ, ρᴵᵍ,
               epi_params.Qᵢᵍ, population.nᵢᵍ_eff, population.mobilityᵍ,
               population.normᵍ, epi_params.βᴬ[1], epi_params.βᴵ[1],
               edgelist, Rᵢⱼ, C, M, G, t)

    # Comptue τᵢᵍ
    @inbounds for indx_e in 1:length(Rᵢⱼ)
        i = edgelist[indx_e, 1]
        j = edgelist[indx_e, 2]

        @simd for g in 1:G
            τᵢᵍ[g, i] += Rᵢⱼ[indx_e] * Pᵢᵍ[g, j]
        end
    end

    # Update probabilities
    @inbounds for i in 1:M

        # Compute secure households
        CHᵢ = 0.0
        if tᶜ == t
            @simd for g in 1:G
                CHᵢ += (ρˢᵍ[g, i, t] + ρᴾᴴᵍ[g, i, t] + ρᴾᴰᵍ[g, i, t] +
                        ρᴴᴿᵍ[g, i, t] + ρᴴᴰᵍ[g, i, t] + ρᴰᵍ[g, i, t] +
                        ρᴿᵍ[g, i, t] + CHᵢᵍ[g, i]) * population.nᵢᵍ[g, i]
            end
            CHᵢ = (1 - ϕ) * κ₀ * (CHᵢ / population.nᵢ[i]) ^ population.σ
        end

        # Update compartmental probabilities
        @simd for g in 1:G
            if tᶜ == t
                ρˢᵍ[g, i, t] += CHᵢᵍ[g, i]
            end

            Πᵢᵍ = (1 - pᵍ_eff[g]) * Pᵢᵍ[g, i] + pᵍ_eff[g] * τᵢᵍ[g, i]

            # Spreading

            ρˢᵍ[g, i, t + 1] = (1 - Πᵢᵍ) * (1 - CHᵢ) * ρˢᵍ[g, i, t]

            ρᴱᵍ[g, i, t + 1] = (1 - ηᵍ[g]) * ρᴱᵍ[g, i, t] +
                Πᵢᵍ * (1 - CHᵢ) * ρˢᵍ[g, i, t]

            ρᴬᵍ[g, i, t + 1] = (1 - αᵍ[g]) * ρᴬᵍ[g, i, t] +
                ηᵍ[g] * ρᴱᵍ[g, i, t]

            ρᴵᵍ[g, i, t + 1] = (1 - μᵍ[g]) * ρᴵᵍ[g, i, t] +
                αᵍ[g] * ρᴬᵍ[g, i, t]

            ρᴾᴴᵍ[g, i, t + 1] = (1 - λᵍ[g]) * ρᴾᴴᵍ[g, i, t] +
                μᵍ[g] * (1 - θᵍ[g]) * γᵍ[g] * ρᴵᵍ[g, i, t]

            ρᴾᴰᵍ[g, i, t + 1] = (1 - ζᵍ[g]) * ρᴾᴰᵍ[g, i, t] +
                μᵍ[g] * θᵍ[g] * ρᴵᵍ[g, i, t]

            ρᴴᴿᵍ[g, i, t + 1] = (1 - χᵍ[g]) * ρᴴᴿᵍ[g, i, t] +
                λᵍ[g] * (1 - ωᵍ[g]) * ρᴾᴴᵍ[g, i, t]

            ρᴴᴰᵍ[g, i, t + 1] = (1 - ψᵍ[g]) * ρᴴᴰᵍ[g, i, t] +
                λᵍ[g] * ωᵍ[g] * ρᴾᴴᵍ[g, i, t]

            ρᴿᵍ[g, i, t + 1] = ρᴿᵍ[g, i, t] + χᵍ[g] * ρᴴᴿᵍ[g, i, t] +
                μᵍ[g] * (1 - θᵍ[g]) * (1 - γᵍ[g]) * ρᴵᵍ[g, i , t]

            ρᴰᵍ[g, i, t + 1] = ρᴰᵍ[g, i, t] + ζᵍ[g] * ρᴾᴰᵍ[g, i, t] +
                ψᵍ[g] * ρᴴᴰᵍ[g, i, t]

            # Reset values
            τᵢᵍ[g, i] = 0.
            Pᵢᵍ[g, i] = 0.

            if tᶜ == t
                aux = ρˢᵍ[g, i, t]
                ρˢᵍ[g, i, t] -= CHᵢᵍ[g, i]
                CHᵢᵍ[g, i] = CHᵢ * aux
            end
        end
    end
end


"""
    compute_P!(Pᵢᵍ::Array{Float64, 2},
                    Pᵢᴬᵍ::Array{Float64, 2},
                    Pᵢᴵᵍ::Array{Float64, 2},
                    Sᵢᵍ::Array{Float64, 2},
                    pᵍ_eff::Array{Float64, 1},
                    ρˢᵍ::Array{Float64, 3},
                    ρᴬᵍ::Array{Float64, 3},
                    ρᴵᵍ::Array{Float64, 3},
                    Qᵢᵍ::Array{Float64, 3},
                    nᵢᵍ_eff::Array{Float64, 2},
                    mobilityᵍ::Array{Float64, 2},
                    normᵍ::Array{Float64, 2},
                    βᴬ::Float64,
                    βᴵ::Float64,
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
                    pᵍ_eff::Array{Float64, 1},
                    ρˢᵍ::Array{Float64, 3},
                    ρᴬᵍ::Array{Float64, 3},
                    ρᴵᵍ::Array{Float64, 3},
                    Qᵢᵍ::Array{Float64, 3},
                    nᵢᵍ_eff::Array{Float64, 2},
                    mobilityᵍ::Array{Float64, 2},
                    normᵍ::Array{Float64, 2},
                    βᴬ::Float64,
                    βᴵ::Float64,
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

    @inbounds for indx_e in 1:length(Rᵢⱼ)
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
                Pᵢᴵᵍ[g, j] += C[g, h] * nᴵᵍ_ij / nᵢᵍ_eff[h, j]
            end
        end
    end

    # Get P and effective ρ
    @inbounds for i in 1:M
        @simd for g in 1:G
            Pᵢᵍ[g, i] = 1 - (1 - βᴬ)^(normᵍ[g, i] * Pᵢᴬᵍ[g, i]) *
                (1 - βᴵ)^(normᵍ[g, i] * Pᵢᴵᵍ[g, i])
        end
    end

    # Compute Q to get the effective R
    @inbounds for indx_e in 1:length(Rᵢⱼ)
        i = edgelist[indx_e, 1]
        j = edgelist[indx_e, 2]

        for g in 1:G
            @simd for h in 1:G
                Qᵢᵍ[g, i, t] += normᵍ[g, i] * C[g, h] * Sᵢᵍ[h, j] *
                    (pᵍ_eff[g] * Rᵢⱼ[indx_e] +
                     (1 - pᵍ_eff[g]) * (i == j ? 1. : 0.))
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

    players  = sum((epi_params.ρᴾᴰᵍ[:, :, t] .+
                    epi_params.ρˢᵍ[:, :, t] .+
                    epi_params.ρᴱᵍ[:, :, t] .+
                    epi_params.ρᴬᵍ[:, :, t] .+
                    epi_params.ρᴵᵍ[:, :, t] .+
                    epi_params.ρᴾᴴᵍ[:, :, t] .+
                    epi_params.ρᴴᴰᵍ[:, :, t] .+
                    epi_params.ρᴴᴿᵍ[:, :, t] .+
                    epi_params.ρᴿᵍ[:, :, t] .+
                    epi_params.ρᴰᵍ[:, :, t]) .* population.nᵢᵍ[:, :])

    infected = sum(epi_params.ρᴵᵍ[:, :, t] .* population.nᵢᵍ[:, :] +
                   epi_params.ρᴬᵍ[:, :, t] .* population.nᵢᵍ[:, :])

    cases    = sum((epi_params.ρᴾᴰᵍ[:, :, t] .+
                    epi_params.ρᴾᴴᵍ[:, :, t] .+
                    epi_params.ρᴴᴰᵍ[:, :, t] .+
                    epi_params.ρᴴᴿᵍ[:, :, t] .+
                    epi_params.ρᴿᵍ[:, :, t] .+
                    epi_params.ρᴰᵍ[:, :, t]) .* population.nᵢᵍ[:, :])

    icus     = sum((epi_params.ρᴴᴿᵍ[:, :, t] .+
                    epi_params.ρᴴᴰᵍ[:, :, t]) .* population.nᵢᵍ[:, :])

    deaths   = sum(epi_params.ρᴰᵍ[:, :, t] .* population.nᵢᵍ[:, :])

    @printf("Time: %d, players: %.2f, infected: %.2f, cases: %.2f, icus: %.2f, deaths: %.2f\n",
            t, players, infected, cases, icus, deaths)
end


"""
    run_epidemic_spreading_mmca!(epi_params::Epidemic_Params,
                                 population::Population_Params;
                                 tᶜ::Int64 = -1,
                                 κ₀::Float64 = 0.0,
                                 ϕ::Float64 = 1.0,
                                 δ::Float64 = 0.0,
                                 t₀::Int64 = 1,
                                 verbose::Bool = false)

Computes the evolution of the epidemic spreading over time, updating the
variables stored in epi_params. It also provides, through optional arguments,
the application of a containmnet on a specific date.

# Arguments

- `epi_params::Epidemic_Params`: Structure that contains all epidemic parameters
  and the epidemic spreading information.
- `population::Population_Params`: Structure that contains all the parameters
  related with the population.

## Optional

- `tᶜ::Int64 = -1`: Timestep of application of containment, or out of timesteps range
  value for no containment.
- `κ⁰::Float64 = 0.0`: Mobility reduction.
- `ϕ::Float64 = 1.0`: Permeability of confined households.
- `δ::Float64 = 0.0`: Social Distancing.
- `t₀::Int64 = 1`: Initial timestep.
- `verbose::Bool = false`: If `true`, prints useful information about the
  evolution of the epidemic process.
"""
function run_epidemic_spreading_mmca!(epi_params::Epidemic_Params,
                                      population::Population_Params;
                                      tᶜ::Int64 = -1,
                                      κ₀::Float64 = 0.0,
                                      ϕ::Float64 = 1.0,
                                      δ::Float64 = 0.0,
                                      t₀::Int64 = 1,
                                      verbose::Bool = false)
    G = population.G
    M = population.M

    # Initialize τᵢ (Π = (1 - p) P + pτ) and Pᵢ for markov chain
    τᵢᵍ = zeros(Float64, G, M)
    Pᵢᵍ = zeros(Float64, G, M)

    # Auxiliar arrays to compute P (avoid the allocation of additional memory)
    Pᵢᴬᵍ = zeros(Float64, G, M)
    Pᵢᴵᵍ = zeros(Float64, G, M)
    Sᵢᵍ = zeros(Float64, G, M)

    run_epidemic_spreading_mmca!(epi_params, population, [tᶜ],
                                 [κ₀], [ϕ], [δ], t₀ = t₀, verbose = verbose)
end


"""
    run_epidemic_spreading_mmca!(epi_params::Epidemic_Params,
                                 population::Population_Params,
                                 tᶜs::Array{Int64, 1},
                                 κ₀s::Array{Float64, 1},
                                 ϕs::Array{Float64, 1},
                                 δs::Array{Float64, 1};
                                 t₀::Int64 = 1,
                                 verbose::Bool = false)

Computes the evolution of the epidemic spreading over time, updating the
variables stored in epi_params. It provides the option of the application
of multiple different containmnets at specific dates.

# Arguments

- `epi_params::Epidemic_Params`: Structure that contains all epidemic parameters
  and the epidemic spreading information.
- `population::Population_Params`: Structure that contains all the parameters
  related with the population.
- `tᶜs::Array{Int64, 1}`: List of timesteps of application of containments.
- `κ⁰s::Array{Float64, 1}`: List of mobility reductions.
- `ϕs::Array{Float64, 1}`: List of permeabilities of confined households.
- `δs::Array{Float64, 1}`: List of social distancings.

## Optional

- `t₀::Int64 = 1`: Initial timestep.
- `verbose::Bool = false`: If `true`, prints useful information about the
  evolution of the epidemic process.
"""
function run_epidemic_spreading_mmca!(epi_params::Epidemic_Params,
                                      population::Population_Params,
                                      tᶜs::Array{Int64, 1},
                                      κ₀s::Array{Float64, 1},
                                      ϕs::Array{Float64, 1},
                                      δs::Array{Float64, 1};
                                      t₀::Int64 = 1,
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

    i = 1

    ## Start loop for time evoluiton
    @inbounds for t in t₀:(T - 1)
        update_prob!(Pᵢᵍ, Pᵢᴬᵍ, Pᵢᴵᵍ, Sᵢᵍ, τᵢᵍ, epi_params, population,
                        κ₀s[i], ϕs[i], δs[i], t, tᶜs[i])

        if t == tᶜs[i] && i < length(tᶜs)
            i += 1
        end

        if verbose
            print_status(epi_params, population, t + 1)
        end
    end
end

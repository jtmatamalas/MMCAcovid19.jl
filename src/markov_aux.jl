#------------------------------------------------------
# EPIDEMIC PARAMS RELATED FUNCTIONS
#------------------------------------------------------

"""
    Epidemic_Params

Struct that contains the parameters related with the epidemic parameters and compartmental evolution.

# Fields

All the parameters contained in this structure are probabilities ranged between 0 and 1.

# Epidemic parameters

- `βᴬ::Float64`: Infectivity of asymptomatic.
- `βᴵ::Float64`: Infectivity of infected.
- `ηᵍ::Array{Float64, 1}`: Latent rate for each strata.
- `αᵍ::Array{Float64, 1}`: Asymptomatic infectious rate for each strata.
- `νᵍ::Array{Float64, 1}`: Isolation factor for each strata.
- `μᵍ::Array{Float64, 1}`: Escape rate for each strata.
- `γᵍ::Array{Float64, 1}`: Fraction of cases requiring ICU for each strata.
- `ωᵍ::Array{Float64, 1}`: Fatality rate of ICU patients for each strata.
- `ψᵍ::Array{Float64, 1}`: Death rate for each strata.
- `χᵍ::Array{Float64, 1}`: ICU discharge rate for each strata.
- `T::Int64`: Number of epidemic timesteps.

# Compartmental evolution

- `ρˢᵍ::Array{Float64, 3}`: Matrix of size ``G \\times M \\times T`` containing infomation about the evolution of fraction of suceptible individuals for each strata and patch.
- `ρᴱᵍ::Array{Float64, 3}`: Matrix of size ``G \\times M \\times T`` containing infomation about the evolution of fraction of exposed individuals for each strata and patch.
- `ρᴬᵍ::Array{Float64, 3}`: Matrix of size ``G \\times M \\times T`` containing infomation about the evolution of fraction of asymptomatic individuals for each strata and patch.
- `ρᴵᵍ::Array{Float64, 3}`: Matrix of size ``G \\times M \\times T`` containing infomation about the evolution of fraction of symptomatic individuals for each strata and patch.
- `ρᴴᵍ::Array{Float64, 3}`: Matrix of size ``G \\times M \\times T`` containing infomation about the evolution of fraction of ICU patients individuals for each strata and patch.
- `ρᴰᵍ::Array{Float64, 3}`: Matrix of size ``G \\times M \\times T`` containing infomation about the evolution of fraction of death individuals for each strata and patch.
- `ρᴿᵍ::Array{Float64, 3}`: Matrix of size ``G \\times M \\times T`` containing infomation about the evolution of fraction of recovered individuals for each strata and patch.

# Auxiliar

- `CHᵢᵍ::Array{Float64, 2}`: Fraction of securely confined individuals for each strata and patch.
- `Qᵢᵍ::Array{Float64, 3}`: Suceptible contacts available for each strata on a given patch.

# Constructor

Use outer constructor [`Epidemic_Params`](#MMCAcovid19.Epidemic_Params-Tuple{Float64,Float64,Array{Float64,1},Array{Float64,1},Array{Float64,1},Array{Float64,1},Array{Float64,1},Array{Float64,1},Array{Float64,1},Array{Float64,1},Int64,Int64,Int64}) for a proper initialization of this struct.
"""
struct Epidemic_Params
    #Epidemic parameters
    βᴬ::Float64
    βᴵ::Float64
    ηᵍ::Array{Float64, 1}
    αᵍ::Array{Float64, 1}
    νᵍ::Array{Float64, 1}
    μᵍ::Array{Float64, 1}
    γᵍ::Array{Float64, 1}
    ωᵍ::Array{Float64, 1}
    ψᵍ::Array{Float64, 1}
    χᵍ::Array{Float64, 1}
    T::Int64

    # Compartments evolution
    ρˢᵍ::Array{Float64, 3}
    ρᴱᵍ::Array{Float64, 3}
    ρᴬᵍ::Array{Float64, 3}
    ρᴵᵍ::Array{Float64, 3}
    ρᴴᵍ::Array{Float64, 3}
    ρᴰᵍ::Array{Float64, 3}
    ρᴿᵍ::Array{Float64, 3}
    CHᵢᵍ::Array{Float64, 2}

    # R_t related arrays
    Qᵢᵍ::Array{Float64, 3}
end


"""
    Epidemic_Params(βᴬ::Float64,
                    βᴵ::Float64,
                    ηᵍ::Array{Float64, 1},
                    αᵍ::Array{Float64, 1},
                    νᵍ::Array{Float64, 1},
                    μᵍ::Array{Float64, 1},
                    γᵍ::Array{Float64, 1},
                    ωᵍ::Array{Float64, 1},
                    ψᵍ::Array{Float64, 1},
                    χᵍ::Array{Float64, 1},
                    M::Int64,
                    G::Int64,
                    T::Int64)

Constructor of the struct  [`Epidemic_Params`](@ref).

# Arguments

- `βᴬ::Float64`: Infectivity of asymptomatic.
- `βᴵ::Float64`: Infectivity of infected.
- `ηᵍ::Array{Float64, 1}`: Vector of size ``G`` with latent rates for each strata.
- `αᵍ::Array{Float64, 1}`: Vector of size ``G`` with asymptomatic infectious rates for each strata.
- `νᵍ::Array{Float64, 1}`: Vector of size ``G`` with isolation factors for each strata.
- `μᵍ::Array{Float64, 1}`: Vector of size ``G`` with escape rates for each strata.
- `γᵍ::Array{Float64, 1}`: Vector of size ``G`` with fractions of cases requiring ICU for each strata.
- `ωᵍ::Array{Float64, 1}`: Vector of size ``G`` with fatality rates of ICU patients for each strata.
- `ψᵍ::Array{Float64, 1}`: Vector of size ``G`` with death rates for each strata.
- `χᵍ::Array{Float64, 1}`: Vector of size ``G`` with ICU discharge rates for each strata.
- `M::Int64`: Number of patches.
- `G::Int64`: Number of strata.
- `T::Int64`: Number of epidemic timesteps.
"""
function Epidemic_Params(βᴬ::Float64,
                         βᴵ::Float64,
                         ηᵍ::Array{Float64, 1},
                         αᵍ::Array{Float64, 1},
                         νᵍ::Array{Float64, 1},
                         μᵍ::Array{Float64, 1},
                         γᵍ::Array{Float64, 1},
                         ωᵍ::Array{Float64, 1},
                         ψᵍ::Array{Float64, 1},
                         χᵍ::Array{Float64, 1},
                         M::Int64,
                         G::Int64,
                         T::Int64)

    # Allocate memory for simulations
    ρˢᵍ = ones(Float64, G, M, T)
    ρᴱᵍ = zeros(Float64, G, M, T)
    ρᴬᵍ = zeros(Float64, G, M, T)
    ρᴵᵍ = zeros(Float64, G, M, T)
    ρᴴᵍ = zeros(Float64, G, M, T)
    ρᴰᵍ = zeros(Float64, G, M, T)
    ρᴿᵍ = zeros(Float64, G, M, T)
    CHᵢᵍ = zeros(Float64, G, M)
    Qᵢᵍ = zeros(Float64, G, M, T)

    return Epidemic_Params(βᴬ, βᴵ, copy(ηᵍ), copy(αᵍ), copy(νᵍ), copy(μᵍ), copy(γᵍ), copy(ωᵍ), copy(ψᵍ), copy(χᵍ), T,
                           ρˢᵍ, ρᴱᵍ, ρᴬᵍ, ρᴵᵍ, ρᴴᵍ, ρᴰᵍ, ρᴿᵍ, CHᵢᵍ, Qᵢᵍ)
end


"""
    reset_epidemic_params!(epi_params::Epidemic_Params)

Reset the ρ's to reuse the structure and avoid additional allocations.

# Arguments
- `epi_params::Epidemic_Params`: structure to reset.
"""
function reset_epidemic_params!(epi_params::Epidemic_Params)
    epi_params.ρˢᵍ .= 1.
    epi_params.ρᴱᵍ .= 0.
    epi_params.ρᴬᵍ .= 0.
    epi_params.ρᴵᵍ .= 0.
    epi_params.ρᴴᵍ .= 0.
    epi_params.ρᴰᵍ .= 0.
    epi_params.ρᴿᵍ .= 0.
    epi_params.CHᵢᵍ .= 0.

    # Rt structures
    epi_params.Qᵢᵍ .= 0.
end


#------------------------------------------------------
# PATCH AND POPULATION RELATED FUNCTIONS
#------------------------------------------------------

"""
    Population_Params

Struct that contains the parameters related with geographical, population and mobility data.

# Fields

- `M::Int64`: Number of patches.
- `G::Int64`: Number of population strata.
- `nᵢ:Array{Float64, 1}`: Population of each patch.
- `nᵢᵍ:Array{Float64, 2}`: Population of each strata on each patch.
- `nᵢ_eff:Array{Float64, 1}`: Effective population on each patch after taking into account mobility.
- `nᵢᵍ_eff:Array{Float64, 1}`: Effective population of each strata on each patch after taking into account mobility.
- `N::Int64`: Total population.
- `Nᵍ::Array{Int64, 1}`: Total population of each strata.
- `pᵍ::Array{Float64, 1}`: Vector with the degree of mobility of each strata.
- `edgelist::Array{Int64, 2}`: Matrix with the directed edgelist between patches, where ``L`` is the number of edges.
- `Rᵢⱼ::Array{Float64, 1}`: Vector with the transition probabilities for each edge in the edgelist.
- `mobilityᵍ::Array{Float64, 1}`: Effective mobility of each strata between patches (used to optimize performance).
- `kᵍ::Array{Float64, 1}`: Average number of contacts of each strata.
- `C::Array{Float64, 2}`: Matrix with the probability of contact between different stratifications.
- `zᵍ::Array{Float64, 1}`: Nomalization factor each strata.
- `normᵍ::Array{Float64, 2}`: Normalization of each strata (used to optimize performance).
- `sᵢ::Array{Float64, 1}`: Surface of each patch.
- `ξ::Float64`: Densisty factor.
- `σ::Float64`: Average household size.

# Constructor

Use outer constructor [`Population_Params`](#MMCAcovid19.Population_Params-Tuple{Int64,Int64,Array{Float64,2},Array{Float64,1},Array{Float64,2},Array{Float64,1},Array{Int64,2},Array{Float64,1},Array{Float64,1},Float64,Float64}) for a proper initialization of this struct.
"""
struct Population_Params
    M::Int64
    G::Int64
    nᵢ::Array{Float64, 1}
    nᵢᵍ::Array{Float64, 2}
    nᵢ_eff::Array{Float64, 1}
    nᵢᵍ_eff::Array{Float64, 2}
    N::Int64
    Nᵍ::Array{Int64, 1}
    pᵍ::Array{Float64, 1}
    edgelist::Array{Int64, 2}
    Rᵢⱼ::Array{Float64, 1}
    mobilityᵍ::Array{Float64, 2}
    kᵍ::Array{Float64, 1}
    C::Array{Float64, 2}
    zᵍ::Array{Float64, 1}
    normᵍ::Array{Float64, 2}
    sᵢ::Array{Float64, 1}
    ξ::Float64
    σ::Float64
end


"""
    Population_Params(M::Int64,
                      G::Int64,
                      nᵢᵍ::Array{Float64, 2},
                      kᵍ::Array{Float64, 1},
                      C::Array{Float64, 2},
                      pᵍ::Array{Float64, 1},
                      edgelist::Array{Int64, 2},
                      Rᵢⱼ::Array{Float64, 1},
                      sᵢ::Array{Float64, 1},
                      ξ::Float64,
                      σ::Float64)

Constructor of the struct [`Population_Params`](@ref).

# Arguments

- `M::Int64`: Number of patches.
- `G::Int64`: Number of population strata.
- `nᵢᵍ::Array{Float64, 2}`: Matrix of size ``G \\times M`` with the population at each strata and patch.
- `kᵍ::Array{Float64, 1}`: Vector of size ``G`` with the average number of contacts of each strata.
- `C::Array{Float64, 2}`: Matrix of size ``G \\times G`` with the probability of contact between different stratifications.
- `pᵍ::Array{Float64, 1}`: Vector of size ``G`` with the degree of mobility of each strata ranged between 0 and 1.
- `edgelist::Array{Int64, 2}`: Matrix of size ``L \\times 2`` containing the directed edgelist between patches, where ``L`` is the number of edges. The IDs of the patches have to go from ``1`` to ``M``.
- `Rᵢⱼ::Array{Float64, 1}`: Vector of size ``L`` containing the transition probabilities for each edge in the edgelist.
- `sᵢ::Array{Float64, 1}`: Vector of size `M` with the surface of each patch.
- `ξ::Float64`: Density factor.
- `σ::Float64`: Average household size.
"""
function Population_Params(M::Int64,
                           G::Int64,
                           nᵢᵍ::Array{Float64, 2},
                           kᵍ::Array{Float64, 1},
                           C::Array{Float64, 2},
                           pᵍ::Array{Float64, 1},
                           edgelist::Array{Int64, 2},
                           Rᵢⱼ::Array{Float64, 1},
                           sᵢ::Array{Float64, 1},
                           ξ::Float64,
                           σ::Float64)

    edgelist, Rᵢⱼ = correct_self_loops(edgelist, Rᵢⱼ, M)

    # Aggregate population count
    Nᵍ = Int.(sum(nᵢᵍ, dims = 2)[:, 1])
    N = sum(Nᵍ)
    nᵢ = sum(nᵢᵍ, dims = 1)[1, :]
    mobilityᵍ = zeros(Float64, G, length(Rᵢⱼ))

    # Init. effective population
    nᵢ_eff = zeros(Float64, M)
    nᵢᵍ_eff = (1 .- pᵍ) .* nᵢᵍ

    # Init. normalization vector
    zᵍ = zeros(Float64, G)
    normᵍ = zeros(Float64, G, M)

    # Compute effective population
    compute_effective_population!(nᵢᵍ_eff, nᵢ_eff, nᵢᵍ, Nᵍ, mobilityᵍ, kᵍ, zᵍ, normᵍ, ξ, pᵍ, sᵢ, edgelist, Rᵢⱼ, M, G)

    return Population_Params(M, G, nᵢ, copy(nᵢᵍ), nᵢ_eff, nᵢᵍ_eff, N, Nᵍ, copy(pᵍ), edgelist, Rᵢⱼ,
                             mobilityᵍ, copy(kᵍ), copy(C), zᵍ, normᵍ, copy(sᵢ), ξ, σ)
end


"""
    update_population_params!(population::Population_Params)

Update population parameters, computing the effective populations and the normalization parameter z if p and k are modified.

# Arguments
- `population::Population_Params`: Structure that contains all the parameters related with the population.
"""
function update_population_params!(population::Population_Params)

    # Reset effective population
    population.nᵢᵍ_eff[:,:] .= (1 .- population.pᵍ) .* population.nᵢᵍ
    population.nᵢ_eff[:] .= zeros(Float64, population.M)
    population.mobilityᵍ[:, :] .= 0.

    # Init. normalization vector
    population.zᵍ .= 0.
    population.normᵍ[:, :] .= 0.

    # Compute effective population
    compute_effective_population!(population.nᵢᵍ_eff,
                                  population.nᵢ_eff,
                                  population.nᵢᵍ,
                                  population.Nᵍ,
                                  population.mobilityᵍ,
                                  population.kᵍ,
                                  population.zᵍ,
                                  population.normᵍ,
                                  population.ξ,
                                  population.pᵍ,
                                  population.sᵢ,
                                  population.edgelist,
                                  population.Rᵢⱼ,
                                  population.M,
                                  population.G)
end


"""
    compute_effective_population!(nᵢᵍ_eff::Array{Float64, 2},
                                  nᵢ_eff::Array{Float64, 1},
                                  nᵢᵍ::Array{Float64, 2},
                                  Nᵍ::Array{Int64, 1},
                                  mobilityᵍ::Array{Float64, 2},
                                  kᵍ::Array{Float64, 1},
                                  zᵍ::Array{Float64, 1},
                                  normᵍ::Array{Float64, 2},
                                  ξ::Float64,
                                  pᵍ::Array{Float64, 1},
                                  sᵢ::Array{Float64, 1},
                                  edgelist::Array{Int64, 2},
                                  Rᵢⱼ::Array{Float64, 1},
                                  M::Int64,
                                  G::Int64)

Compute the effective population at each patch.
"""
function compute_effective_population!(nᵢᵍ_eff::Array{Float64, 2},
                                       nᵢ_eff::Array{Float64, 1},
                                       nᵢᵍ::Array{Float64, 2},
                                       Nᵍ::Array{Int64, 1},
                                       mobilityᵍ::Array{Float64, 2},
                                       kᵍ::Array{Float64, 1},
                                       zᵍ::Array{Float64, 1},
                                       normᵍ::Array{Float64, 2},
                                       ξ::Float64,
                                       pᵍ::Array{Float64, 1},
                                       sᵢ::Array{Float64, 1},
                                       edgelist::Array{Int64, 2},
                                       Rᵢⱼ::Array{Float64, 1},
                                       M::Int64,
                                       G::Int64)

    # Compute the effective population for each strata and patch
    for indx_e in 1:size(edgelist)[1]
        i = edgelist[indx_e, 1]
        j = edgelist[indx_e, 2]

        for g in 1:G
            nᵢᵍ_eff[g, j] += pᵍ[g] * Rᵢⱼ[indx_e] * nᵢᵍ[g, i]
        end
    end

    # Compute the aggerated effective populatoin
    for i in 1:M
        for g in 1:G
            nᵢ_eff[i] += nᵢᵍ_eff[g, i]
        end
    end

    # Correction for populations without effective population
    for i in 1:M
        if nᵢ_eff[i] == 0.
            nᵢ_eff[i] = 1e-7
        end
        for g in 1:G
            if nᵢᵍ_eff[g, i] == 0.
                nᵢᵍ_eff[g, i] = 1e-7
            end
        end
    end

    # Compute the normalization vector
    for g in 1:G
        zᵍ[g] = kᵍ[g] * Nᵍ[g]/sum((2 .- exp.(-ξ .* nᵢ_eff ./ sᵢ)) .* nᵢᵍ_eff[g, :]);
    end

    # Update the precomuted matrices
    for indx_e in 1:length(Rᵢⱼ)
        i = edgelist[indx_e, 1]
        j = edgelist[indx_e, 2]
        for g in 1:G
            mobilityᵍ[g, indx_e] = nᵢᵍ[g, i] * ((1 - pᵍ[g]) * (i == j ? 1. : 0.) + pᵍ[g] * Rᵢⱼ[indx_e])
        end
    end

    for i in 1:M
        for g in 1:G
            normᵍ[g, i] = zᵍ[g] * (2 - exp(-ξ * nᵢ_eff[i] / sᵢ[i]))
        end
    end
end


"""
    set_initial_infected!(epi_params::Epidemic_Params,
                          population::Population_Params,
                          A₀::Array{Float64, 2},
                          I₀::Array{Float64, 2})

Set the initial number of infected individuals on a population. They can be splitted as asymptomatic (A) or Symptomatic (I) individuals.

# Arguments
- `epi_params::Epidemic_Params`: Structure that contains all epidemic parameters and the epidemic spreading information.
- `population::Population_Params`: Structure that contains all the parameters related with the population.
- `A₀::Array{Float64, 2}`: Matrix of size ``G \\times M`` containing the number of asymptomatic infected individuals of each strata on each patch.
- `I₀::Array{Float64, 2}`: Matrix of size ``G \\times M`` containing the number of symptomatic infected individualsof each strata on each patch.

"""
function set_initial_infected!(epi_params::Epidemic_Params,
                               population::Population_Params,
                               A₀::Array{Float64, 2},
                               I₀::Array{Float64, 2})

    # Initial asymptomatic population
    @. epi_params.ρᴬᵍ = A₀ / population.nᵢᵍ

    # Initial sypmtomatic population
    @. epi_params.ρᴵᵍ = I₀ / population.nᵢᵍ

    # Control over division by zero
    epi_params.ρᴬᵍ[isnan.(epi_params.ρᴬᵍ)] .= 0
    epi_params.ρᴵᵍ[isnan.(epi_params.ρᴵᵍ)] .= 0

    # Update the fraction of suceptible individual
    @. epi_params.ρˢᵍ = 1 - epi_params.ρᴵᵍ - epi_params.ρᴬᵍ

end


#------------------------------------------------------
# COMPUTE R FUNCTIONS
#------------------------------------------------------

"""
    compute_R_eff(epi_params::Epidemic_Params,
                  population::Population_Params,
                  T::Int64,
                  M::Int64,
                  G::Int64;
                  τ::Int64 = 20)

Compute the effective reproduction number R.

# Arguments

- `epi_params::Epidemic_Params`: Structure that contains all epidemic parameters and the epidemic spreading information.
- `population::Population_Params`: Structure that contains all the parameters related with the population.

# Optional

- `τ::Int64 = 20`: kernel length.

# Return

Returns a Data Frame containing information about the evolution of the effective reproduction number R.
"""
function compute_R_eff(epi_params::Epidemic_Params,
                       population::Population_Params,
                       τ::Int64 = 20)

    M = population.M
    G = population.G
    T = epi_params.T

    αᵍ = epi_params.αᵍ
    μᵍ = epi_params.μᵍ
    νᵍ = epi_params.νᵍ

    # Setup the kernel
    tʷ = 0:τ
    wᵍ = ones(Float64, G, 1 + τ)
    for g in 1:G
        @. wᵍ[g, :] = epi_params.βᴬ * (1 - αᵍ[g]) ^ tʷ +
                      epi_params.βᴵ * νᵍ[g] * (αᵍ[g] / (μᵍ[g] - αᵍ[g]) * ((1 - αᵍ[g]) ^ tʷ - (1 - μᵍ[g]) ^ tʷ))
    end

    # Init. results dataframe
    df = DataFrame()
    df.strata = repeat(1:G, outer = (T - τ) * M)
    df.src_indx = repeat(1:M, inner = G, outer = (T - τ))
    df.time = repeat(1:(T - τ), inner = G * M)

    # Compute R
    R = zeros(Float64, G, M, (T - τ))
    for t in 1:(T - τ - 1)
        @simd for s in tʷ
            @. R[:, :, t] += epi_params.Qᵢᵍ[:, :, t + s + 1] * wᵍ[:, s + 1]
        end
    end

    df.R = reshape(R, G * M * (T - τ))

    return df
end


#------------------------------------------------------
# OUTPUT FUNCTIONS
#------------------------------------------------------

"""
    store_compartment(epi_params::Epidemic_Params,
                      population::Population_Params,
                      compartment::Char,
                      sufix::String,
                      folder::String)

Store the evolution of the given epidemic compartment for each patch, strata and timestep.

# Arguments

- `epi_params::Epidemic_Params`: Structure that contains all epidemic parameters and the epidemic spreading information.
- `population::Population_Params`: Structure that contains all the parameters related with the population.
- `compartment::Char`: A character in the set {S, E, A, I, H, D, R}.
- `suffix::String`: String used to identify the experiment.
- `folder::String`: String containing the path where the results will be stored.
"""
function store_compartment(epi_params::Epidemic_Params,
                           population::Population_Params,
                           compartment::Char,
                           suffix::String,
                           folder::String)

    M = population.M
    G = population.G
    T = epi_params.T

    # Init. dataframe
    df = DataFrame()
    df.strata = repeat(1:G, outer = T * M)
    df.src_indx = repeat(1:M, inner = G, outer = T)
    df.time = repeat(1:T, inner = G * M)

    # Store number of cases
    if compartment == 'S'
        df.cases = reshape(epi_params.ρˢᵍ .* population.nᵢᵍ, G * M * T)
    elseif compartment == 'E'
        df.cases = reshape(epi_params.ρᴱᵍ .* population.nᵢᵍ, G * M * T)
    elseif compartment == 'A'
        df.cases = reshape(epi_params.ρᴬᵍ .* population.nᵢᵍ, G * M * T)
    elseif compartment == 'I'
        df.cases = reshape(epi_params.ρᴵᵍ .* population.nᵢᵍ, G * M * T)
    elseif compartment == 'H'
        df.cases = reshape(epi_params.ρᴴᵍ .* population.nᵢᵍ, G * M * T)
    elseif compartment == 'D'
        df.cases = reshape(epi_params.ρᴰᵍ .* population.nᵢᵍ, G * M * T)
    elseif compartment == 'R'
        df.cases = reshape(epi_params.ρᴿᵍ .* population.nᵢᵍ, G * M * T)
    end

    CSV.write(@sprintf("%s/output_%s_%s.csv", folder, compartment, suffix), df)
end


"""
    store_R_eff(epi_params::Epidemic_Params,
                population::Population_Params,
                suffix::String,
                folder::String;
                τ::Int64 = 20)

Compute and store the effective reproduction number R

# Arguments

- `epi_params::Epidemic_Params`: Structure that contains all epidemic parameters and the epidemic spreading information.
- `population::Population_Params`: Structure that contains all the parameters related with the population.
- `suffix::String`: String used to identify the experiment.
- `folder::String`: String containing the path where the results will be stored.

# Optional
- `τ::Int64 = 20`: kernel length.
"""
function store_R_eff(epi_params::Epidemic_Params,
                     population::Population_Params,
                     suffix::String,
                     folder::String;
                     τ::Int64 = 20)

    df = compute_R_eff(epi_params, population, τ)

    # Write the results
    CSV.write(@sprintf("%s/output_Reff_%s.csv", folder, suffix), df)
end


#------------------------------------------------------
# EXTRA FUNCTIONS
#------------------------------------------------------

"""
    correct_self_loops(edgelist::Array{Int64, 2},
                       Rᵢⱼ::Array{Float64, 1},
                       M::Int64)

Repair weighted matrices to ensure that all nodes are represented adding an additional self-loop.

# Arguments
- `edgelist::Array{Int64, 2}`: Matrix containing the directed edgelist between patches. The IDs of the patches have to go from 1 to M.
- `Rᵢⱼ::Array{Float64, 1}`: Vector containing the transition probabilities for each edge in the edgelist.
- `M::Int64`: Number of patches.

# Return

Returns the corrected list of edges and transition probabilities.
"""
function correct_self_loops(edgelist::Array{Int64, 2},
                            Rᵢⱼ::Array{Float64, 1},
                            M::Int64)

    self_loops = falses(M)
    k_in = zeros(Int64, M)
    k_out = zeros(Int64, M)

    for indx_e in 1:size(edgelist)[1]
        i = edgelist[indx_e, 1]
        j = edgelist[indx_e, 2]

        k_in[j] += 1
        k_out[i] += 1

        if i == j
            self_loops[i] = true
        end
    end

    to_add = collect(1:M)[(k_out .== 0) .& .!(self_loops)]

    new_edges = zeros(Int64, length(to_add), 2)
    new_edges[:,1] = to_add
    new_edges[:,2] = to_add

    edgelist = vcat(edgelist, new_edges)
    Rᵢⱼ = vcat(Rᵢⱼ, ones(Float64, length(to_add)))

    to_add = collect(1:M)[(k_out .!= 0) .& .!(self_loops)]

    new_edges = zeros(Int64, length(to_add), 2)
    new_edges[:,1] = to_add
    new_edges[:,2] = to_add

    edgelist = vcat(edgelist, new_edges)
    Rᵢⱼ = vcat(Rᵢⱼ, zeros(Float64, length(to_add)))

    return (edgelist, Rᵢⱼ)
end

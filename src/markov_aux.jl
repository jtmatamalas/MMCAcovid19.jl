#------------------------------------------------------
# EPIDEMIC PARAMS RELATED FUNCTIONS
#------------------------------------------------------

"""
    Epidemic_Params

Struct that contains the parameters related with the epidemic parameters and
compartmental evolution

# Fields

All the parameters contained in this structure are probabilities ranged between
0 and 1.

## Epidemic parameters

- `βᴬ::Float64`: Infectivity of asymptomatic.
- `βᴵ::Float64`: Infectivity of infected.
- `νᵍ::Array{Float64, 1}`: Isolation factor for each strata.
- `ηᵍ::Array{Float64, 1}`: Latent rate for each strata.
- `αᵍ::Array{Float64, 1}`: Asymptomatic infectious rate for each strata.
- `μᵍ::Array{Float64, 1}`: Escape rate for each strata.
- `χᵍ::Array{Float64, 1}`: ICU discharge rate for each strata.
- `ψᵍ::Array{Float64, 1}`: Death rate for each strata.
- `ωᵍ::Array{Float64, 1}`: Fatality rate of ICU patients for each strata.
- `γᵍ::Array{Float64, 1}`: Fraction of cases requiring ICU for each strata.

## Compartmental evolution

- `ρˢᵍ::Array{Float64, 3}`: ``G \times M \times T`` array containing infomation
  about the evolution of fraction of suceptible individuals for each strata and
  patch.
- `ρᴱᵍ::Array{Float64, 3}`: ``G \times M \times T`` array containing infomation
  about the evolution of fraction of exposed individuals for each strata and
  patch.
- `ρᴬᵍ::Array{Float64, 3}`: ``G \times M \times T`` array containing infomation
  about the evolution of fraction of asymptomatic individuals for each strata and
  patch.
- `ρᴴᵍ::Array{Float64, 3}`: ``G \times M \times T`` array containing infomation
  about the evolution of fraction of ICU patients individuals for each strata and
  patch.
- `ρᴵᵍ::Array{Float64, 3}`: ``G \times M \times T`` array containing infomation
  about the evolution of fraction of symptomatic individuals for each strata and
  patch.
- `ρᴰᵍ::Array{Float64, 3}`: ``G \times M \times T`` array containing infomation
  about the evolution of fraction of death individuals for each strata and
  patch.
- `ρᴿᵍ::Array{Float64, 3}`: ``G \times M \times T`` array containing infomation
  about the evolution of fraction of recovered individuals for each strata and
  patch.

## Auxiliar

- `CHᵢᵍ::Array{Float64, 2}`: Fraction of securely confined individuals from each
  strata and patch.
- `Qᵢᵍ::Array{Float64, 3}`: Suceptible contacts available of each strata on a
  given patch.

"""
# Define the structure of parameters
struct Epidemic_Params
    #Epidemic parameters
    βᴬ::Float64
    βᴵ::Float64
    νᵍ::Array{Float64, 1}
    ηᵍ::Array{Float64, 1}
    αᵍ::Array{Float64, 1}
    μᵍ::Array{Float64, 1}
    χᵍ::Array{Float64, 1}
    ψᵍ::Array{Float64, 1}
    ωᵍ::Array{Float64, 1}
    γᵍ::Array{Float64, 1}

    # Compartments evolution
    ρˢᵍ::Array{Float64, 3}
    ρᴱᵍ::Array{Float64, 3}
    ρᴬᵍ::Array{Float64, 3}
    ρᴴᵍ::Array{Float64, 3}
    ρᴵᵍ::Array{Float64, 3}
    ρᴰᵍ::Array{Float64, 3}
    ρᴿᵍ::Array{Float64, 3}
    CHᵢᵍ::Array{Float64, 2}

    # R_t related arrays
    Qᵢᵍ::Array{Float64, 3}

    # Constructor
    """
        function Epidemic_Params( βᴬ::Float64,
                                  βᴵ::Float64,
                                  νᵍ::Array{Float64, 1},
                                  ηᵍ::Array{Float64, 1},
                                  αᵍ::Array{Float64, 1},
                                  μᵍ::Array{Float64, 1},
                                  χᵍ::Array{Float64, 1},
                                  ψᵍ::Array{Float64, 1},
                                  ωᵍ::Array{Float64, 1},
                                  γᵍ::Array{Float64, 1},
                                  M::Int64,
                                  G::Int64,
                                  timesteps::Int64)

    Constructor of the struct Epidemic_Params

    # Arguments

    - `βᴬ::Float64`: Infectivity of asymptomatic.
    - `βᴵ::Float64`: Infectivity of infected.
    - `νᵍ::Array{Float64, 1}`: Isolation factor for each strata.
    - `ηᵍ::Array{Float64, 1}`: Latent rate for each strata.
    - `αᵍ::Array{Float64, 1}`: Asymptomatic infectious rate for each strata.
    - `μᵍ::Array{Float64, 1}`: Escape rate for each strata.
    - `χᵍ::Array{Float64, 1}`: ICU discharge rate for each strata.
    - `ψᵍ::Array{Float64, 1}`: Death rate for each strata.
    - `ωᵍ::Array{Float64, 1}`: Fatality rate of ICU patients for each strata.
    - `γᵍ::Array{Float64, 1}`: Fraction of cases requiring ICU for each strata.
    - `M::Int64`: Number of patches.
    - `G::Int64`: Number of stratifications.
    - `timesteps::Int64`: Number of epidemic timesteps.
    """
    function Epidemic_Params( βᴬ::Float64,
                              βᴵ::Float64,
                              νᵍ::Array{Float64, 1},
                              ηᵍ::Array{Float64, 1},
                              αᵍ::Array{Float64, 1},
                              μᵍ::Array{Float64, 1},
                              χᵍ::Array{Float64, 1},
                              ψᵍ::Array{Float64, 1},
                              ωᵍ::Array{Float64, 1},
                              γᵍ::Array{Float64, 1},
                              M::Int64, # Num. patches
                              G::Int64, # Num. Statifications
                              timesteps::Int64)

        # Allocate memory for simulations
        ρˢᵍ = ones(Float64, G, M, timesteps)
        ρᴱᵍ = zeros(Float64, G, M, timesteps)
        ρᴬᵍ = zeros(Float64, G, M, timesteps)
        ρᴵᵍ = zeros(Float64, G, M, timesteps)
        ρᴴᵍ = zeros(Float64, G, M, timesteps)
        ρᴰᵍ = zeros(Float64, G, M, timesteps)
        ρᴿᵍ = zeros(Float64, G, M, timesteps)
        CHᵢᵍ = zeros(Float64, G, M)
        Qᵢᵍ = zeros(Float64, G, M, timesteps)

        return new(βᴬ, βᴵ, νᵍ, ηᵍ, αᵍ, μᵍ, χᵍ, ψᵍ, ωᵍ, γᵍ,
                   ρˢᵍ, ρᴱᵍ, ρᴬᵍ, ρᴵᵍ, ρᴴᵍ, ρᴰᵍ, ρᴿᵍ, CHᵢᵍ,
                   Qᵢᵍ)
    end
end


"""
    reset_epi_param!(epi_param::Epidemic_Params)

Reset the ρs to reuse the structure and avoid additional allocatiosn.

# Arguments
- `epi_param::Epidemic_Params`: structure to reset
"""

function reset_epi_param!(epi_param::Epidemic_Params)
    epi_param.ρˢᵍ .= 1.
    epi_param.ρᴱᵍ .= 0.
    epi_param.ρᴬᵍ .= 0.
    epi_param.ρᴴᵍ .= 0.
    epi_param.ρᴵᵍ .= 0.
    epi_param.ρᴰᵍ .= 0.
    epi_param.ρᴿᵍ .= 0.
    epi_param.CHᵢᵍ .= 0.

    # Rt structures
    epi_param.Qᵢᵍ .= 0.
end

#------------------------------------------------------
# PATCH AND POPULATION RELATED FUNCTIONS
#------------------------------------------------------

"""
    Population_Params

Struct that contains the parameters related with demographics and mobility

# Fields

- `nᵢ:Array{Float64, 1}`: Population of each patch.
- `nᵢᵍ:Array{Float64, 2}`: Population of each strata on each patch.
- `nᵢ_eff:Array{Float64, 1}`: Effective population on each patch after taking
  into account mobility.
- `nᵢᵍ_eff:Array{Float64, 1}`: Effective population of each strata on each patch
  after taking into account mobility.
- `N::Int64`: Total population.
- `Nᵍ::Array{Int64, 1}`: Total population of each strata.
- `mobilityᵍ::Array{Float64, 1}`: Effective mobility of each strata between
  patches. (Used to optimize the speed).
- `sᵢ::Array{Float64, 1}`: Surface of each patch.
- `kᵍ::Array{Float64, 1}`: average number of contacts of each strata.
- `zᵍ::Array{Float64, 1}`: nomalization factor each strata.
- `normᵍ::Array{Float64, 2}`: normalization of each strata (Used to optimize the speed).
- `σ::Float64`: Average household size.
- `ξ::Float64`: Densisty factor.
"""
struct Population_Params
    nᵢ::Array{Float64, 1}
    nᵢᵍ::Array{Float64, 2}
    nᵢ_eff::Array{Float64, 1}
    nᵢᵍ_eff::Array{Float64, 2}
    N::Int64
    Nᵍ::Array{Int64, 1}
    mobilityᵍ::Array{Float64, 2}
    sᵢ::Array{Float64, 1}
    kᵍ::Array{Float64, 1}
    zᵍ::Array{Float64, 1}
    normᵍ::Array{Float64, 2}
    σ::Float64
    ξ::Float64

    """
        Population_Params(nᵢᵍ::Array{Float64, 2},
                          kᵍ::Array{Float64, 1},
                          pᵍ::Array{Float64, 1},
                          edgelist::Array{Int64, 2},
                          Rᵢⱼ::Array{Float64, 1},
                          M::Int64,
                          G::Int64,
                          ξ::Float64)

    Constructor of the struct Population_Params

    # Arguments

    - `nᵢᵍ::Array{Float64, 2}`: population of each strata on each patch.
    - `kᵍ::Array{Float64, 1}`: average number of contacts of each strata.
    - `pᵍ::Array{Float64, 1}`: Degree of mobility of each strata ranged between 0 and 1.
    - `edgelist::Array{Float64, 2}`: ``M \\times 2`` matrix containing the directed
    edgelist between patches. The IDs of the patches have to go from 1 to M.
    - `Rᵢⱼ::Array{Float64, 1}`: Vector containing the transition probabilities for
    each edge in the edgelist.
    - `M::Int64`: Number of patches.
    - `G::Int64`: Number of stratifications of the population.
    - `ξ::Float64`: Density factor.
    """
    function Population_Params(nᵢᵍ::Array{Float64, 2},
                               kᵍ::Array{Float64, 1},
                               pᵍ::Array{Float64, 1},
                               edgelist::Array{Int64, 2},
                               Rᵢⱼ::Array{Float64, 1},
                               M::Int64,
                               G::Int64,
                               ξ::Float64)

        # Aggregate population count
        Nᵍ = Int.(sum(nᵢᵍ, dims = 2)[:, 1])
        N = sum(Nᵍ)
        nᵢ = sum(nᵢᵍ, dims = 1)[1,:]
        mobilityᵍ = zeros(Float64, G, length(Rᵢⱼ))

        # Init. effective population
        nᵢ_eff = zeros(Float64, M)
        nᵢᵍ_eff = (1 .- pᵍ) .* nᵢᵍ

        # Init. normalization vector
        zᵍ = zeros(Float64, G)
        normᵍ = zeros(Float64, G, M)

        # Compute effective population
        compute_effective_population!(nᵢᵍ_eff, nᵢ_eff, nᵢᵍ, Nᵍ, mobilityᵍ, kᵍ, zᵍ, normᵍ, ξ, pᵍ, sᵢ, edgelist, Rᵢⱼ, M, G)

        return new(nᵢ, nᵢᵍ, nᵢ_eff, nᵢᵍ_eff, N, Nᵍ, mobilityᵍ, sᵢ, copy(kᵍ), zᵍ, normᵍ, σ, ξ)
    end
end


"""
    update_population_param!(population::Population_Params,
                             pᵍ::Array{Float64, 1},
                             edgelist::Array{Int64, 2},
                             Rᵢⱼ::Array{Float64, 1},
                             M::Int64,
                             G::Int64)

Update population parameters computing the effective populations and the
normalization parameter z if p are k are modified

# Arguments
- `population::Population_Params`: Structure that contains all the parameters
  related with the population
- `pᵍ::Array{Float64, 1}`: Degree of mobility of each strata ranged between 0 and 1.
- `edgelist::Array{Float64, 2}`: ``M \\times 2`` matrix containing the directed
  edgelist between patches. The IDs of the patches have to go from 1 to M
- `Rᵢⱼ::Array{Float64, 1}`: Vector containing the transition probabilities for
  each edge in the edgelist.
- `M::Int64`: Number of patches.
- `G::Int64`: Number of stratifications of the population.
"""
function update_population_param!(population::Population_Params,
                                  pᵍ::Array{Float64, 1},
                                  edgelist::Array{Int64, 2},
                                  Rᵢⱼ::Array{Float64, 1},
                                  M::Int64,
                                  G::Int64)

    # Reset effective population
    population.nᵢᵍ_eff[:,:] .= (1 .- pᵍ) .* population.nᵢᵍ
    population.nᵢ_eff[:] .= zeros(Float64, M)
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
                                  pᵍ,
                                  population.sᵢ,
                                  edgelist,
                                  Rᵢⱼ,
                                  M,
                                  G)
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

    Compute the effective population of each patch.
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
    setInitialInfected!(epi_param::Epidemic_Params,
                        population_struct::Population_Params,
                        A₀::Array{Float64, 2},
                        I₀::Array{Float64, 2})

Set the initial number of infected individuals on a population. They can be
splitted as asymptomatic (A) or Symptomatic (I) individuals.

# Arguments
- `epi_params::Epidemic_Params`: Structure that contains all realted epidemic
  parameters.
- `population::Population_Params`: Structure that contains all the parameters
  related with the population
- `A₀::Array{Float64, 2}`: ``G \times M`` matr containing the number of
  asymptomatic infected individuals of each strata on each patch.
- `I₀::Array{Float64, 2}`: ``G \times M`` matr containing the number of
  symptomatic infected individualsof each strata on each patch.

"""
function setInitialInfected!(epi_param::Epidemic_Params,
                             population_struct::Population_Params,
                             A₀::Array{Float64, 2},
                             I₀::Array{Float64, 2})

    # Init. asymptomatic population
    @. epi_param.ρᴬᵍ = A₀ / population_struct.nᵢᵍ

    # Init. sypmtomatic population
    @. epi_param.ρᴵᵍ = I₀ / population_struct.nᵢᵍ

    # Control over division by zero
    epi_param.ρᴬᵍ[isnan.(epi_param.ρᴬᵍ)] .= 0
    epi_param.ρᴵᵍ[isnan.(epi_param.ρᴵᵍ)] .= 0

    # Update the fraction of suceptible individual
    @. epi_param.ρˢᵍ = 1 - epi_param.ρᴵᵍ - epi_param.ρᴬᵍ

end

#------------------------------------------------------
# OUTPUT FUNCTIONS
#------------------------------------------------------

"""
    store_compartment(epi_param::Epidemic_Params,
                      population::Population_Params,
                      T::Int64,
                      M::Int64,
                      G::Int64,
                      compartment::Char,
                      sufix::String,
                      folder::String)

Store the number

# Arguments

- `epi_params::Epidemic_Params`: Structure that contains all realted epidemic
  parameters.
- `population::Population_Params`: Structure that contains all the parameters
  related with the population
- `T::Int64`: Number of timesteps
- `M::Int64`: Number of patches
- `G::Int64`: Number of statifications
- `compartment::Char`: A char in the set {S, E, A, I, H, D, R}
- `suffix::String`: String used to identify the experiment
- `folder::String`: String containing the path where the results will be stored.
"""
function store_compartment(epi_param::Epidemic_Params,
                           population::Population_Params,
                           T::Int64,
                           M::Int64,
                           G::Int64,
                           compartment::Char,
                           sufix::String,
                           folder::String)

    # Init. dataframe
    df = DataFrame()
    df.strata = repeat(1:G, outer = T * M)
    df.src_indx = repeat(1:M, inner = G, outer = T)
    df.time = repeat(1:T, inner = G * M)

    # Store number of cases
    if compartment == 'S'
        df.cases = reshape(epi_param.ρˢᵍ .* population.nᵢᵍ, G * M * T)
    elseif compartment == 'E'
        df.cases = reshape(epi_param.ρᴱᵍ .* population.nᵢᵍ, G * M * T)
    elseif compartment == 'A'
        df.cases = reshape(epi_param.ρᴬᵍ .* population.nᵢᵍ, G * M * T)
    elseif compartment == 'I'
        df.cases = reshape(epi_param.ρᴵᵍ .* population.nᵢᵍ, G * M * T)
    elseif compartment == 'H'
        df.cases = reshape(epi_param.ρᴴᵍ .* population.nᵢᵍ, G * M * T)
    elseif compartment == 'D'
        df.cases = reshape(epi_param.ρᴰᵍ .* population.nᵢᵍ, G * M * T)
    elseif compartment == 'R'
        df.cases = reshape(epi_param.ρᴿᵍ .* population.nᵢᵍ, G * M * T)
    end

    CSV.write(@sprintf("%s/output_%s_%s.csv", folder, compartment, sufix), df)
end

"""
    store_R0(epi_param::Epidemic_Params,
                    population::Population_Params,
                    T::Int64,
                    M::Int64,
                    G::Int64,
                    suffix::String,
                    folder::String;
                    τ::Int64 = 20)

Compute and store the R

# Arguments

- `epi_params::Epidemic_Params`: Structure that contains all realted epidemic
  parameters.
- `population::Population_Params`: Structure that contains all the parameters
  related with the population
- `T::Int64`: Number of timesteps
- `M::Int64`: Number of patches
- `G::Int64`: Number of statifications
- `suffix::String`: String used to identify the experiment
- `folder::String`: String containing the path where the results will be stored.

## Optional
- `τ::Int64 = 20`: kernel length
"""
function store_R0(epi_param::Epidemic_Params,
                  population::Population_Params,
                  T::Int64,
                  M::Int64,
                  G::Int64,
                  sufix::String,
                  folder::String,
                  τ = 20)

    # Setup the kernel
    tʷ = 0:τ
    wᵍ = ones(Float64, G, 1 + τ)
    for g in 1:G
        @. wᵍ[g, :] = epi_param.βᴬ * (1 - αᵍ[g]) ^ tʷ +
                      epi_param.βᴵ * νᵍ[g] * (αᵍ[g] / (μᵍ[g] - αᵍ[g]) * ((1 - αᵍ[g]) ^ tʷ - (1 - μᵍ[g]) ^ tʷ))
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
            @. R[:, :, t] += epi_param.Qᵢᵍ[:, :, t + s + 1] * wᵍ[:, s + 1]
        end
    end

    df.R = reshape(R, G * M * (T - τ))

    # Write the results
    CSV.write(@sprintf("%s/output_R_%s.csv", folder, suffix), df)
end

#------------------------------------------------------
# EXTRA FUNCTIONS
#------------------------------------------------------

"""
    correct_self_loops(edgelist::Array{Int64, 2},
                       Rᵢⱼ::Array{Float64, 1},
                       M::Int64)

Repair weigted matrices to ensure that all nodes are represented
adding an additional self-loop

# Arguments
- `edgelist::Array{Float64, 2}`: ``M \\times 2`` matrix containing the directed
  edgelist between patches. The IDs of the patches have to go from 1 to M
- `Rᵢⱼ::Array{Float64, 1}`: Vector containing the transition probabilities for
  each edge in the edgelist.
- `M::Int64`: Number of patches.
"""
function correct_self_loops(edgelist::Array{Int64, 2},
                             weights::Array{Float64, 1},
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
    weights = vcat(Rᵢⱼ, ones(Float64, length(to_add)))

    to_add = collect(1:M)[(k_out .!= 0) .& .!(self_loops)]

    new_edges = zeros(Int64, length(to_add), 2)
    new_edges[:,1] = to_add
    new_edges[:,2] = to_add

    edgelist = vcat(edgelist, new_edges)
    weights = vcat(Rᵢⱼ, zeros(Float64, length(to_add)))

    return (edgelist, Rᵢⱼ)
end

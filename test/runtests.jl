using Test
using DataFrames
using MMCAcovid19

## -----------------------------------------------------------------------------
## Population parameters
## -----------------------------------------------------------------------------

G = 3
M = 5

nᵢᵍ = [ 4995.0   9875.0  14970.0   30010.0   40326.0
       30107.0  59630.0  90009.0  179745.0  239983.0
       15145.0  29827.0  45086.0   90266.0  120026.0]

C = [0.5980 0.3849 0.0171
     0.2440 0.7210 0.0350
     0.1919 0.5705 0.2376]

edgelist = [1  1; 1  2; 1  3; 1  5; 2  1; 2  2; 2  3; 2  4;
            2  5; 3  1; 3  2; 3  3; 3  4; 3  5; 4  1; 4  3;
            4  4; 4  5; 5  1; 5  2; 5  3; 5  5]

Rᵢⱼ = [0.3288; 0.0905; 0.0995; 0.4812; 0.3916; 0.2213; 0.1052; 0.2775;
       0.0044; 0.0233; 0.5205; 0.0117; 0.0807; 0.3638; 0.5156; 0.0579;
       0.0218; 0.4047; 0.3081; 0.2862; 0.0621; 0.3436]

kᵍ = [11.8, 13.3, 6.6]
kᵍ_h = [3.15, 3.17, 3.28]
kᵍ_w = [1.72, 5.18, 0.0]
pᵍ = [0.0, 1.0, 0.05]
sᵢ = [10.6, 23.0, 26.6, 5.7, 61.6]
ξ = 0.01
σ = 2.5

population = Population_Params(G, M, nᵢᵍ, kᵍ, kᵍ_h, kᵍ_w, C, pᵍ, edgelist, Rᵢⱼ, sᵢ, ξ, σ)

## Test population parameters
@testset "Population_Params" begin

    @test population.G == G
    @test population.M == M
    @test population.nᵢᵍ == nᵢᵍ
    @test population.kᵍ == kᵍ
    @test population.kᵍ_h == kᵍ_h
    @test population.kᵍ_w == kᵍ_w
    @test population.C == C
    @test population.pᵍ == pᵍ
    @test population.edgelist == edgelist
    @test population.Rᵢⱼ == Rᵢⱼ
    @test population.sᵢ == sᵢ
    @test population.ξ == ξ
    @test population.σ == σ

    @test length(population.nᵢ) == M
    @test length(population.Nᵍ) == G
    @test length(population.nᵢ_eff) == M
    @test length(population.zᵍ) == G
    @test size(population.nᵢᵍ_eff) == (G, M)
    @test size(population.normᵍ) == (G, M)
    @test size(population.mobilityᵍ) == (G, length(Rᵢⱼ))

    @test population.N == Int64(sum(nᵢᵍ))
    @test population.nᵢ ≈ sum(nᵢᵍ, dims=1)[1, :]  atol=0.0001
    @test population.Nᵍ ≈ sum(nᵢᵍ, dims=2)[:, 1]  atol=0.0001

    @test sum(population.nᵢᵍ_eff) ≈ sum(nᵢᵍ)  atol=0.0001
    @test population.nᵢ_eff ≈ sum(population.nᵢᵍ_eff, dims=1)[1, :]  atol=0.0001

end


## -----------------------------------------------------------------------------
## Epidemic parameters
## -----------------------------------------------------------------------------

βᴵ = 0.075
βᴬ = 0.5 * βᴵ
ηᵍ = [1/2.444, 1/2.444, 1/2.444]
αᵍ = [1/5.671, 1/2.756, 1/2.756]
μᵍ = [1/1.0, 1/3.915, 1/3.915]
θᵍ = [0.0, 0.008, 0.047]
γᵍ = [0.0003, 0.003, 0.026]
ζᵍ = [1/7.084, 1/7.084, 1/7.084]
λᵍ = [1/4.084, 1/4.084, 1/4.084]
ωᵍ = [0.3, 0.3, 0.3]
ψᵍ = [1/7.0, 1/7.0, 1/7.0]
χᵍ = [1/20.0, 1/20.0, 1/20.0]

T = 200

epi_params = Epidemic_Params(βᴵ, βᴬ, ηᵍ, αᵍ, μᵍ, θᵍ, γᵍ,
                             ζᵍ, λᵍ, ωᵍ, ψᵍ, χᵍ, G, M, T)

## Test population parameters
@testset "Epidemic_Params" begin

    @test epi_params.βᴵ[1] == βᴵ
    @test epi_params.βᴬ[1] == βᴬ
    @test epi_params.ηᵍ == ηᵍ
    @test epi_params.αᵍ == αᵍ
    @test epi_params.μᵍ == μᵍ
    @test epi_params.θᵍ == θᵍ
    @test epi_params.γᵍ == γᵍ
    @test epi_params.ζᵍ == ζᵍ
    @test epi_params.λᵍ == λᵍ
    @test epi_params.ωᵍ == ωᵍ
    @test epi_params.ψᵍ == ψᵍ
    @test epi_params.χᵍ == χᵍ
    @test epi_params.T == T

    @test size(epi_params.ρˢᵍ) == (G, M, T)
    @test size(epi_params.ρᴱᵍ) == (G, M, T)
    @test size(epi_params.ρᴬᵍ) == (G, M, T)
    @test size(epi_params.ρᴵᵍ) == (G, M, T)
    @test size(epi_params.ρᴾᴴᵍ) == (G, M, T)
    @test size(epi_params.ρᴾᴰᵍ) == (G, M, T)
    @test size(epi_params.ρᴴᴿᵍ) == (G, M, T)
    @test size(epi_params.ρᴴᴰᵍ) == (G, M, T)
    @test size(epi_params.ρᴰᵍ) == (G, M, T)
    @test size(epi_params.ρᴿᵍ) == (G, M, T)
    @test size(epi_params.CHᵢᵍ) == (G, M)
    @test size(epi_params.Qᵢᵍ) == (G, M, T)

end


## -----------------------------------------------------------------------------
## Initialization of the epidemics
## -----------------------------------------------------------------------------

E₀ = zeros(G, M)

A₀ = zeros(G, M)
A₀[2, 5] = 2.0
A₀[3, 3] = 1.0

I₀ = zeros(G, M)
I₀[2, 5] = 1.0

set_initial_infected!(epi_params, population, E₀, A₀, I₀)


## Test Initialization of the epidemics
@testset "Initialization_Epidemics" begin

    @test sum(population.nᵢᵍ .* epi_params.ρˢᵍ[:, :, 1]) ≈
        population.N - (sum(E₀) + sum(A₀) + sum(I₀))  atol=0.0001
    @test sum(population.nᵢᵍ .* epi_params.ρᴱᵍ[:, :, 1]) ≈ sum(E₀)  atol=0.0001
    @test sum(population.nᵢᵍ .* epi_params.ρᴬᵍ[:, :, 1]) ≈ sum(A₀)  atol=0.0001
    @test sum(population.nᵢᵍ .* epi_params.ρᴵᵍ[:, :, 1]) ≈ sum(I₀)  atol=0.0001
    @test sum(population.nᵢᵍ .* epi_params.ρᴾᴴᵍ[:, :, 1]) ≈ 0.0  atol=0.0001
    @test sum(population.nᵢᵍ .* epi_params.ρᴾᴰᵍ[:, :, 1]) ≈ 0.0  atol=0.0001
    @test sum(population.nᵢᵍ .* epi_params.ρᴴᴿᵍ[:, :, 1]) ≈ 0.0  atol=0.0001
    @test sum(population.nᵢᵍ .* epi_params.ρᴴᴰᵍ[:, :, 1]) ≈ 0.0  atol=0.0001
    @test sum(population.nᵢᵍ .* epi_params.ρᴰᵍ[:, :, 1]) ≈ 0.0  atol=0.0001
    @test sum(population.nᵢᵍ .* epi_params.ρᴿᵍ[:, :, 1]) ≈ 0.0  atol=0.0001
    @test sum(population.nᵢᵍ .* epi_params.CHᵢᵍ[:, :]) ≈ 0.0  atol=0.0001

end


## -----------------------------------------------------------------------------
## Epidemic spreading without containment
## -----------------------------------------------------------------------------

run_epidemic_spreading_mmca!(epi_params, population, verbose = false)


# Test epidemic spreading without containment
@testset "Epidemic_Spreading_No_Contaiment" begin

    @test sum(population.nᵢᵍ .* epi_params.ρˢᵍ[:, :, T]) ≈ 30673.5333  atol=0.0001
    @test sum(population.nᵢᵍ .* epi_params.ρᴱᵍ[:, :, T]) ≈ 0.0  atol=0.0001
    @test sum(population.nᵢᵍ .* epi_params.ρᴬᵍ[:, :, T]) ≈ 0.0  atol=0.0001
    @test sum(population.nᵢᵍ .* epi_params.ρᴵᵍ[:, :, T]) ≈ 0.0  atol=0.0001
    @test sum(population.nᵢᵍ .* epi_params.ρᴾᴴᵍ[:, :, T]) ≈ 0.0  atol=0.0001
    @test sum(population.nᵢᵍ .* epi_params.ρᴾᴰᵍ[:, :, T]) ≈ 0.0  atol=0.0001
    @test sum(population.nᵢᵍ .* epi_params.ρᴴᴿᵍ[:, :, T]) ≈ 7.2739  atol=0.0001
    @test sum(population.nᵢᵍ .* epi_params.ρᴴᴰᵍ[:, :, T]) ≈ 0.0  atol=0.0001
    @test sum(population.nᵢᵍ .* epi_params.ρᴰᵍ[:, :, T]) ≈ 20301.7580  atol=0.0001
    @test sum(population.nᵢᵍ .* epi_params.ρᴿᵍ[:, :, T]) ≈ 949017.4345  atol=0.0001
    @test sum(population.nᵢᵍ .* epi_params.CHᵢᵍ[:, :]) ≈ 0.0  atol=0.0001

end


## -----------------------------------------------------------------------------
## Reset the parameters
## -----------------------------------------------------------------------------

reset_params!(epi_params, population)


## Test Reset of the epidemics
@testset "Reset_Parameters" begin

  @test sum(population.nᵢᵍ .* epi_params.ρˢᵍ[:, :, 1]) ≈ population.N  atol=0.0001
  @test sum(population.nᵢᵍ .* epi_params.ρᴱᵍ[:, :, 1]) ≈ 0.0  atol=0.0001
  @test sum(population.nᵢᵍ .* epi_params.ρᴬᵍ[:, :, 1]) ≈ 0.0  atol=0.0001
  @test sum(population.nᵢᵍ .* epi_params.ρᴵᵍ[:, :, 1]) ≈ 0.0  atol=0.0001
  @test sum(population.nᵢᵍ .* epi_params.ρᴾᴴᵍ[:, :, 1]) ≈ 0.0  atol=0.0001
  @test sum(population.nᵢᵍ .* epi_params.ρᴾᴰᵍ[:, :, 1]) ≈ 0.0  atol=0.0001
  @test sum(population.nᵢᵍ .* epi_params.ρᴴᴿᵍ[:, :, 1]) ≈ 0.0  atol=0.0001
  @test sum(population.nᵢᵍ .* epi_params.ρᴴᴰᵍ[:, :, 1]) ≈ 0.0  atol=0.0001
  @test sum(population.nᵢᵍ .* epi_params.ρᴰᵍ[:, :, 1]) ≈ 0.0  atol=0.0001
  @test sum(population.nᵢᵍ .* epi_params.ρᴿᵍ[:, :, 1]) ≈ 0.0  atol=0.0001
  @test sum(population.nᵢᵍ .* epi_params.CHᵢᵍ[:, :]) ≈ 0.0  atol=0.0001

end


## -----------------------------------------------------------------------------
## Epidemic spreading with single containment
## -----------------------------------------------------------------------------

tᶜ = 30
κ₀ = 0.65
ϕ = 0.174
δ = 0.207

set_initial_infected!(epi_params, population, E₀, A₀, I₀)
run_epidemic_spreading_mmca!(epi_params, population; tᶜ = tᶜ, κ₀ = κ₀, ϕ = ϕ,
                             δ = δ, verbose = false)


## Test epidemic spreading with single containment
@testset "Epidemic_Spreading_With_Single_Contaiment" begin

    @test sum(population.nᵢᵍ .* epi_params.ρˢᵍ[:, :, T]) ≈ 455606.5390  atol=0.0001
    @test sum(population.nᵢᵍ .* epi_params.ρᴱᵍ[:, :, T]) ≈ 1.0411  atol=0.0001
    @test sum(population.nᵢᵍ .* epi_params.ρᴬᵍ[:, :, T]) ≈ 1.4510  atol=0.0001
    @test sum(population.nᵢᵍ .* epi_params.ρᴵᵍ[:, :, T]) ≈ 2.1264  atol=0.0001
    @test sum(population.nᵢᵍ .* epi_params.ρᴾᴴᵍ[:, :, T]) ≈ 0.0232  atol=0.0001
    @test sum(population.nᵢᵍ .* epi_params.ρᴾᴰᵍ[:, :, T]) ≈ 0.0985  atol=0.0001
    @test sum(population.nᵢᵍ .* epi_params.ρᴴᴿᵍ[:, :, T]) ≈ 0.3454  atol=0.0001
    @test sum(population.nᵢᵍ .* epi_params.ρᴴᴰᵍ[:, :, T]) ≈ 0.0167  atol=0.0001
    @test sum(population.nᵢᵍ .* epi_params.ρᴰᵍ[:, :, T]) ≈ 293.3980  atol=0.0001
    @test sum(population.nᵢᵍ .* epi_params.ρᴿᵍ[:, :, T]) ≈ 15836.2800 atol=0.0001
    @test sum(population.nᵢᵍ .* epi_params.CHᵢᵍ[:, :]) ≈ 528258.6803  atol=0.0001

end


## -----------------------------------------------------------------------------
## Epidemic spreading with multiple containments
## -----------------------------------------------------------------------------

tᶜs = [30, 60, 90, 120]
κ₀s = [0.65, 0.75, 0.65, 0.55]
ϕs = [0.174, 0.174, 0.174, 0.174]
δs = [0.207, 0.207, 0.207, 0.207]

reset_params!(epi_params, population)
set_initial_infected!(epi_params, population, E₀, A₀, I₀)
run_epidemic_spreading_mmca!(epi_params, population, tᶜs, κ₀s, ϕs, δs; verbose = false)


## Test epidemic spreading with multiple containments
@testset "Epidemic_Spreading_With_Multiple_Contaiments" begin

    @test sum(population.nᵢᵍ .* epi_params.ρˢᵍ[:, :, T]) ≈ 537528.5243  atol=0.0001
    @test sum(population.nᵢᵍ .* epi_params.ρᴱᵍ[:, :, T]) ≈ 4.1081  atol=0.0001
    @test sum(population.nᵢᵍ .* epi_params.ρᴬᵍ[:, :, T]) ≈ 5.2122  atol=0.0001
    @test sum(population.nᵢᵍ .* epi_params.ρᴵᵍ[:, :, T]) ≈ 6.9785  atol=0.0001
    @test sum(population.nᵢᵍ .* epi_params.ρᴾᴴᵍ[:, :, T]) ≈ 0.0647  atol=0.0001
    @test sum(population.nᵢᵍ .* epi_params.ρᴾᴰᵍ[:, :, T]) ≈ 0.2486  atol=0.0001
    @test sum(population.nᵢᵍ .* epi_params.ρᴴᴿᵍ[:, :, T]) ≈ 0.3839  atol=0.0001
    @test sum(population.nᵢᵍ .* epi_params.ρᴴᴰᵍ[:, :, T]) ≈ 0.0374  atol=0.0001
    @test sum(population.nᵢᵍ .* epi_params.ρᴰᵍ[:, :, T]) ≈ 269.0122  atol=0.0001
    @test sum(population.nᵢᵍ .* epi_params.ρᴿᵍ[:, :, T]) ≈ 14590.1274 atol=0.0001
    @test sum(population.nᵢᵍ .* epi_params.CHᵢᵍ[:, :]) ≈ 447595.3021  atol=0.0001

end


## -----------------------------------------------------------------------------
## Effective reproduction number R
## -----------------------------------------------------------------------------

τ = 21

Rᵢᵍ_eff, R_eff = compute_R_eff(epi_params, population, τ)


## Test effective reproduction number R
@testset "Effective_Reproduction_Number" begin

    R_eff_sample = [4.26608; 4.24034; 3.61001; 0.70638; 0.70111; 0.67445;
                    0.54267; 0.54263; 0.56545; 0.69465; 0.69531; 0.72088;
                    0.86301; 0.86293; 0.86287; 0.86281; 0.86276; 0.86272]

    Rᵢᵍ_eff_end = [0.525514  0.525466  0.525477  0.525607  0.525486
                   1.000240  1.000230  1.000210  1.000300  1.000220
                   0.695367  0.695334  0.695326  0.695431  0.695339]

    @test nrow(Rᵢᵍ_eff) == G * M * (T - τ)
    @test ncol(Rᵢᵍ_eff) == 4
    @test nrow(R_eff) == T - τ
    @test ncol(R_eff) == 2
    @test reshape(Rᵢᵍ_eff.R_eff, (G, M, T - τ))[:, :, T - τ] ≈ Rᵢᵍ_eff_end  atol=0.0001
    @test R_eff.R_eff[1:10:end] ≈ R_eff_sample  atol=0.0001

end

"""
Module that contains all the required material to analyze the epidemic spreading
of COVID-19 using MMCA reported in
"""
module MMCAcovid19

using Statistics, DataFrames, CSV, Printf

export run_epidemic_spreading_mmca!,
    Epidemic_Params,
    epidemic_params,
    Population_Params,
    population_params,
    set_initial_infected!,
    compute_R_eff,
    store_compartment,
    store_R_eff

# Load source from files
include("markov_aux.jl")
include("markov.jl")

end

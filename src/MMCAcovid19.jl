"""
Module that contains all the required material to analyze the epidemic spreading
of COVID-19 using MMCA reported in
"""
module MMCAcovid19

using Statistics, DataFrames, CSV, Printf

export run_epidemic_spreading_mmca!,
    Epidemic_Params,
    Population_Params,
    set_initial_infected!,
    reset_params!,
    compute_R_eff,
    store_compartment,
    store_R_eff

# Load source from files
include("markov_aux.jl")
include("markov.jl")

end

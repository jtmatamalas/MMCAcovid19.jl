"""
Module that contains all the required material to analyze the epidemic spreading
of COVID-19 using MMCA reported in
"""
module MMCAcovid19

using Statistics, DataFrames, CSV, Printf

export run_epidemic_spreading_mmca!,
    Epidemic_Params,
    reset_epi_param!,
    Population_Params,
    update_population_param!,
    setInitialInfected!,
    store_compartment,
    store_R0,
    correct_self_loops

# Load source from files
include("markov_aux.jl")
include("markov.jl")

end

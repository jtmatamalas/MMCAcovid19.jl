push!(LOAD_PATH, "../src/")
using Documenter, MMCAcovid19

makedocs(sitename = "MMCAcovid19.jl",
         pages = [
             "Overview" => "index.md",
             "Getting started" => "getting_started.md",
             "Library" => "library.md"
         ]
         )

deploydocs(
    repo = "github.com/jtmatamalas/MMCAcovid19.jl.git"
)

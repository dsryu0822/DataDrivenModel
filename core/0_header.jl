using Base.Threads: @threads, nthreads # Base.Threads.nthreads()
using Dates
device = gethostname() # In julia v1.11, it could be replaced by `Sys.username()`
@info "$(now()) - @$device $(nthreads()) threads"

import Pkg
function load_packages(packages)
    packages = string.(packages)
    try
        @time "$packages load" eval(Meta.parse("using $(join(packages, ", "))"))
    catch e
        # required = setdiff(packages, getproperty.(values(Pkg.dependencies()), :name))
        required = setdiff(packages, keys(Pkg.installed()))
        print("\r\x1b[F")
        if required |> !isempty
            @info "Installing required packages: $(required)"
            Pkg.add.(required)
        end
        @time "$packages load" eval(Meta.parse("using $(join(packages, ", "))"))
    end
end
load_packages([:CSV, :JLD2, :ProgressMeter])


# packages = [:Combinatorics, :LinearAlgebra, :SparseArrays, :DataFrames,
#             :PrettyTables, :Symbolics, :CSV, :DecisionTree, :Random,
#             :StatsBase, :Dates, :ProgressMeter, :Plots, :LaTeXStrings,
#             :Colors, :ColorSchemes, :JLD2, :Graphs,
#             :DifferentialEquations, :OrdinaryDiffEqLowOrderRK, :Sundials, :DiffEqBase
#             ] .|> string
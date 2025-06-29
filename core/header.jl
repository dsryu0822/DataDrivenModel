using Base.Threads: @threads, nthreads # Base.Threads.nthreads()
import Pkg
packages = [:Combinatorics, :LinearAlgebra, :SparseArrays, :DataFrames,
            :PrettyTables, :Symbolics, :CSV, :DecisionTree, :Random,
            :StatsBase, :Dates, :ProgressMeter, :Plots, :LaTeXStrings,
            :JLD2, :Graphs] .|> string
required = setdiff(packages, getproperty.(values(Pkg.dependencies()), :name))
if required |> !isempty
    @info "Installing required packages: $(required)"
    Pkg.add.(installed)
end
@time "All packages load" eval(Meta.parse("using $(join(packages, ", "))"))
mm = Plots.mm
cm = Plots.cm

device = gethostname() # In julia v1.11, it could be replaced by `Sys.username()`
if Sys.iswindows()
    if device != "Sickbook"
        cd("G:/DDM")
    end
    # device = ENV["COMPUTERNAME"];
elseif Sys.islinux()
    cd("/home/$(ENV["LOGNAME"])/g/DDM")
    # device = ENV["LOGNAME"];
end
@info "$(now()) - @$device $(nthreads()) threads"

include("../core/factorio.jl")
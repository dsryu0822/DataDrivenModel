using Base.Threads: @threads, nthreads # Base.Threads.nthreads()
using ProgressMeter
packages = [:Combinatorics, :LinearAlgebra, :SparseArrays, :DataFrames,
            :PrettyTables, :Symbolics, :CSV, :DecisionTree, :Random,
            :StatsBase, :Dates, :Plots]
@showprogress @threads for package in packages
    @eval using $(package)
end
mm = Plots.mm
cm = Plots.cm

if Sys.iswindows()
    cd("G:/DDM")
    # device = ENV["COMPUTERNAME"];
elseif Sys.islinux()
    cd("/home/$(ENV["LOGNAME"])/g/DDM")
    # device = ENV["LOGNAME"];
end
device = gethostname() # In julia v1.11, it could be replaced by `Sys.username()`
@info "$(now()) - @$device $(nthreads()) threads"

include("../core/factorio.jl")
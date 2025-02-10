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
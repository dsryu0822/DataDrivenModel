using ProgressMeter, CSV, DecisionTree, Random, StatsBase, Dates, Plots;
using Base.Threads: @threads, nthreads # Base.Threads.nthreads()

if Sys.iswindows()
    cd("G:/DDM")
    device = ENV["COMPUTERNAME"];
elseif Sys.islinux()
    cd("/home/$(ENV["LOGNAME"])/g/DDM")
    device = ENV["LOGNAME"];
end
@info "$(now()) - $device $(nthreads()) threads"

include("../core/factorio.jl")
using DecisionTree, Random, StatsBase, Dates; 
using Base.Threads: @threads, nthreads # Base.Threads.nthreads()

if Sys.iswindows()
    cd("G:/DDM")
    device = ENV["COMPUTERNAME"];
elseif Sys.islinux()
    cd("/home/$(ENV["LOGNAME"])/g/DDM")
    device = ENV["LOGNAME"];
end
@info "$(now()) - $device $(nthreads()) threads"

include("../core/DDM.jl")
include("../core/factorio.jl")
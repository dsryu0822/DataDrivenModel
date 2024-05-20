include("../core/DDM.jl")
include("../core/factorio.jl")
include("../core/visual.jl")
using DecisionTree, Random, StatsBase
using Base.Threads: @threads # Base.Threads.nthreads()

data = CSV.read("G:/DDM/bifurcation/buck/02501.csv", DataFrame)

Random.seed!(0)
schedules = DataFrame(idx = 1:100, ic = rand(nrow(data), 100))

vrbl = [:dV, :dI], [:V, :I]
cnfg = (; N = 1)
dt = 1e-7
θ1 = 1e+1; θ2 = 1e+0; θ3 = 1e+0; min_rank = 2;

@threads for dr in eachrow(schedules) # dr = eachrow(schedules)[end]
    filename = "G:/DDM/forecast/buck_rcvd/$(lpad(dr.idx, 5, '0')).csv"
    isfile(filename) && continue
    data = CSV.read(replace(filename, "buck_rcvd" => "buck"), DataFrame);

    add_subsystem!(data, vrbl, cnfg; (; θ1, θ2, θ3, min_rank)...); # 30 sec
    f_ = [SINDy(df, vrbl...; cnfg...) for df in groupby(data, :subsystem)];
    Dtree = dryad(data, [:V, :I, :Vr])
    # Dtree = dryad(data, last(vrbl))

    ic = collect(data[1, last(vrbl)])
    try
        ŷ = DataFrame(solve(f_, ic, dt, data.Vr, Dtree, data.Vr), last(vrbl))
        CSV.write(filename, ŷ)
    catch
        @error "Error: $(lpad(dr.idx, 5, '0'))"
    end
end

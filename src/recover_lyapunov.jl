using DecisionTree, Random, StatsBase, Dates; @info now(); sec = Second(1);
using Base.Threads: @threads, nthreads # Base.Threads.nthreads()
include("../core/DDM.jl")
include("../core/factorio.jl")
include("../core/visual.jl")


# canonicalize(2001*100*50sec)
# canonicalize((201*50÷20)*50sec)

# function initialize_soft()
    # schedules = CSV.read("G:/DDM/lyapunov/soft_schedules_cache.csv", DataFrame)
    # schedules[!, :λ1] .= .0; schedules[!, :λ2] .= .0; schedules[!, :λ3] .= .0;
    # for i = 1:2001
    #     dr = schedules[[i], :]
    #     CSV.write("G:/DDM/lyapunov/soft/$(lpad(i, 5, '0')).csv", dr, bom = true)
    # end
# end

function lyapunov_soft()
    schedules = CSV.read("G:/DDM/lyapunov/soft_schedules_cache.csv", DataFrame)
    vrbl = [:dt, :du, :dv], [:t, :u, :v]
    cnfg = (; f_ = [cospi, sign], λ = 1e-2)
    dt = 1e-6; θ1 = 1e-8; θ2 = 1e-12; θ3 = 1e-5; min_rank = 21;

    sync = 0
    @showprogress @threads for dr = eachrow(schedules)[1:100:end]
        # dr = eachrow(schedules)[1]
        
        filename = "G:/DDM/lyapunov/soft/$(lpad(dr.idx, 5, '0')).csv"
        !isfile(filename) && continue
        note = CSV.read(filename, DataFrame)
        data = CSV.read(replace(filename, "lyapunov" => "bifurcation"), DataFrame)
    
        # add_subsystem!(data, vrbl, cnfg; (; θ1, θ2, θ3, min_rank)...); # 30 sec
        f_ = [SINDy(df, vrbl...; cnfg...) for df in groupby(data, :subsystem)]
        Dtree = dryad(data, last(vrbl)); # print_tree(Dtree)
        data = nothing; GC.gc()

        try
            if sync ≤ nthreads()
                sleep(dr.idx/10); sync += 1 # It's technical trick for threading
                J_ = jacobian.(f_)
            end
        catch
            @error "Error in $(dr.idx)"
            continue
        end
        J_ = jacobian.(f_)
        
        # J = substitute(J_[1], Dict(t => note.t[1]))
        # U = Float64[0 0 -1; 0 1 0; -1 0 0]
        U, J = rand(3, 3), rand(3, 3)
        for _t = round(Int64, note.t[end]):60
            # _t = 50
            # println(dr)
            data = DataFrame(solve(f_, collect(note[end, 3:5]), dt, _t:dt:(_t+1), Dtree), last(vrbl))
            λ = collect(note[end, 6:end])
            for i in 1:nrow(data)
                U, V = gram_schmidt(U)
                λ += V |> eachcol .|> norm .|> log
                
                U = RK4(J, U, dt)
                s = apply_tree(Dtree, collect(data[i, 1:3]))
                J = Float64.(substitute(J_[s], Dict(t => data.t[i])))
            end
            push!(note, [dr.idx, dr.d, data[end, 1:3]..., λ...])
            CSV.write(filename, note, bom = true)
        end
    end
    open("G:/DDM/lyapunov/done!", "w") do io
        println(io, Dates.now())
    end
end
lyapunov_soft()
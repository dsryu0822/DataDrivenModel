using DecisionTree, Random, StatsBase, Dates; @info now()
using Base.Threads: @threads # Base.Threads.nthreads()
include("../core/DDM.jl")
include("../core/factorio.jl")
include("../core/visual.jl")

function lyapunov_soft()
    schedules = CSV.read("G:/DDM/lyapunov/soft_schedules_cache.csv", DataFrame)
    vrbl = [:dt, :du, :dv], [:t, :u, :v]
    cnfg = (; f_ = [cospi, sign], λ = 1e-2)
    dt = 1e-6; θ1 = 1e-8; θ2 = 1e-12; θ3 = 1e-5; min_rank = 21;

    result = DataFrame(d = Float64[], λ1 = Float64[], λ2 = Float64[], λ3 = Float64[])
    @showprogress for dr = eachrow(schedules)
    # for dr = eachrow(schedules)[1]
        try
            filename = "G:/DDM/bifurcation/soft/$(lpad(dr.idx, 5, '0')).csv"
            !isfile(filename) && continue
        
            data = CSV.read(filename, DataFrame)
            # add_subsystem!(data, vrbl, cnfg; (; θ1, θ2, θ3, min_rank)...); # 30 sec
            f_ = [SINDy(df, vrbl...; cnfg...) for df in groupby(data, :subsystem)]
            Dtree = dryad(data, last(vrbl)); # print_tree(Dtree)
            J_ = jacobian.(f_)

            ic = [50, dr.u, dr.v]
            data = DataFrame(solve(f_, ic, dt, 50:dt:150, Dtree), last(vrbl));

            λ = zeros(3);
            s = apply_tree(Dtree, collect(data[1, 1:3]))
            J = Float64.(substitute(J_[s], Dict(t => data.t[1])))
            U, _ = qr(J); U = Matrix(U)
            for i in 2:nrow(data)
                U, V = gram_schmidt(U)
                λ += V |> eachcol .|> norm .|> log          

                U = RK4(J, U, dt)
                s = apply_tree(Dtree, collect(data[i, 1:3]))
                J = Float64.(substitute(J_[s], Dict(t => data.t[i])))
            end

            λ ./= dt*nrow(data)
            push!(result, [dr.d, λ...])
            open("G:/DDM/lyapunov/soft.log", "a") do logfile
                println(logfile, now(), ",good,$(lpad(dr.idx, 5, '0'))")
            end
        catch
            open("G:/DDM/lyapunov/soft.log", "a") do logfile
                println(logfile, now(), ",error,$(lpad(dr.idx, 5, '0'))")
            end
        end
    end

    sort!(result, :d)
    CSV.write("G:/DDM/lyapunov/soft_lyapunov.csv", result, bom = true)
    plot(legend = :none, size = [600, 300])
    plot!(result.d, result.λ1, lw = 2, color = 1)
    plot!(result.d, result.λ2, lw = 2, color = 2)
    plot!(result.d, result.λ3, lw = 2, color = 3)
    png("G:/DDM/lyapunov/soft_lyapunov.png")
end
lyapunov_soft()
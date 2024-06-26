using DecisionTree, Random, StatsBase, Dates; @info now(); sec = Second(1);
using Base.Threads: @threads, nthreads # Base.Threads.nthreads()
include("../core/DDM.jl")
include("../core/factorio.jl")
include("../core/visual.jl")

if Sys.iswindows()
    cd("G:/DDM")
elseif Sys.islinux()
    cd("/home/$(ENV["LOGNAME"])/g/DDM")
end

# canonicalize(2001*100*50sec)
# canonicalize((201*50÷20)*50sec)

# function initialize_soft()
    # schedules = CSV.read("G:/DDM/lyapunov/soft_schedules_cache.csv", DataFrame)
    # schedules[!, :λ1] .= .0; schedules[!, :λ2] .= .0; schedules[!, :λ3] .= .0;
    # for i = 1:2001
    #     dr = schedules[[i], :]
    #     if !isfile("G:/DDM/lyapunov/soft/$(lpad(i, 5, '0')).csv")
    #         CSV.write("G:/DDM/lyapunov/soft/$(lpad(i, 5, '0')).csv", dr, bom = true)
    #     end
    # end
# end

# J_soft(t, u, v, d) = [
#                0                                0                                 0
#                0                                0                                 1
#      -π*sinpi(t) ifelse(abs(u) ≥ d/2, -160000, 0) ifelse(abs(u) ≥ d/2, -172.363, 0)
# ]

# print(f_[1])
# print(f_[2])
function lyapunov_soft()
    schedules = CSV.read("lyapunov/soft_schedules_cache.csv", DataFrame)
    vrbl = [:dt, :du, :dv], [:t, :u, :v]
    cnfg = (; f_ = [cospi, sign], λ = 1e-2)
    dt = 1e-6; θ1 = 1e-8; θ2 = 1e-12; θ3 = 1e-5; min_rank = 21;

    sync = 0
    # for dr = eachrow(schedules)[[1563]]
    # chaos1    -    1:500
    # chaos2    -  501:1000
    # chaos3    - 1001:1500
    # sicklinux - 1501:end 
    @showprogress @threads for dr = eachrow(schedules)[1:end]
        # dr = eachrow(schedules)[1563]
        
        filename = "lyapunov/soft/$(lpad(dr.idx, 5, '0')).csv"
        !isfile(filename) && continue
        note = CSV.read(filename, DataFrame); note[!, Not(1)] .= Float64.(note[!, Not(1)])
        data = CSV.read(replace(filename, "soft" => "big"), DataFrame)
        
        # data = factory_soft(DataFrame, dr.idx, dr.d; ic = collect(note[end, 3:5]), tspan = [0, 20], dt)
        # add_subsystem!(data, vrbl, cnfg; (; θ1, θ2, θ3, min_rank)...); # 30 sec

        f_ = [SINDy(df, vrbl...; cnfg...) for df in groupby(data, :subsystem)]
        Dtree = dryad(data, last(vrbl)); # print_tree(Dtree)
        data = nothing; GC.gc()

        try
            if sync ≤ nthreads()
                # sleep(dr.idx/10); sync += 1 # It's technical trick for threading
                # sleep(nthreads()*rand()) # It's technical trick for threading
                sleep(dr.idx % 109) # It's technical trick for threading
                J_ = jacobian.(f_)
            end
        catch
            @error "Error in $(dr.idx)"
            continue
        end
        J_ = jacobian.(f_)
        
        s = apply_tree(Dtree, collect(note[end, 3:5]))
        J = substitute(J_[s], Dict(t => note.t[end]))
        # J = J_soft(collect(note[1, 3:5])..., dr.d)
        
        U = reshape(collect(note[end, 6:(end-3)]), 3, 3)
        if U |> isone
            U, _ = qr(J); U = Matrix(U)
        end
        @info "$(now()) lyapunov start"
        for _t = round(Int64, note.t[end]):10:40
            
            data = DataFrame(solve(f_, collect(note[end, 3:5]), dt, _t:dt:(_t+10), Dtree), last(vrbl))
            # λ = collect(note[end, ((end-2):end)])
            # x = collect(note[end, 3:5])
            λ = zeros(3)
            for i in 1:nrow(data)
                U, V = gram_schmidt(U)
                λ += V |> eachcol .|> norm .|> log
                
                U = RK4(J, U, dt)
                # s = abs(data.u[i]) ≥ dr.d/2 ? 2 : 1
                s = apply_tree(Dtree, collect(data[i, 1:3]))
                J = Float64.(substitute(J_[s], Dict(t => data.t[i])))
                # 이게 더 중요한 것 같은데??

                # s = apply_tree(Dtree, x)
                # x, _ = RK4(f_[s], x, dt)
                # J = substitute(J_[s], Dict(t => x[1]))
                # s = apply_tree(Dtree, x)
                # J = Float64.(substitute(J_[s], Dict(t => x[1])))
                # x, _ = RK4(f_[s], x, dt)
            end; λ

            #############################
            # data = factory_soft(DataFrame, dr.idx, dr.d; ic = collect(note[end, 3:5]), tspan = [_t, _t+10], dt)[:, 1:3]
            # λ = zeros(3);
            # for i in 1:nrow(data)
            #     U, V = gram_schmidt(U)
            #     λ += V |> eachcol .|> norm .|> log

            #     U = RK4(J, U, dt)
            #     J = J_soft(data.t[i], data.u[i], data.v[i], dr.d)
            # end; λ
            #############################

            push!(note, [dr.idx, dr.d, data[end, 1:3]..., U..., λ...])
            CSV.write(filename, note, bom = true)
            CSV.write(string(@__DIR__) * "/../archive/$(lpad(dr.idx, 5, '0')).csv", note, bom = true)
        end
    end
    open("lyapunov/$(ENV["COMPUTERNAME"]) done!", "w") do io
        println(io, Dates.now())
    end
end
lyapunov_soft()

# plot(_data.t[1:1000:end], _data.v[1:1000:end]); png("temp1")
# plot( data.t[1:1000:end],  data.v[1:1000:end]); png("temp2")


# X = [0	-0.464422967	-0.885613985;	0	-0.88561352	0.464422082;	1	5.52e-17	1.10e-10]
# Y = [0	-0.493610779	-0.869683424;	0	-0.86968293	0.49360991;	1	9.16E-19	1.83E-12]

# qr(X).R |> eachrow .|> norm .|> log
# qr(Y).R |> eachrow .|> norm .|> log

# _λ
# λ
# findall(s1 .!= s2)

# function prepare_soft()
#     schedules = CSV.read("lyapunov/soft_schedules_cache.csv", DataFrame)
#     vrbl = [:dt, :du, :dv], [:t, :u, :v]
#     cnfg = (; f_ = [cospi, sign], λ = 1e-2)
#     dt = 1e-6; θ1 = 1e-8; θ2 = 1e-12; θ3 = 1e-5; min_rank = 21;

#     @showprogress @threads for dr = eachrow(schedules)
#         filename = "lyapunov/soft/$(lpad(dr.idx, 5, '0')).csv"
#         note = CSV.read(filename, DataFrame); note[!, Not(1)] .= Float64.(note[!, Not(1)])
#         data = factory_soft(DataFrame, dr.idx, dr.d; ic = collect(note[end, 3:5]), tspan = [0, 10], dt)
#         add_subsystem!(data, vrbl, cnfg; (; θ1, θ2, θ3, min_rank)...); # 30 sec
#         CSV.write("G:/DDM/lyapunov/big/$(lpad(dr.idx, 5, '0')).csv", data, bom = true)
#     end
# end
# prepare_soft()
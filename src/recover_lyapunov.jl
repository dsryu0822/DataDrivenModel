include("../core/header.jl")

# sec = Second(1)
# canonicalize(2001*100*50sec)
# canonicalize((201*50÷20)*50sec)


function lyapunov_exponent(_data, J_, bf_param;
    # U = Matrix(qr(J_(collect(_data[1, :])..., bf_param)).Q),
    U = I(ncol(_data)))

    λ = zeros(size(U, 1))
    for k = 1:nrow(_data)
        J = J_(collect(_data[k, :])..., bf_param)
        U, V = gram_schmidt(U)
        if 2k > nrow(_data)
            λ += V |> eachcol .|> norm .|> log
        end
        U = RK4(J, U, dt)
    end
    return sort(2λ / (last(_data.t) - first(_data.t)), rev=true)
end

# function initialize_soft()
    # schedules = CSV.read("lyapunov/soft_schedules_cache.csv", DataFrame)
    # schedules[!, :λ1] .= .0; schedules[!, :λ2] .= .0; schedules[!, :λ3] .= .0;
    # for i = 1:2001
    #     dr = schedules[[i], :]
    #     if !isfile("lyapunov/soft/$(lpad(i, 5, '0')).csv")
    #         CSV.write("lyapunov/soft/$(lpad(i, 5, '0')).csv", dr, bom = true)
    #     end
    # end
# end


function lyapunov_soft()
    schedules = CSV.read("lyapunov/soft_schedules_cache.csv", DataFrame)
    vrbl = [:dt, :du, :dv], [:t, :u, :v]
    cnfg = (; f_ = [cospi, sign], λ = 1e-2)
    dt = 1e-6; θ1 = 1e-8; θ2 = 1e-12; θ3 = 1e-5; min_rank = 21;

    tend = 90
    sync = 0
    # chaos1    -    1:500
    # chaos2    -  501:1000
    # chaos3    - 1001:1500
    # sicklinux - 1501:end
    @showprogress for dr = eachrow(schedules)
        # dr = eachrow(schedules)[1563]
        
        filename = "lyapunov/soft/$(lpad(dr.idx, 5, '0')).csv"
        note = CSV.read(filename, DataFrame); note[!, Not(1)] .= Float64.(note[!, Not(1)])
        length(note.t) ≥ 11 && continue
        data = CSV.read(replace(filename, "soft" => "big"), DataFrame)
        
        # data = factory_soft(DataFrame, dr.d; ic = collect(note[end, 3:5]), tspan = [0, 20], dt)
        # add_subsystem!(data, vrbl, cnfg; (; θ1, θ2, θ3, min_rank)...); # 30 sec

        f_ = [SINDy(df, vrbl...; cnfg...) for df in groupby(data, :subsystem)]
        Dtree = dryad(data, last(vrbl)); # print_tree(Dtree)
        data = nothing; GC.gc()

        
        sleep(rand()); sync += 1
        try
            if sync ≤ nthreads()
                sleep((dr.idx % 109) + rand(1:nthreads())) # It's technical trick for threading
                J_ = jacobian.(f_)
            end
        catch
            @error "Error in $(dr.idx)"
            continue
        end
        J_ = jacobian.(f_)
        
        s = apply_tree(Dtree, collect(note[end, 3:5]))
        J = substitute(J_[s], Dict(t => note.t[end]))
        
        U = reshape(collect(note[end, 6:(end-3)]), 3, 3)
        if U |> isone
            U, _ = qr(J); U = Matrix(U)
        end
        for _t = round(Int64, note.t[end]):10:90
            
            data = DataFrame(solve(f_, collect(note[end, 3:5]), dt, _t:dt:(_t+10), Dtree), last(vrbl))
            λ = zeros(3)
            for i in 1:nrow(data)
                U, V = gram_schmidt(U)
                λ += V |> eachcol .|> norm .|> log
                
                U = RK4(J, U, dt)
                s = apply_tree(Dtree, collect(data[i, 1:3]))
                J = Float64.(substitute(J_[s], Dict(t => data.t[i])))
            end; λ

            push!(note, [dr.idx, dr.d, data[end, 1:3]..., U..., λ...])
            CSV.write(filename, note, bom = true)
            CSV.write(string(@__DIR__) * "/../archive/$(lpad(dr.idx, 5, '0')).csv", note, bom = true)
        end
    end
    open("lyapunov/$(device) done!", "w") do io
        println(io, Dates.now())
    end
end
lyapunov_soft()



##########################################################################
#                                                                        #
#                           Hindmarsh-Rose model                         #
#                                                                        #
##########################################################################
schedules = CSV.read("bifurcation/hrnm_schedules.csv", DataFrame)
schedules[!, :λ1] .= .0; schedules[!, :λ2] .= .0; schedules[!, :λ3] .= .0; schedules[!, :λ4] .= .0;
vrbl = [:dt, :dx, :dy, :dz], [:t, :x, :y, :z]
eval(Meta.parse("@variables $(join(string.(last(vrbl)), " "))"))
cnfg = (; N = 3, f_ = [cos])

@showprogress @threads for dr = eachrow(schedules)
        filename = "lyapunov/hrnm_traj/$(lpad(dr.idx, 5, '0')).csv"
        data = CSV.read(filename, DataFrame)
        
        f_ = [SINDy(df, vrbl...; cnfg...) for df in groupby(data, :subsystem)]
        J_ = jacobian.(f_)
        Dtree = dryad(data, last(vrbl)); # print_tree(Dtree)
        s = apply_tree(Dtree, collect(note[end, 3:5]))

        λ = lyapunov_exponent(data[:, last(vrbl)], J_, dr.f)
        dr[[:λ1, :λ2, :λ3, :λ4]] .= λ
end
CSV.write("lyapunov/!linux hrnm_lyapunov 1e-3 t = [0, 4000].csv", schedules, bom = true)

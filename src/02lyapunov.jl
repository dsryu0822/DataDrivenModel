include("../core/header.jl")

function lyapunov_exponent(_data::DataFrame, J_, DT::Root{Float64, Int64}, bf_param;
    U = I(ncol(_data)), T = (last(_data.t) - first(_data.t)))

    λ = zeros(size(U, 1))
    for k = 1:nrow(_data)
        s = apply_tree(DT, collect(_data[k, :]))
        J = J_[s](collect(_data[k, :]))
        U, V = gram_schmidt(U)
        λ += V |> eachcol .|> norm .|> log
        U = RK4(J, U, dt)
    end
    return sort(λ / T, rev=true)
end

##########################################################################
#                                                                        #
#                            Soft impact model                           #
#                                                                        #
##########################################################################
# function J_(t, u, v, d)
#     return [   0                                0                                 0
#                0                                0                                 1
#      -π*sinpi(t) ifelse(abs(u) ≥ d/2, -160000, 0) ifelse(abs(u) ≥ d/2, -172.363, 0) ]
# end
# idx_tgt = parse.(Int64, first.(readdir("data/soft"), 5))
# schedules = CSV.read("schedules/soft.csv", DataFrame)[idx_tgt, :]
# schedules = schedules[1:1:end, :]
schedules = CSV.read("schedules/soft.csv", DataFrame)[1:1:end, :]
schedules[!, :λ1] .= .0; schedules[!, :λ2] .= .0; schedules[!, :λ3] .= .0;
vrbl = [:dt, :du, :dv], [:t, :u, :v]
cnfg = (; f_ = [cospi], λ = 5e-1) # λ = 5e-1 → 1e-2
dt = 1e-5; θ = 1e-6;

@showprogress @threads for dr = eachrow(schedules)
    filename = "data/soft/$(lpad(dr.idx, 5, '0')).csv"
    data = CSV.read(filename, DataFrame)
    
    f_ = [SINDy(df, vrbl...; cnfg...) for df in groupby(data, :subsystem)]
    Dtree = dryad(data, last(vrbl)); # print_tree(Dtree)
    J_ = []
    while true
        try
            J_ = jacobian.(Function, f_)
            break
        catch
            print(".")
            sleep(rand())
        end
    end

    # Initial condition이 [first(data)...][1:3]이 되면 오히려 이상해짐
    # data = DataFrame(solve(f_, [first(data)...][1:3], dt, 0:dt:100, Dtree), last(vrbl))
    data = DataFrame(solve(f_, [eps(), .05853, .47898], dt, 0:dt:150, Dtree), last(vrbl))
    λ = lyapunov_exponent(data[:, last(vrbl)], J_, Dtree, dr.bp)
    dr[[:λ1, :λ2, :λ3]] .= λ
    CSV.write("output/...$(device)ing lpnv_soft.csv", schedules, bom = true)
end
CSV.write("output/!$(device) lpnv_soft.csv", schedules, bom = true)


##########################################################################
#                                                                        #
#                             Gear system                                #
#                                                                        #
##########################################################################
∃ = parse.(Int64, replace.(readdir("data/gear"), ".csv" => ""))
schedules = CSV.read("schedules/gear.csv", DataFrame)[∃, :]
schedules[!, :λ1] .= .0; schedules[!, :λ2] .= .0; schedules[!, :λ3] .= .0; schedules[!, :λ4] .= .0;
vrbl = [:dx, :dv, :dΩ, :dθ], [:x, :v, :Ω, :θ]
cnfg = (; N = 1, f_ = [cos], C = 2,  λ = 1e-2)
# λ = 1e-2 works for 762 bps/ λ = 1e-4 works for 138 bps
dt = 1e-2; tspan = [0, 1000]; θ = 1e-4; dos = 1
smse(SINDy_) = sum(getproperty.(SINDy_, :MSE))

# dr = eachrow(schedules)[399]
# bfcn = DataFrame(hrzn = [], vrtc = [])
@showprogress @threads for dr = eachrow(schedules)
    filename = "data/gear/$(lpad(dr.idx, 5, '0')).csv"
    data = CSV.read(filename, DataFrame)
    
    # add_subsystem!(data, vrbl, cnfg; θ1, θ2, θ3, min_rank); # 30 sec
    f1_ = [SINDy(df, vrbl...; N = 1, f_ = [cos], C = 2,  λ = 1e-2) for df in groupby(data, :subsystem)]
    f2_ = [SINDy(df, vrbl...; N = 1, f_ = [cos], C = 2,  λ = 1e-4) for df in groupby(data, :subsystem)]
    f_ = ifelse(smse(f1_) < smse(f2_), f1_, f2_)
    # f_ = [SINDy(df, vrbl...; cnfg...) for df in groupby(data, :subsystem)]
    Dtree = dryad(data, last(vrbl)); # print_tree(Dtree)
    J_ = []
    while true
        try
            J_ = jacobian.(Function, f_)
            break
        catch
            print(".")
        end
    end

    data = DataFrame(solve(f_, [0.1, 0.1, 0.1, eps()], dt, 0:dt:last(tspan), Dtree), last(vrbl))
    λ = lyapunov_exponent(data[:, last(vrbl)], J_, Dtree, dr.bp, T = last(tspan))
    dr[[:λ1, :λ2, :λ3, :λ4]] .= λ
    CSV.write("output/...$(device)ing gear_lyapunov_rcvd.csv", schedules, bom = true)
end
CSV.write("output/!$(device) lpnv_soft.csv", schedules, bom = true)
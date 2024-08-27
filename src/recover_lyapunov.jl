include("../core/header.jl")
# sec = Second(1)
# canonicalize(2001*100*50sec)
# canonicalize((201*50÷20)*50sec)

function lyapunov_exponent(_data::DataFrame, J_::Vector{Matrix{Num}}, DT::Root{Float64, Int64}, bf_param;
    U = I(ncol(_data)), T = (last(_data.t) - first(_data.t)))

    λ = zeros(size(U, 1))
    for k = 1:nrow(_data)
        s = apply_tree(DT, collect(_data[k, :]))
        J = Float64.(substitute(J_[s], Dict(t => _data.t[k])))
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
schedules = CSV.read("bifurcation/soft_schedules.csv", DataFrame)[1:10:end, :]
schedules[!, :λ1] .= .0; schedules[!, :λ2] .= .0; schedules[!, :λ3] .= .0;
vrbl = [:dt, :du, :dv], [:t, :u, :v]
cnfg = (; f_ = [cospi, sign], λ = 1e-2)
dt = 1e-6; θ1 = 1e-8; θ2 = 1e-12; θ3 = 1e-5; min_rank = 21;

@showprogress @threads for dr = eachrow(schedules)
    filename = "lyapunov/soft_traj/$(lpad(dr.idx, 5, '0')).csv"
    data = CSV.read(filename, DataFrame)
    
    # add_subsystem!(data, vrbl, cnfg; θ1, θ2, θ3, min_rank); # 30 sec
    f_ = [SINDy(df, vrbl...; cnfg...) for df in groupby(data, :subsystem)]
    Dtree = dryad(data, last(vrbl)); # print_tree(Dtree)
    J_ = []
    while true
        try
            J_ = jacobian.(f_)
            break
        catch
            print(".")
        end
    end

    data = DataFrame(solve(f_, [eps(), .05853, .47898], dt, 0:dt:150, Dtree), last(vrbl))
    λ = lyapunov_exponent(data[:, last(vrbl)], J_, Dtree, dr.d)
    dr[[:λ1, :λ2, :λ3]] .= λ
    CSV.write("lyapunov/!$(device)ing soft_lyapunov_rcvd.csv", schedules, bom = true)
end
CSV.write("lyapunov/!$(device) soft_lyapunov_rcvd.csv", schedules, bom = true)
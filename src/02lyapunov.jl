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
idx_tgt = parse.(Int64, first.(readdir("output/soft"), 5))
schedules = CSV.read("schedules/soft.csv", DataFrame)[idx_tgt, :]
schedules = schedules[1:10:end, :]
schedules[!, :λ1] .= .0; schedules[!, :λ2] .= .0; schedules[!, :λ3] .= .0;
vrbl = [:dt, :du, :dv], [:t, :u, :v]
cnfg = (; f_ = [cospi, sign], λ = 5e-1) # λ = 5e-1 → 1e-2
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

    # data = DataFrame(solve(f_, [first(data)...][1:3], dt, 0:dt:100, Dtree), last(vrbl))
    data = DataFrame(solve(f_, [eps(), .05853, .47898], dt, 0:dt:100, Dtree), last(vrbl))
    λ = lyapunov_exponent(data[:, last(vrbl)], J_, Dtree, dr.bp)
    dr[[:λ1, :λ2, :λ3]] .= λ
    CSV.write("output/...$(device)ing lpnv_soft.csv", schedules, bom = true)
end
CSV.write("output/!$(device) lpnv_soft.csv", schedules, bom = true)
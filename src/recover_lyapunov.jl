include("../core/header.jl")
# sec = Second(1)
# canonicalize(2001*100*50sec)
# canonicalize((201*50÷20)*50sec)

function lyapunov_exponent(_data::DataFrame, J_::Vector{Matrix{Num}}, DT::Root{Float64, Int64}, bf_param;
    U = I(ncol(_data)), T = (last(_data.t) - first(_data.t)))
    last_vrbl = Symbol.(names(_data))

    λ = zeros(size(U, 1))
    for k = 1:nrow(_data)
        s = apply_tree(DT, collect(_data[k, :]))
        String_Dict = "substitute($(J_[s]), Dict($(join(["$vb => $(_data[k, vb])" for vb in last_vrbl], ", "))))"
        J = Float64.(eval(Meta.parse(String_Dict)))
        U, V = gram_schmidt(U)
        λ += V |> eachcol .|> norm .|> log
        U = RK4(J, U, dt)
    end
    return sort(λ / T, rev=true)
end

# ##########################################################################
# #                                                                        #
# #                            Soft impact model                           #
# #                                                                        #
# ##########################################################################
# schedules = CSV.read("bifurcation/soft_schedules.csv", DataFrame)[1:10:end, :]
# schedules[!, :λ1] .= .0; schedules[!, :λ2] .= .0; schedules[!, :λ3] .= .0;
# vrbl = [:dt, :du, :dv], [:t, :u, :v]
# cnfg = (; f_ = [cospi, sign], λ = 1e-2)
# dt = 1e-6; θ1 = 1e-8; θ2 = 1e-12; θ3 = 1e-5; min_rank = 21;

# @showprogress @threads for dr = eachrow(schedules)
#     filename = "lyapunov/soft_traj/$(lpad(dr.idx, 5, '0')).csv"
#     data = CSV.read(filename, DataFrame)
    
#     # add_subsystem!(data, vrbl, cnfg; θ1, θ2, θ3, min_rank); # 30 sec
#     f_ = [SINDy(df, vrbl...; cnfg...) for df in groupby(data, :subsystem)]
#     Dtree = dryad(data, last(vrbl)); # print_tree(Dtree)
#     J_ = []
#     while true
#         try
#             J_ = jacobian.(f_)
#             break
#         catch
#             print(".")
#         end
#     end

#     data = DataFrame(solve(f_, [eps(), .05853, .47898], dt, 0:dt:150, Dtree), last(vrbl))
#     λ = lyapunov_exponent(data[:, last(vrbl)], J_, Dtree, dr.bp)
#     dr[[:λ1, :λ2, :λ3]] .= λ
#     CSV.write("lyapunov/!$(device)ing soft_lyapunov_rcvd.csv", schedules, bom = true)
# end
# CSV.write("lyapunov/!$(device) soft_lyapunov_rcvd.csv", schedules, bom = true)


# ##########################################################################
# #                                                                        #
# #                           Hindmarsh-Rose model                         #
# #                                                                        #
# ##########################################################################
# schedules = CSV.read("bifurcation/hrnm_schedules.csv", DataFrame)
# schedules[!, :λ1] .= .0; schedules[!, :λ2] .= .0; schedules[!, :λ3] .= .0; schedules[!, :λ4] .= .0;
# vrbl = [:dt, :dx, :dy, :dz], [:t, :x, :y, :z]
# cnfg = (; N = 3, f_ = [cos])
# dt = 1e-3; θ1 = 1e-2; θ2 = 1e-27; θ3 = 1e-1; min_rank = 32;

# # for dr = eachrow(schedules)[1:250]
# #     try
# #         filename = "bifurcation/hrnm/$(lpad(dr.idx, 5, '0')).csv"
# #         data = CSV.read(filename, DataFrame)
# #         MM = maximum(data.subsystem)
# #         if MM > 3
# #             print("$(dr.idx)... $(MM) → ")
# #             add_subsystem!(data, vrbl, cnfg; θ1, θ2, θ3, min_rank)
# #             CSV.write(filename, data)
# #             println(maximum(data.subsystem))
# #         end
# #     catch
# #         @error "$(dr.dix)"
# #     end
# # end
# # scatter(tmp)

# @showprogress @threads for dr = eachrow(schedules)
#     try
#     filename = "bifurcation/hrnm/$(lpad(dr.idx, 5, '0')).csv"
#     data = CSV.read(filename, DataFrame)

#     f_ = [SINDy(df, vrbl...; cnfg...) for df in groupby(data, :subsystem)]
#     Dtree = dryad(data, last(vrbl)); # print_tree(Dtree)
#     J_ = []
#     while true
#         try
#             J_ = jacobian.(f_)
#             break
#         catch
#             print(".")
#         end
#     end

#     data = DataFrame(solve(f_, [eps(), eps(), eps(), 0.1], dt, 0:dt:10000, Dtree), last(vrbl))
#     λ = lyapunov_exponent(data[:, last(vrbl)], J_, Dtree, dr.bp)
#     dr[[:λ1, :λ2, :λ3, :λ4]] .= λ
#     # CSV.write("lyapunov/!$(device)ing hrnm_lyapunov.csv", schedules, bom = true)
#     CSV.write("lyapunov/6채우기 hrnm_lyapunov.csv", schedules, bom = true)
# catch
#     @error "Error in $(dr.idx)"
# end
# end
# CSV.write("lyapunov/6채우기 hrnm_lyapunov.csv", schedules, bom = true)
# CSV.write("lyapunov/!$(device) hrnm_lyapunov.csv", schedules, bom = true)

##########################################################################
#                                                                        #
#                             Gear system                                #
#                                                                        #
##########################################################################
# data = factory_gear(DataFrame, -.09; dt = 1e-2, tspan = [500, 5000])
# plot(data.x, data.v)
# idx_sampled = diff([0; mod.(data.θ, 2π)]) .< 0
# data = data[idx_sampled, :]
# scatter(data.x, data.v, ms = 0.5, legend = :none)
schedules = CSV.read("bifurcation/gear_schedules.csv", DataFrame)
schedules[!, :λ1] .= .0; schedules[!, :λ2] .= .0; schedules[!, :λ3] .= .0; schedules[!, :λ4] .= .0;
vrbl = [:dx, :dv, :dΩ, :dθ], [:x, :v, :Ω, :θ]
cnfg = (; N = 1, f_ = [cos], C = 2,  λ = 1e-2)
dt = 1e-2; tend = 10000; θ1 = 7e-9; θ2 = 1e-12; θ3 = 1e-5; min_rank = 31; dos = 1

bfcn = DataFrame(hrzn = [], vrtc = [])
@showprogress @threads for dr = eachrow(schedules)
    filename = "bifurcation/gear/$(lpad(dr.idx, 5, '0')).csv"
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

    data = DataFrame(solve(f_, [0.1, 0.1, 0.1, eps()], dt, 0:dt:tend, Dtree), last(vrbl))
    λ = lyapunov_exponent(data[:, last(vrbl)], J_, Dtree, dr.bp, T = tend)
    dr[[:λ1, :λ2, :λ3, :λ4]] .= λ
    CSV.write("lyapunov/!$(device)ing gear_lyapunov_rcvd.csv", schedules, bom = true)
    
    data = data[data.Ω .> 5000, :]
    idx_sampled = diff([0; mod.(data.Ω, 2π)]) .< 0
    sampledx = data.v[idx_sampled]
    hrzn, vrtc = fill(dr.bp, length(sampledx)), sampledx
    append!(bfcn, DataFrame(; hrzn, vrtc))
    
    CSV.write("lyapunov/!$(device)ing gear_lyapunov_rcvd.csv", schedules, bom = true)
    CSV.write("lyapunov/!$(device)ing gear_bifurcation_rcvd.csv", bfcn, bom = true)
end
CSV.write("lyapunov/!$(device) gear_lyapunov_rcvd.csv", schedules, bom = true)
CSV.write("lyapunov/!$(device) gear_bifurcation_rcvd.csv", bfcn, bom = true)
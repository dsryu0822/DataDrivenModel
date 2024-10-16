include("../core/header.jl")

function lyapunov_exponent(_data::DataFrame, J_::Function, bf_param;
    U = I(ncol(_data)), T = (last(_data.t) - first(_data.t)))

    λ = zeros(size(U, 1))
    for k = 1:nrow(_data)
        J = J_(collect(_data[k, :])..., bf_param)
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
# # schedules = CSV.read("bifurcation/soft_schedules.csv", DataFrame)
# # schedules[!, :λ1] .= .0; schedules[!, :λ2] .= .0; schedules[!, :λ3] .= .0;
# vrbl = [:dt, :du, :dv], [:t, :u, :v]
# cnfg = (; f_ = [cospi, sign], λ = 1e-2)
# dt = 1e-6; θ1 = 1e-8; θ2 = 1e-12; θ3 = 1e-5; min_rank = 21;
# # dt = 1e-5 is good to bifurcation diagram but not for lyapunov spectrum, 1e-6 is needed
# function J_(t, u, v, d)
#     return [   0                                0                                 0
#                0                                0                                 1
#      -π*sinpi(t) ifelse(abs(u) ≥ d/2, -160000, 0) ifelse(abs(u) ≥ d/2, -172.363, 0) ]
# end

# schedules = CSV.read("lyapunov/linux_soft.csv", DataFrame)
# @showprogress @threads for dr = eachrow(schedules)[iszero.(schedules.λ1)]
# # @showprogress @threads for dr = eachrow(schedules)
#     # filename = "bifurcation/soft/$(lpad(dr.idx, 5, '0')).csv"
#     # data = CSV.read(filename, DataFrame)
#     data = factory_soft(DataFrame, dr.bp, tspan = [0, 150]; dt)
#     # add_subsystem!(data, vrbl, cnfg; θ1, θ2, θ3, min_rank); # 30 sec
#     # CSV.write(filename, data)

#     λ = lyapunov_exponent(data[:, last(vrbl)], J_, dr.bp)
#     dr[[:λ1, :λ2, :λ3]] .= λ
#     # idx_sampled = diff(abs.(data.u) .> (dr.bp/2)) .> 0
#     # sampledv = data[Not(1), :v][idx_sampled]
#     # append!(hrzn, fill(dr.bp, length(sampledv)))
#     # append!(vrtc, sampledv)
#     CSV.write("lyapunov/...$(device)ing soft t = [0, 150].csv", schedules, bom = true)
# end
# CSV.write("lyapunov/...$(device)ing soft t = [0, 150].csv", schedules, bom = true)


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
schedules = schedules[.!isfile.(["bifurcation/gear/$(lpad(dr.idx, 5, '0')).csv" for dr in eachrow(schedules)]), :]
schedules[!, :λ1] .= .0; schedules[!, :λ2] .= .0; schedules[!, :λ3] .= .0; schedules[!, :λ4] .= .0;
vrbl = [:dx, :dv, :dΩ, :dθ], [:x, :v, :Ω, :θ]
cnfg = (; N = 1, f_ = [cos], C = 2,  λ = 1e-4)
# λ = 1e-2 works for 762 bps/ λ = 1e-4 works for 138 bps
dt = 1e-2; tend = 1000; θ1 = 7e-9; θ2 = 1e-12; θ3 = 1e-5; min_rank = 31; dos = 1
function J_(x, v, Ω, θ, Fe)
    dfdx = ifelse(abs(x) > 1, 1, 0)
    return [ 0 1 0 0
             -(1 + k1*cos(Ω))*dfdx -2ζ (-Fe*H*sin(Ω) + k1*sin(Ω)*dfdx) -Fe*H*sin(θ)
             0 0 0 0
             0 0 0 0 ]
end

# bfcn = DataFrame(hrzn = [], vrtc = [])
@showprogress @threads for dr = eachrow(schedules)[381:420]
    data = factory_gear(DataFrame, dr.bp; tspan = [0, tend] .+ 500)
    # data = CSV.read("bifurcation/gear/$(lpad(dr.idx, 5, '0')).csv", DataFrame)
    add_subsystem!(data, vrbl, cnfg; θ1, θ2, θ3, min_rank, dos)
    CSV.write("bifurcation/gear/$(lpad(dr.idx, 5, '0')).csv", data, bom = true)


    # λ = lyapunov_exponent(data[:, last(vrbl)], J_, dr.bp, T = tend)
    # dr[[:λ1, :λ2, :λ3, :λ4]] .= λ

    # data = data[data.Ω .> 5000, :]
    # idx_sampled = diff([0; mod.(data.Ω, 2π)]) .< 0
    # sampledx = data.v[idx_sampled]
    # hrzn, vrtc = fill(dr.bp, length(sampledx)), sampledx
    # append!(bfcn, DataFrame(; hrzn, vrtc))
end
# CSV.write("lyapunov/gear_bifurcation.csv", bfcn, bom = true)
# CSV.write("lyapunov/gear_lyapunov.csv", schedules, bom = true)

# scatter(bfcn.hrzn, bfcn.vrtc, legend = false, alpha = .5, ms = .1, xlabel = "Fe", ylabel = "v")
# png("lyapunov/gear_bifurcation.png")

# plot(legend = :none)
# plot!(schedules.Fe, schedules.λ1)
# plot!(schedules.Fe, schedules.λ2)
# plot!(schedules.Fe, schedules.λ3)
# plot!(schedules.Fe, schedules.λ4)
# png("lyapunov/gear_lyapunov.png")

# ##########################################################################
# #                                                                        #
# #                           Hindmarsh-Rose model                         #
# #                                                                        #
# ##########################################################################
# # schedules = CSV.read("bifurcation/hrnm_schedules.csv", DataFrame)
# schedules[!, :λ1] .= .0; schedules[!, :λ2] .= .0; schedules[!, :λ3] .= .0; schedules[!, :λ4] .= .0;
# vrbl = [:dt, :dx, :dy, :dz], [:t, :x, :y, :z]
# cnfg = (; N = 3, f_ = [cos])
# dt = 1e-3; θ1 = 1e-2; θ2 = 1e-27; θ3 = 1e-1; min_rank = 32;
# function J_(t, x, y, z, _f)
#     _a = 1.0
#     _b = 3.0
#     _c = 1.0
#     _d = 5.0
#     _k = 0.9
#     # _f = 0.1
#     _ω = 1.0
#     _α = 0.1
#     _β = 0.8 
#     return [                0                             0   0    0
#              -_ω*_f*sin(_ω*t) (-3*_a*(x^2) + 2*_b*x + _k*z)   1 _k*x
#                             0                       -2*_d*x  -1    0
#                             0                            _β   0  -_α]
# end

# hrzn, vrtc = Dict(), Dict()
# @showprogress @threads for dr = eachrow(schedules)
#         # filename = "bifurcation/hrnm/$(lpad(dr.idx, 5, '0')).csv"
#         # filename = "lyapunov/hrnm_traj/$(lpad(dr.idx, 5, '0')).csv"
#         # data = CSV.read(filename, DataFrame)
#         # data = factory_hrnm(DataFrame, dr.bp; tspan = [0, 1000], dt)
#         data = factory_hrnm(DataFrame, dr.bp; tspan = [0, 10000], dt)
#         # add_subsystem!(data, vrbl, cnfg; θ1, θ2, θ3, min_rank)
#         # CSV.write(replace(filename, "bifurcation/hrnm" => "lyapunov/hrnm_traj"), data)
#         # CSV.write(filename, data)

#         λ = lyapunov_exponent(data[:, last(vrbl)], J_, dr.bp)
#         dr[[:λ1, :λ2, :λ3, :λ4]] .= λ

#         # data = data[data.t .> 9000, :]
#         # idx_sampled = abs.(diff(data.dz)) .> 0.1
#         # sampledx = data[Not(1), :x][idx_sampled]
#         # push!(hrzn, dr.idx => fill(dr.bp, length(sampledx)))
#         # push!(vrtc, dr.idx => sampledx)
# end
# bfcn = DataFrame(hrzn = vcat(values(hrzn)...), vrtc = vcat(values(vrtc)...))
# CSV.write("lyapunov/hrnm_lyapunov.csv", schedules, bom = true)
# CSV.write("lyapunov/hrnm_bifurcation.csv", bfcn, bom = true)

# @time eg_period = factory_hrnm(DataFrame, 0.14; tspan = [0, 10000])
# lyapunov_exponent(eg_period[:, last(vrbl)], J_, 0.14)
# plot(eg_period.x[1:100:end], eg_period.y[1:100:end], eg_period.z[1:100:end])
# scatter!(eg_period.x[[end]], eg_period.y[[end]], eg_period.z[[end]], shape = :x)


# ##########################################################################
# #                                                                        #
# #                               Lozi map                                 #
# #                                                                        #
# ##########################################################################
# b = 0.3
# hrzn = []; vrtc = []
# for a = 0.5:0.001:1.8
#     Lozi(x) = [1 + x[2] - a*abs(x[1]), b*x[1]]
#     v_ = [[.0, .0]]
#     for k in 1:1000
#         push!(v_, Lozi(v_[end]))
#     end
#     push!(hrzn, fill(a, 500))
#     push!(vrtc, first.(v_[(end-500):end]))
# end
# p1 = plot(legend = :none, xlabel = "a", ylabel = "x", title = "Lozi map")
# for ak in eachindex(hrzn)
#     scatter!(p1, hrzn[ak], vrtc[ak], color = :black, ms = 0.5)
# end
# plot!(p1, xlims = [0.5, 1.5]);
# png("Lozi_map.png")


@info "----------------------------------------------------------------------------------"

# ##########################################################################
# #                                                                        #
# #                            Lorenz attractor                            #
# #                                                                        #
# ##########################################################################

# schedules = CSV.read("lyapunov/lorenz_schedules_cache.csv", DataFrame)
# schedules[!, :λ1] .= .0; schedules[!, :λ2] .= .0; schedules[!, :λ3] .= .0;
# dt = 1e-4
# vrbl = [:dx, :dy, :dz], [:x, :y, :z]
# J_(x,y,z,ρ) = [
#      -10   10  0
#       ρ-z -1  -x
#         y  x  -8/3
# ]
# @showprogress @threads for dr = eachrow(schedules)
#         # filename = "bifurcation/hrnm/$(lpad(dr.idx, 5, '0')).csv"
#         filename = "lyapunov/lorenz/$(lpad(dr.idx, 5, '0')).csv"
#         # data = CSV.read(filename, DataFrame)
#         data = factory_lorenz(DataFrame, dr.ρ; tspan = [0, 2000.], dt = 1e-4)

#         # data = factory_hrnm(DataFrame, dr.bp; tspan = [0, 1500], dt = 1e-4); data = data[1000(nrow(data) ÷ 1500):end , :]
#         # data = factory_hrnm(DataFrame, dr.bp, ic = [dr.t, dr.x, dr.y, dr.z], tspan = [0, 1000]; dt = 1e-4)
#         # add_subsystem!(data, vrbl, cnfg; θ1, θ2, θ3, min_rank)
#         # CSV.write(filename, data)
#         # CSV.write(replace(filename, "bifurcation/hrnm" => "lyapunov/hrnm_traj"), data)

#         λ = lyapunov_exponent(data[:, last(vrbl)], J_, dr.ρ)

#         # data = data[(nrow(data) ÷ 2:end), :]
#         # idx_sampled = abs.(diff(data.dz)) .> 0.1
#         # sampledx = data[Not(1), :x][idx_sampled]
#         # append!(hrzn, fill(dr.bp, length(sampledx)))
#         # append!(vrtc, sampledx)
#         dr[[:λ1, :λ2, :λ3]] .= λ
# end
# # plot(schedules.λ1)
# using Plots
# plot(legend = :none)
# plot!(schedules.λ1)
# plot!(schedules.λ2)
# plot!(schedules.λ3)
# png("lyapunov/lorenz_laypunov.png")
# CSV.write("lyapunov/lorenz.csv", schedules, bom = true)
# temp = CSV.read("lyapunov/lorenz.csv", DataFrame)
# plot(legend = :none, xlabel = "ρ", ylabel = "λ")
# plot!(temp.ρ, temp.λ1)
# plot!(temp.ρ, temp.λ2)
# plot!(temp.ρ, temp.λ3)
# png("lyapunov/lorenz_laypunov2.png")

# ##########################################################################
# #                                                                        #
# #                           DC-DC buck converter                         #
# #                                                                        #
# ##########################################################################
# schedules = CSV.read("bifurcation/buck_schedules.csv", DataFrame)
# schedules[!, :λ1] .= .0; schedules[!, :λ2] .= .0;
# vrbl = [:dV, :dI], [:V, :I]
# cnfg = (; N = 1)
# dt = 1e-7; θ1 = 1e+1; θ2 = 1e+0; θ3 = 1e+0; min_rank = 2;
# function J_(V, I, E)
#     _R = 22 # 22
#     _L = 20_m # 20m
#     _C = 47_μ # 22μ
#     return [ -1/(_R*_C) (1/_C)
#                 -(1/_L)      0 ]
# end
# # _data = data[:, last(vrbl)]
# @showprogress @threads for dr = eachrow(schedules)
#     # filename = "bifurcation/buck/$(lpad(dr.idx, 5, '0')).csv"
#     # data = CSV.read(filename, DataFrame)
#     data = factory_buck(DataFrame, dr.E, tspan = [0, .1])
#     # add_subsystem!(data, vrbl, cnfg; θ1, θ2, θ3, min_rank); # 30 sec
#     # CSV.write(filename, data)

#     λ = lyapunov_exponent(data[:, last(vrbl)], J_, dr.E, T = .1)
#     dr[[:λ1, :λ2]] .= λ
#     # idx_sampled = diff(data.Vr) .< 0
#     # sampledV = data[Not(1), :V][idx_sampled]
#     # append!(hrzn, fill(dr.E, length(sampledV)))
#     # append!(vrtc, sampledV)
# end
# # plot(data.I[1:100:100000])
# CSV.write("lyapunov/buck.csv", schedules, bom = true)

# ##########################################################################
# #                                                                        #
# #                   Murali-Lakshmanan-Chua Circuit                       #
# #                                                                        #
# ##########################################################################
# schedules = CSV.read("bifurcation/mlcc_schedules.csv", DataFrame)[1:50:end, :]
# schedules[!, :λ1] .= .0; schedules[!, :λ2] .= .0; schedules[!, :λ3] .= .0; schedules[!, :λ4] .= .0;
# vrbl = [:dx1, :dx2, :dx3, :dx4], [:x1, :x2, :x3, :x4]
# dt = 1e-2; tend = 2000;
# function J_(x1, x2, x3, x4, β)
#     L = 19.1_m  ; Ga1 = -0.0009302325
#     C = 12.5_n  ; Ga2 = -0.000240577
#     ν = 8300    ; F = 0.3535533654213462

#     G = √(C/(L*β))
#     R = 1/G
#     a₁ = Ga1*R
#     a₂ = Ga2*R
#     f = F*β
#     ω = 2π*ν*C/G
#     return [ 0                           1  0              0
#              0 ifelse(abs(x1) ≤ 1, a₁, a₂)  1              0
#              0                          -β -β -f*ω*cos(ω*x4)
#              0                           0  0              0 ]
# end

# hrzn, vrtc = Dict(), Dict()
# @showprogress for dr = eachrow(schedules)
#         data = factory_mlcc(DataFrame, dr.bp; dt, tspan = [0, tend])

#         λ = lyapunov_exponent(data[:, last(vrbl)], J_, dr.bp, T = tend)
#         dr[[:λ1, :λ2, :λ3, :λ4]] .= λ

#         data = data[data.x4 .> (tend-200), :]
#         idx_sampled = abs.(diff(data.dx2)) .> 0.01
#         sampled = data[Not(1), :x3][idx_sampled]
#         push!(hrzn, dr.idx => fill(dr.bp, length(sampled)))
#         push!(vrtc, dr.idx => sampled)
# end
# bfcn = DataFrame(hrzn = vcat(values(hrzn)...), vrtc = vcat(values(vrtc)...))
# # CSV.write("lyapunov/mlcc_lyapunov.csv", schedules, bom = true)
# # CSV.write("lyapunov/mlcc_bifurcation.csv", bfcn, bom = true)
# scatter(bfcn.hrzn, bfcn.vrtc, legend = false, alpha = .5, ms = .1, xlabel = "β", ylabel = "x3")
# # png("lyapunov/mlcc_bifurcation.png")
# plot(legend = :none); plot!(schedules.λ1); plot!(schedules.λ1); plot!(schedules.λ1)
# # png("lyapunov/mlcc_lyapunov.png")
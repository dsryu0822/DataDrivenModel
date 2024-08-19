include("../core/header.jl")
using Plots

function lyapunov_exponent(_data::DataFrame, J_::AbstractMatrix, bf_param;
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

#         # data = factory_hrnm(DataFrame, dr.f; tspan = [0, 1500], dt = 1e-4); data = data[1000(nrow(data) ÷ 1500):end , :]
#         # data = factory_hrnm(DataFrame, dr.f, ic = [dr.t, dr.x, dr.y, dr.z], tspan = [0, 1000]; dt = 1e-4)
#         # add_subsystem!(data, vrbl, cnfg; θ1, θ2, θ3, min_rank)
#         # CSV.write(filename, data)
#         # CSV.write(replace(filename, "bifurcation/hrnm" => "lyapunov/hrnm_traj"), data)

#         λ = lyapunov_exponent(data[:, last(vrbl)], J_, dr.ρ)

#         # data = data[(nrow(data) ÷ 2:end), :]
#         # idx_sampled = abs.(diff(data.dz)) .> 0.1
#         # sampledx = data[Not(1), :x][idx_sampled]
#         # append!(hrzn, fill(dr.f, length(sampledx)))
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
# #                            Soft impact model                           #
# #                                                                        #
# ##########################################################################
# schedules = CSV.read("bifurcation/soft_schedules.csv", DataFrame)
# schedules[!, :λ1] .= .0; schedules[!, :λ2] .= .0; schedules[!, :λ3] .= .0;
# vrbl = [:dt, :du, :dv], [:t, :u, :v]
# cnfg = (; f_ = [cospi, sign], λ = 1e-2)
# dt = 1e-5; θ1 = 1e-8; θ2 = 1e-12; θ3 = 1e-5; min_rank = 21;
# function J_(t, u, v, d)
#     return [   0                                0                                 0
#                0                                0                                 1
#      -π*sinpi(t) ifelse(abs(u) ≥ d/2, -160000, 0) ifelse(abs(u) ≥ d/2, -172.363, 0) ]
# end
# @showprogress @threads for dr = eachrow(schedules)
#     # filename = "bifurcation/soft/$(lpad(dr.idx, 5, '0')).csv"
#     # data = CSV.read(filename, DataFrame)
#     data = factory_soft(DataFrame, dr.d, tspan = [0, 1000])
#     # add_subsystem!(data, vrbl, cnfg; θ1, θ2, θ3, min_rank); # 30 sec
#     # CSV.write(filename, data)

#     λ = lyapunov_exponent(data[:, last(vrbl)], J_, dr.d)
#     # idx_sampled = diff(abs.(data.u) .> (dr.d/2)) .> 0
#     # sampledv = data[Not(1), :v][idx_sampled]
#     # append!(hrzn, fill(dr.d, length(sampledv)))
#     # append!(vrtc, sampledv)
# end
# CSV.write("lyapunov/!linux soft t = [0, 1000].csv", schedules, bom = true)

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
# #                           Hindmarsh-Rose model                         #
# #                                                                        #
# ##########################################################################
# schedules = CSV.read("bifurcation/hrnm_schedules.csv", DataFrame)
# # schedules = schedules[.!isfile.(["lyapunov/hrnm_traj/$(lpad(idx, 5, '0')).csv" for idx in 1:nrow(schedules)]), :]
# # schedules = CSV.read("lyapunov/hrnm_schedules_cache.csv", DataFrame)[1:10:401, :]
# schedules[!, :λ1] .= .0; schedules[!, :λ2] .= .0; schedules[!, :λ3] .= .0; schedules[!, :λ4] .= .0;
# vrbl = [:dt, :dx, :dy, :dz], [:t, :x, :y, :z]
# cnfg = (; N = 3, f_ = [cos])
# dt = 1e-3; θ1 = 1e-2; θ2 = 1e-27; θ3 = 1e-1; min_rank = 32;
# function J_(t, x, y, z, _f)
#     return [                0                             0   0    0
#              -_ω*_f*sin(_ω*t) (-3*_a*(x^2) + 2*_b*x + _k*z)   1 _k*x
#                             0                       -2*_d*x  -1    0
#                             0                            _β   0  -_α]
# end

# bfcn = DataFrame(hrzn = [], vrtc = [])
# @showprogress @threads for dr = eachrow(schedules)[136]
#         # filename = "bifurcation/hrnm/$(lpad(dr.idx, 5, '0')).csv"
#         # filename = "lyapunov/hrnm_traj/$(lpad(dr.idx, 5, '0')).csv"
#         # data = CSV.read(filename, DataFrame)
#         data = factory_hrnm(DataFrame, dr.f; tspan = [0, 10000], dt)
#         # add_subsystem!(data, vrbl, cnfg; θ1, θ2, θ3, min_rank)
#         # CSV.write(replace(filename, "bifurcation/hrnm" => "lyapunov/hrnm_traj"), data)
#         # CSV.write(filename, data)

#         λ = lyapunov_exponent(data[:, last(vrbl)], J_, dr.f)
#         dr[[:λ1, :λ2, :λ3, :λ4]] .= λ

#         # data = data[data.t .> 3000, :]
#         # idx_sampled = abs.(diff(data.dz)) .> 0.1
#         # sampledx = data[Not(1), :x][idx_sampled]
#         # hrzn, vrtc = fill(dr.f, length(sampledx)), sampledx
#         # append!(bfcn, DataFrame(; hrzn, vrtc))
# end
# scatter(bfcn.hrzn, bfcn.vrtc, legend = false, alpha = .5, ms = .1, xlabel = "f", ylabel = "x")
# png("lyapunov/!linux hrnm_bifurcation 1e-4 t = [0, 4000].csv")
# CSV.write("lyapunov/!linux hrnm_bifurcation 1e-4 t = [0, 4000].csv", bfcn, bom = true)
# CSV.write("lyapunov/!linux hrnm_lyapunov 1e-4 t = [0, 10000].csv", schedules, bom = true)
# CSV.write("lyapunov/hrnm_bifurcation_test.csv", DataFrame(; hrzn, vrtc), bom = true)
# scatter(hrzn, vrtc, legend = false, alpha = .5, ms = .1, xlabel = "f", ylabel = "x")
# CSV.write("lyapunov/hrnm_bifurcation_test.csv", bfcn, bom = true)
# CSV.write("lyapunov/hrnm 1e-3 test.csv", schedules, bom = true)
# plot(data.z[1:100:end], color = data.subsystem[1:100:end])

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


##########################################################################
#                                                                        #
#                             Gear system                                #
#                                                                        #
##########################################################################
# data = factory_gear(DataFrame, -.06; dt = 1e-2, tspan = [500, 5000])
# plot(data.x, data.v)
# idx_sampled = diff([0; mod.(data.θ, 2π)]) .< 0
# data = data[idx_sampled, :]
# scatter(data.x, data.v, ms = 0.5, legend = :none)
schedules = CSV.read("bifurcation/gear_schedules.csv", DataFrame)
schedules[!, :λ1] .= .0; schedules[!, :λ2] .= .0; schedules[!, :λ3] .= .0; schedules[!, :λ4] .= .0;
vrbl = [:dx, :dv, :dΩ, :dθ], [:x, :v, :Ω, :θ]
dt = 1e-2; tend = 10000;
function J_(x, v, Ω, θ, Fe)
    dfdx = ifelse(abs(x) > 1, 1, 0)
    return [ 0 1 0 0
             -(1 + k1*cos(Ω))*dfdx -2ζ (-Fe*H*sin(Ω) + k1*sin(Ω)*dfdx) -Fe*H*sin(θ)
             0 0 0 0
             0 0 0 0 ]
end

bfcn = DataFrame(hrzn = [], vrtc = [])
@showprogress for dr = eachrow(schedules)
    data = factory_gear(DataFrame, dr.Fe; tspan = [0, tend])

    λ = lyapunov_exponent(data[:, last(vrbl)], J_, dr.Fe, T = tend)
    dr[[:λ1, :λ2, :λ3, :λ4]] .= λ

    data = data[data.Ω .> 5000, :]
    idx_sampled = diff([0; mod.(data.Ω, 2π)]) .< 0
    sampledx = data.v[idx_sampled]
    hrzn, vrtc = fill(dr.Fe, length(sampledx)), sampledx
    append!(bfcn, DataFrame(; hrzn, vrtc))
end
CSV.write("lyapunov/gear_bifurcation.csv", bfcn, bom = true)
CSV.write("lyapunov/gear_lyapunov.csv", schedules, bom = true)

scatter(bfcn.hrzn, bfcn.vrtc, legend = false, alpha = .5, ms = .1, xlabel = "Fe", ylabel = "v")
png("lyapunov/gear_bifurcation.png")

plot(legend = :none)
plot!(schedules.Fe, schedules.λ1)
plot!(schedules.Fe, schedules.λ2)
plot!(schedules.Fe, schedules.λ3)
plot!(schedules.Fe, schedules.λ4)
png("lyapunov/gear_lyapunov.png")

println("Done!")
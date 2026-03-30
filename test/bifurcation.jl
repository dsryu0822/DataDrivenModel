include("../core/header.jl")

# function factory_lorenz96(A::Number; ic = [0,0,1,1,1,1,1], tspan = [0., 10.], dt = 1e-4)
#     function sys(v::AbstractVector)
#         T, x1, x2, x3, x4, x5, x6 = v

#         dT = 1
#         dx1 = x6*(x2 - x5) - x1
#         dx2 = x1*(x3 - x6) - x2
#         dx3 = x2*(x4 - x1) - x3
#         dx4 = x3*(x5 - x2) - x4
#         dx5 = x4*(x6 - x3) - x5
#         dx6 = x5*(x1 - x4) - x6
#         return [dT; [dx1, dx2, dx3, dx4, dx5, dx6] .+ (A*sin(2T) + 2)]
#     end

        
#     t_ = first(tspan):dt:last(tspan)
#     len_t_ = length(t_)
    
#     t, tk = .0, 0
#     v = ic; DIM = length(v)
#     traj = zeros(len_t_+2, 2DIM)
#     while tk ≤ len_t_
#         t,_,_ = v
#         v, dv = RK4(sys, v, dt)

#         if t ≥ first(t_)
#             tk += 1
#             traj[tk+1,         1:DIM ] =  v
#             traj[tk  , DIM .+ (1:DIM)] = dv
#         end
#     end
#     return traj[2:(end-2), :]
# end
# factory_lorenz96(T::Type, args...; kargs...) =
# DataFrame(factory_lorenz96(args...; kargs...), ["t", "x1", "x2", "x3", "x4", "x5", "x6", "dt", "dx1", "dx2", "dx3", "dx4", "dx5", "dx6"])

# A_ = 1:0.01:4
# @showprogress @threads for k in eachindex(A_)
#     A = A_[k]
#     traj = factory_lorenz96(DataFrame, A, tspan = [0., 2500], dt = 1e-3)
#     _traj = traj[150001:end, :]
#     CSV.write("G:/lorenz96/lorenz96_$(lpad(k, 4, '0')).csv", _traj)
# end

# function arglocalmax(x)
#     bits = circshift(x, 1) .< x .> circshift(x, -1)
#     bits[1] = false
#     bits[end] = false
#     return findall(bits)
# end
# scatter(repeat([3.2], length(_traj.x1[almx1])), _traj.x1[almx1])

# plt_bfcn = plot(legend = :none)
# @showprogress for k in eachindex(A_)
#     _traj = CSV.read("G:/lorenz96/lorenz96_$(lpad(k, 4, '0')).csv", DataFrame)
#     x1 = _traj.x1[_traj.t .> 1000]
#     almx1 = arglocalmax(x1)
#     scatter!(plt_bfcn, repeat([A_[k]], length(x1[almx1])), x1[almx1], color = :black, ms = 1, ma = 0.05)
# end
# plt_bfcn
# png(plt_bfcn, "bifurcation_plot.png")

# A_[121]
# A_[161]
traj0 = CSV.read("G:/lorenz96/lorenz96_0121.csv", DataFrame)
traj1 = CSV.read("G:/lorenz96/lorenz96_0161.csv", DataFrame)

vrbl = ([:dt, :dx1, :dx2, :dx3, :dx4, :dx5, :dx6], [:t, :x1, :x2, :x3, :x4, :x5, :x6])
cnfg = cook(last(vrbl), poly = 0:2, trig = 0:2, format = sin)
@time f0 = SINDy(traj0, vrbl, cnfg, λ = 0.01)
@time f1 = SINDy(traj1, vrbl, cnfg, λ = 0.01)
@info "SINDy done!"

function syntheticSINDy(Ξ, sysms::Tuple, recipe::AbstractDataFrame)
    Ysyms, Xsyms = sysms
    bit_sparse = all.(map(x -> iszero.(x), eachrow(Ξ)))
    recipeF = recipe[.!bit_sparse, :]
    _Ξ = Ξ[.!bit_sparse, :]
    mse = 0
    aic = -Inf
    r2 = 1
    return STLSQresult(recipe, recipeF, Ξ, _Ξ, mse, aic, r2, Ysyms, Xsyms)
end

# f0.matrix
# α = 0.1
# α*2.2 + (1-α)*2.6

α_ = range(-3.5, 4, length(1:0.01:4))
@showprogress @threads for k in eachindex(α_)
    α = α_[k]
    g = syntheticSINDy(α*f0.matrix + (1-α)*f1.matrix, vrbl, cnfg)
    @time traj = solve(g, [0,0,1,1,1,1,1], 0:1e-3:25)
    # _traj = DataFrame(traj[1500001:end, :], ["t", "x1", "x2", "x3", "x4", "x5", "x6"])
    _traj = DataFrame(traj, ["t", "x1", "x2", "x3", "x4", "x5", "x6"])
    CSV.write("../lorenz96SINDy/lorenz96SINDy_$(lpad(k, 4, '0')).csv", _traj)
end

include("../core/header.jl")

function factory_tcg(α::Number; ic = [0, .1, 0, 0], tspan = (0., 1.), dt = 1e-4)
    ξ1 = ξ2 = .1
    γ = .1
    δ = .1
    Ω1 = .7
    F1 = .2
    function sys(txyz::AbstractVector, nonsmooth::Real)
        t,x,y,z = txyz

        ṫ = 1
        ẋ = y
        ẏ = -(ξ1*ẋ) - (x*(1 - (inv(sqrt(x^2 + α^2))))) + ((z^2)*x) + (F1*cos(Ω1*t))
        ż = -(ξ2*z) - γ*nonsmooth + δ
        return [ṫ, ẋ, ẏ, ż]
    end

    t_ = first(tspan):dt:last(tspan)
    len_t_ = length(t_)
    
    t, tk = .0, 0
    v = ic; DIM = length(v)
    traj = zeros(len_t_+2, 2DIM)
    while tk ≤ len_t_
        t,x,y,z = v
        nonsmooth = abs(x)
        v, dv = RK4(sys, v, dt, nonsmooth)

        if t ≥ first(t_)
            tk += 1
            traj[tk+1,         1:DIM ] =  v
            traj[tk  , DIM .+ (1:DIM)] = dv
        end
    end
    return traj[2:(end-2), :]
end
factory_tcg(T::Type, args...; kargs...) = 
DataFrame(factory_tcg(args...; kargs...), ["t", "x", "y", "z", "dt", "dx", "dy", "dz"])

@time data = factory_tcg(DataFrame, 0.5, ic = [0, 1, 0, 0], tspan = (0., 100.), dt = 1e-2)[:, Not(:dt)]
data[!, :label] = ((data.x .> 0) .+ 1)
@time CSV.write("output.csv", data)
_data = data[1:(nrow(data) ÷ 10000):end, :]; _data = _data[5(nrow(_data)÷10):end, :]
scatter(diff(_data.dz))
plot(plot(_data.y, _data.x, _data.z, xlabel = "y", ylab = "x", zlab = "z", color = :black),
     plot(_data.x, _data.y, color = 2),
     plot(_data.x, _data.z, color = 3),
     plot(_data.y, _data.z, color = 1), size = (800, 800), legend = :none)
plot(_data.t, _data.z)

function localminima(x)
    _x = circshift(x, 1)
    x_ = circshift(x, -1)
    return x[(x .< _x) .& (x .< x_)][2:end-1]
end

vrtl = []
@showprogress for α = 0:0.001:0.5
    data = factory_tcg(DataFrame, α, ic = [0, 1, 0, 0], tspan = (0., 1000.), dt = 1e-1)[:, Not(:dt)]
    _data = data[1:(nrow(data) ÷ 10000):end, :]; _data = _data[7(nrow(_data)÷10):end, :]
    push!(vrtl, localminima(_data.z))
end
bfcn = plot()
for (k, α) in enumerate(0:0.001:0.5)
    try
        scatter!(bfcn, fill(α, length(vrtl[k])), vrtl[k], color = :black, ms = 1, label = "")
    catch e
        @warn "Error plotting for α = $α: $e"
    end
end
bfcn
vrtl[1]
scatter(localminima(data.x))
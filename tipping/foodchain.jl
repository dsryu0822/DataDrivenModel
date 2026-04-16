include("../core/header.jl")
include("tippingutils.jl")

import DifferentialEquations as DE
import Sundials
import DiffEqBase


function define(T::Type, sindy::STLSQresult)
    header = join(setdiff(sindy.rname, sindy.lname), ", ") * " = u; " * join(sindy.lname, ", ") * " = du" 
    if sindy.method == "SINDyPI"
        header = "function f(out, du, u, p, t)\n" * header
    end
    factors = sindy.recipe.term[sindy.recipe.func .== identity]
    term = join.([factors[idx_term] for idx_term in sindy.recipe.vecv], "*")
    # du = sindy.lname .* " = "
    du = ["out[$j]" for j in eachindex(sindy.lname)] .* " = "
    for j in eachindex(du)
        for i in eachindex(term)
            if sindy.matrix[i, j] |> !iszero
                du[j] *= " + "*string(round(sindy.matrix[i, j], sigdigits = 5))*term[i]
            end
        end
    end
    du = du .* " - " .* sindy.lname
    body = replace(join(du, "\n"), "=  + " => "= ", "+ -" => "- ")
    footer = "end"
    function_string = join([header, body, footer], "\n")
    if T == Function
        ex = Meta.parse(function_string)
        return eval(ex)
    else
        if T != String
            @warn "T should be either Function or String. Returning String by default."
        end
        return function_string
    end
end
define(f) = define(String, f)

"""
    termprod(config::AbstractDataFrame, term::AbstractString, num::Integer)

    - `config`: a configuration DataFrame.
    - `term`: a string representing the new term to be added to the configuration.
    - `num`: an integer representing the index of the new term in the configuration.
"""
function termprod(config::AbstractDataFrame, term::AbstractString, num::Integer)
    _config = deepcopy(config)
    _config.term  = replace.("$term⋅" .* _config.term, "⋅1" => "")
    for i in eachindex(_config.func)
        push!(_config.vecv[i], num)
        if _config.func[i] == constant
            _config.func[i] = identity
            _config.vecv[i] = [num]
        elseif _config.func[i] == identity
            _config.func[i] = prod
            _config.funh[i] = eachrow
        end
    end
    return _config
end
"""
    cookPI(variable; kargs...)

    - `variable`: a tuple of two vectors of strings.
    - `kargs`: see `cook`.

"""
function cookPI(variable; kargs...)
    vnames = last(variable)
    dvnames = first(variable)

    config = cook(vnames; kargs...)
    masking = [zeros(Bool, nrow(config), length(dvnames))]
    config_ = [config]
    for l in eachindex(dvnames)
        newconfig = termprod(config, dvnames[l], l+length(vnames))
        masking_tmp = ones(Bool, nrow(config), length(dvnames)); masking_tmp[2:end, l] .= false
        push!(masking, masking_tmp)
        push!(config_, newconfig)
    end
    _config = [config_...;]
    _config.index = 1:nrow(_config)
    return _config, [masking...;]
end
function SINDyPI(df, variable, config; kargs...)
    _variable = (first(variable), [reverse(variable)...;])
    _config, mask = config
    return SINDy(df, _variable, _config; mask = mask, method = "SINDyPI", kargs...)
end

frac(x, y) = x ./ y
function factory_foodchain(a0::Number; ic = [0., 1, 1, 1], tspan = 0:1e-2:10)
      b0, w0, d0, w1, d1,  w2, d2, a1, w3, d3,  c3 = (
    0.05,  1, 10,  2, 10, 1.5, 10,  1,  2, 20, 0.7 )
    function sys(v::AbstractVector)
        _,v1,v2,v3 = v

        dt = 1
        dv1 = a0*v1 - b0*v1^2 - w0*frac(v1*v2, d0 + v1)
        dv2 = w1*frac(v1*v2, d1 + v1) - w2*frac(v2*v3, d2 + v2) - a1*v2
        dv3 = w3*frac(v2*v3, d3 + v2) - c3*v3
        return [dt, dv1, dv2, dv3]
    end
        
    len_t_ = length(tspan); h = tspan.step.hi
    v = ic; DIM = length(v)
    traj = zeros(2+len_t_, 2DIM); traj[1, 1:DIM] = v

    tk = 0
    while tk ≤ len_t_
        t,_,_,_ = v
        v, dv = RK4(sys, v, h)

        if t ≥ first(tspan)
            tk += 1
            traj[tk+1,         1:DIM ] =  v
            traj[tk  , DIM .+ (1:DIM)] = dv
        end
    end
    return traj[1:(end-2), :]
end
factory_foodchain(T::Type, args...; kargs...) =
DataFrame(factory_foodchain(args...; kargs...), ["t", "v1", "v2", "v3", "dt", "dv1", "dv2", "dv3"])[:, Not(:dt)]

# hrzn = []
# vrcl = []
# for a0 = 0.15:0.001:2.0
#     data = factory_foodchain(DataFrame, a0, tspan = 0:1e-2:1000)[50001:end, :]
#     bits = arglmax(data.v1)
#     push!(hrzn, [a0 for _ in bits])
#     push!(vrcl, data.v1[bits])
# end
# scatter([hrzn...;], [vrcl...;], ms = 1, msw = 0, color = :black, alpha = 0.5, legend = :none, xlims = [1.5, 2.0], ylims = [25, 40])
# png("temp")

traj0 = factory_foodchain(DataFrame, 1.6, tspan = 0:1e-2:1000)[50000:end, :]
traj1 = factory_foodchain(DataFrame, 1.9, tspan = 0:1e-2:1000)[50000:end, :]
traj2 = factory_foodchain(DataFrame, 2.0, tspan = 0:1e-2:1000)[50000:end, :]
plot(
    plot(traj0.v1, traj0.v2, traj0.v3, color = :black),
    plot(traj1.v1, traj1.v2, traj1.v3, color = :black),
    plot(traj2.v1, traj2.v2, traj2.v3, color = :black),
    legend = :none, layout = (1, 3), size = [1200, 400]
)

vrbl = reverse(half(names(traj0)[Not(1)]))
cnfg = cookPI(vrbl; poly = 0:3)

f0 = SINDyPI(traj0, vrbl, cnfg; λ = 1e-7); f0 |> print
f1 = SINDyPI(traj1, vrbl, cnfg; λ = 1e-7); f1 |> print
f2 = SINDyPI(traj2, vrbl, cnfg; λ = 1e-7); f2 |> print

tspan = (0, 1000)
deargs = (; reltol = 1e-6, initializealg = DiffEqBase.BrownFullBasicInit(), saveat = 0:1e-2:1000)
prob0 = DE.DAEProblem(define(Function, f0), collect(traj0[1, 5:end]), collect(traj0[1, 2:4]), tspan, differential_vars = ones(Bool, 3))
sol0 = DE.solve(prob0, Sundials.IDA(); deargs...)

prob1 = DE.DAEProblem(define(Function, f1), collect(traj1[1, 5:end]), collect(traj1[1, 2:4]), tspan, differential_vars = ones(Bool, 3))
sol1 = DE.solve(prob1, Sundials.IDA(); deargs...)

prob2 = DE.DAEProblem(define(Function, f2), collect(traj2[1, 5:end]), collect(traj2[1, 2:4]), tspan , differential_vars = ones(Bool, 3))
sol2 = DE.solve(prob2, Sundials.IDA(); deargs...)

plot(
    plot(eachrow(stack(sol0.u))..., ticks = [], color = :blue),
    plot(eachrow(stack(sol1.u))..., ticks = [], color = :blue),
    plot(eachrow(stack(sol2.u))..., ticks = [], color = :blue),
    legend = :none, layout = (1, 3), size = [1200, 400]
)

α_ = 1:-0.01:-4
hrzn = []
vrtl = []
@showprogress for k in eachindex(α_)
    α = α_[k]
    g = syntheticSINDy((1-α)*f1.matrix + α*f2.matrix, vrbl, cnfg[1])
    probg = DE.DAEProblem(define(Function, g), collect(traj0[1, 5:end]), collect(traj0[1, 2:4]), (0, 1000), differential_vars = [true, true, true])
    solg = DE.solve(probg, Sundials.IDA(); deargs...)
    traj = stack(solg.u)'[solg.t .≥ 500, :]
    points = arglmax(traj[:, 1])
    push!(vrtl, traj[points, 1])
    push!(hrzn, repeat([α], length(points)))
end
scatter(hrzn, vrtl, ms = 1, msw = 0, color = :blue, alpha = 0.5, legend = :none, xlims = [-4, 1], ylims = [25, 40])
png("temp")



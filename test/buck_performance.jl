include("buck_.jl")

function recover(ip, SINDy, Classifier)
    v_ = [ip]
    d_ = [foo(v_[end])]
    dt = 10^(-7)
    for (tk, t) in ProgressBar(enumerate(dt:dt:0.001))
        v_[end][end] = Classifier([v_[end][1:2]; DATA.Vr[tk]]) |> argmax |> Float64
        v, d = RK4(SINDy, v_[end], dt)
        push!(v_, v)
        push!(d_, d)
    end
    prdt = stack(v_)
    dprdt = stack(d_)
    return prdt
end
function summaryplot(DATA, TRAJ, t_)
    a1 = plot(DATA.V[_t], DATA.I[_t], color = :black, alpha = 0.5, label = "Data")
    a2 = plot(a1, TRAJ.V[_t], TRAJ.I[_t], color = :blue, ls = :dash, label = "Predicted")
    b1 = plot(DATA.t[_t], DATA.V[_t], color = :black, alpha = 0.5, label = "Data")
    b2 = plot(b1, DATA.t[_t], TRAJ.V[_t], color = :blue, ls = :dash, label = "Predicted")
    c1 = plot(DATA.t[_t], DATA.I[_t], color = :black, alpha = 0.5, label = "Data")
    c2 = plot(c1, DATA.t[_t], TRAJ.I[_t], color = :blue, ls = :dash, label = "Predicted")
    
    bc = plot(b2, c2, layout = (2,1), legend = :best)
    abc = plot(a2, bc, size = (1600, 900)); png(abc, "abc.png")
    return abc
end
breakaway(DATA, TRAJ, _t; ε = 10^(-3)) = findfirst(abs.(DATA.V[_t] - TRAJ.V[_t]) .> ε)

ANN, loss, acry, epch = load_ANN("C:/Temp/SSS_000023.jld2")
traj = recover(ip, subeqs, ANN)
_t = axes(traj, 2)
TRAJ = DataFrame(traj', [:V, :I, :now])

summaryplot(DATA, TRAJ, _t)
breakaway(DATA, TRAJ, _t)

acry_ = []
brkw_ = []
for file in readdir("C:/Temp")[Not(1)]
    @info file
    ANN, loss, acry, epch = load_ANN("C:/Temp/" * file)
    traj = recover(ip, subeqs, ANN)
    _t = axes(traj, 2)

    push!(acry_, acry)
    push!(brkw_, breakaway(DATA, TRAJ, _t))
end
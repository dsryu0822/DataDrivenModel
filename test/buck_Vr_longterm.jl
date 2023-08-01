using ProgressBars
packages = [:DataFrames, :CSV, :LinearAlgebra, :Plots, :Flux, :Clustering, :JLD2, :LaTeXStrings, :CUDA]
@time for package in ProgressBar(packages)
    @eval using $(package)
end
println(join(packages, ", "), " loaded!")

include("../src/DDM.jl")
include("../src/ML.jl")
include("../src/nonsmooth.jl")
include("../src/ODEdata.jl")
default(size = (600,600), color = :black, legend = :none)

## Data Load
DATA = CSV.read("G:/buck/buck_000006.csv", DataFrame)
Y = DATA[end-100000:end,[:dV, :dI]] |> Matrix # .|> Float32
X = DATA[end-100000:end,[ :V,  :I]] |> Matrix # .|> Float32
# XY = [X Y]

## Clustering
@showtime dbs = dbscan(col_normalize(Y)', 0.001);
nsubsys = length(dbs.clusters); println(nsubsys, " clusters found!")
# plot(X[:,1], X[:,2], color = dbs.assignments, alpha = 0.5)
# scatter(Y[:,1], Y[:,2], color = dbs.assignments, alpha = 0.5, msw = 0)
bit_ = [dbs.assignments .== k for k in 1:nsubsys]
Θ_ = [poly_basis(X[s,:], 2)     for s in bit_]
Y_ = [Y[vcat(findall(s)...),:] for s in bit_]
Ξ_ = [STLSQ(Θ_[s], Y_[s], 0.01) for s in 1:nsubsys]

function foo(v)
    v = deepcopy(v)
    s = Int64(pop!(v))
    Θv = poly_basis(v, 2)'
    return vec([(Θv*Ξ_[s]) 0])
end


# a1 = plot(_DATA.V, _DATA.I, alpha = 0.1, legend = :best, label = "Trajectory", xlabel = L"V", ylabel = L"I")
# scatter!(a1,
#     _DATA.V[Not(1)][.!iszero.(diff(dbs.assignments))]
#   , _DATA.I[Not(1)][.!iszero.(diff(dbs.assignments))]
#   , color = :red, shape = :+, label = "Jumping points")


subsystem = (DATA.dI .> 0) .+ 1
# _DATA = DATA[Not(end), :]
# _DATA[!, :now] = subsystem[Not(end)]
# _DATA[!, :next] = subsystem[Not(1)]

_DATA = DATA
_DATA[!, :now] = subsystem
__DATA = _DATA[1:1:end,:]
__next = __DATA.now[Not(1)]
__DATA = __DATA[Not(end), :]
__DATA[!, :next] = __next
_DATA = __DATA

# _DATA[_DATA.now .!= _DATA.next, :]
# idx_change = findall(_DATA.now .!= _DATA.next)
# # _DATA[51950:51970, :]
# # _DATA[58740:58760, :]
# idx_data = [1:100:nrow(_DATA);
#             idx_change;
#             # idx_change .+ 1;
#             # idx_change .- 1;
#             # idx_change .+ 2;
#             # idx_change .- 2;
#             # idx_change .+ 3;
#             # idx_change .- 3;
#             # idx_change .+ 4;
#             # idx_change .- 4;
#             # idx_change .+ 5;
#             # idx_change .- 5;
#             5000001 .- idx_change] |> sort |> unique
# idx_data = idx_data[1 .≤ idx_data .≤ 5000000]
# _DATA = _DATA[idx_data, :]
# _DATA = _DATA[end-50000:100:end, :]

## Classification
target = deepcopy(_DATA.now)
data = Float32.(Matrix(select(_DATA, [:V, :I, :Vr]))')
subs = Flux.onehotbatch(target, 1:nsubsys)
trng = Flux.DataLoader((data,subs), shuffle = true, batchsize = 10000)

p = size(data, 1)
SSSf = Chain( # SubSystemSelector
    Flux.Scale(p),
    Fourier(100),
    Dense(p+1 => 50, relu),
    Dense(50 => 50, relu),
    Dense(50 => 50, relu),
    Dense(50 => nsubsys),
    softmax
)
Loss(x,y) = Flux.crossentropy(SSSf(x), y)
loss_ = [Loss(data, subs)]
acry_ = [sum(target .== argmax.(eachcol(SSSf(data)))) / size(data, 2)]
optimizer = ADAM()

norm_variables = Float32.(norm.(eachrow(Matrix(data)))); norm_variables[end] = 1.0
SSSf[1].scale ./= norm_variables

if !isfile("data/SSSf.jld2")
    print("data/SSSf.jld2 not exists, ")
    
    ps = Flux.params(SSSf)
    println("Training... First Loss: ", loss_[1])
    @time for epch in 1:100_000
        if epch < 10
            @time Flux.train!(Loss, ps, trng, optimizer)
        else
            Flux.train!(Loss, ps, trng, optimizer)
        end
        loss = Loss(data, subs)
        acry = sum(target .== argmax.(eachcol(SSSf(data)))) / size(data, 2)
        if acry < acry_[end]
        # if loss < loss_[end]
            jldsave("C:/Temp/SSSf.jld2"; SSSf, loss, acry) # , loss_, acry_)
            println()
            print("epoch ", lpad(length(loss_), 5)
                    , ": loss = ", rpad(trunc(loss, digits = 6), 8)
                    , ", acry = ", rpad(trunc(100acry, digits = 4), 8)
                    , " saved!")
            push!(loss_, loss)
            push!(acry_, acry)
        else
            epch % 100 == 0 && print(", ", lpad(epch, 5))
            push!(loss_, loss_[end])
            push!(acry_, acry_[end])
        end
    end
    cp("C:/Temp/SSSf.jld2", "data/SSSf.jld2")
else
    @info "Loading SSSf..."
    fSSSf = jldopen("data/SSSf.jld2")
    SSSf = fSSSf["SSSf"]
    println("loss = ", fSSSf["loss"])
    println("acry = ", fSSSf["acry"])
    close(fSSS)
end

### Classifier testing
# a1_1 = scatter(a1, 
#     _DATA.V[Not(1, 2)][.!iszero.(diff(argmax.(eachcol(SSSf(data)))))]
#   , _DATA.I[Not(1, 2)][.!iszero.(diff(argmax.(eachcol(SSSf(data)))))]
#   , color = 2, shape = :+, label = "Detected Jumping points")
# png(a1_1, "a1_1.png")

## ODE recovery
dt = 10^(-7); tspan = 0:dt:0.01
zeros(5, length(tspan))
v_ = [Float64[data[1:2, 1]; argmax(subs[:,1])]]
d_ = [foo(v_[end])]
for (tk, t) in ProgressBar(enumerate(tspan))
    v, d = RK4(foo, v_[end], dt)
    push!(v_, v)
    push!(d_, d)
    # v_[end][end] = [v_[end][1:2]; _DATA.Vr[tk]]  |> SSSf |> vec |> argmax |> Float64
    v_[end][end] = [v_[end][1:2]; t]  |> SSSf |> vec |> argmax |> Float64
    # v_[end][end] = [v_[end][1:2]; d_[end][1:2]; t]  |> SSSf |> vec |> argmax |> Float64
end

### Trajectories
prdt = stack(v_)
a1_2 = plot(prdt[1,:], prdt[2,:], legend = :best, label = "Recovered system", xlabel = L"V", ylabel = L"I")
png(a1_2, "a1_2.png")

a1_2 = plot(prdt[1,:], prdt[2,:], legend = :best, label = "Recovered system", xlabel = L"V", ylabel = L"I")
v1 = plot(DATA.V[1:100001], legend = :best, color = :black, label = "V")
plot!(v1, prdt[1,:], ylabel = L"V", color = :blue, style = :dash, label = "predicted V")
i1 = plot(DATA.I[1:100001], ylabel = L"I", xlabel = L"t", legend = :none)
plot!(i1, prdt[2,:], color = :blue, style = :dash)
a3 = plot(
    v1, i1
    , layout = (2,1), xformatter = x -> x*dt, size = (800,600)
)
png(a3, "a3.png")

### Fourier check
plot(Flux.params(SSSf[2])[1], linetype = :bar)
plot(Flux.params(SSSf[2])[2], linetype = :bar)
[Flux.params(SSSf[2])[1] Flux.params(SSSf[2])[2]]
plot(DATA.t)

y = (1:50) .* DATA.t' .* 1000000exp(Flux.params(SSSf[2])[3][1])
plot(vec(sum(
    (Flux.params(SSSf[2])[1] .* cospi.(y)) +
    (Flux.params(SSSf[2])[2] .* sinpi.(y)), dims = 1)))

y = (1:50) .* DATA.t' .* 1000000exp(Flux.params(SSSf[2])[3][1])
plot(vec(sum(
        (randn(50) .* cospi.(y)) +
        (randn(50) .* sinpi.(y)), dims = 1)))

# ----------------------------------------------------------------------------------

v_ = [[XY[1, 1:2]; dbs.assignments[1]]]
d_ = [foo(v_[end])]
dt = 10^(-7)
for (tk, t) in ProgressBar(enumerate(dt:dt:0.015))
    v_[end][end] = v_[end][1] > DATA.Vr[tk] ? 1 : 2
    push!(v_, RK4(foo, v_[end], dt))
    push!(d_, foo(v_[end]))
end
prdt = stack(v_)
dprdt = stack(d_)


a1_2 = plot(prdt[1,:], prdt[2,:], legend = :best, label = "Recovered system", xlabel = L"V", ylabel = L"I")
v1 = plot(DATA.V[1:1500], legend = :best, color = :black, label = "V")
plot!(v1, prdt[1,:], ylabel = L"V", color = :blue, style = :dash, label = "predicted V")
plot!(v1, DATA.Vr[1:1500], color = :red, label = "Vr(t)")
i1 = plot(prdt[2,:], ylabel = L"I", xlabel = L"t", legend = :none)
plot!(i1, DATA.I[1:1500], color = :blue, style = :dash)
a3 = plot(
    v1, i1
    , layout = (2,1), xformatter = x -> x*dt, size = (800,600)
); png("a3 true.png")

[prdt[1,1:150] DATA.V[1:150] dprdt[2,1:150] DATA.dI[1:150] DATA.Vr[1:150] abs.(prdt[1,1:150] - DATA.V[1:150])]

abs.(dprdt[1,:] - DATA.dV[1:length(dprdt[1,:])])
plot((prdt[1,:] .> DATA.Vr[1:length(dprdt[1,:])]) .- (DATA.V[1:length(dprdt[1,:])] .> DATA.Vr[1:length(dprdt[1,:])]))
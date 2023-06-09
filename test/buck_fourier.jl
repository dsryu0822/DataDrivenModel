using ProgressBars
packages = [:DataFrames, :CSV, :LinearAlgebra, :Plots, :Flux, :Clustering, :JLD2, :LaTeXStrings, :CUDA]
@showtime for package in ProgressBar(packages)
    @eval using $(package)
end
println(join(packages, ", "), " loaded!")

include("../src/DDM.jl")
include("../src/ML.jl")
include("../src/nonsmooth.jl")
include("../src/ODEdata.jl")
default(size = (600,600), color = :black, legend = :none)

struct Fourier
    an
    bn
    L
end
Fourier(m::Integer) = Fourier(zeros(Float32, m), zeros(Float32, m), Float32[0.0])
function (layer::Fourier)(x)
    x = x .|> Float32
    m = length(layer.an)
    nxL⁻¹ = (1:m) * x[end, :]' .* exp(layer.L[1])
    return [x[1:end-1, :]; sum(
        [cospi.(nxL⁻¹) .* layer.an
       ; sinpi.(nxL⁻¹) .* layer.bn], dims = 1)]
end
Flux.@functor Fourier

function Base.show(io::IO, l::Fourier)
    print(io, "Fourier(", length(l.an), ")")
end
# vcat(collect(Flux.params(Fourier(50)))...)

## Data Load
DATA = CSV.read("data/buck.csv", DataFrame)
Y = select(DATA, [:dV, :dI]) |> Matrix # .|> Float32
X = select(DATA, [ :V,  :I]) |> Matrix # .|> Float32
XY = [X Y]

## Clustering
dbs = dbscan(col_normalize(Y)', 0.01); nsubsys = length(dbs.clusters); println(nsubsys, " clusters found!")
plot(DATA.V, DATA.I, color = dbs.assignments, alpha = 0.5)
bit_ = [dbs.assignments .== k for k in 1:nsubsys]
Θ_ = [poly_basis(X[s,:], 2)     for s in bit_]
Y_ = [Y[vcat(findall(s)...),:] for s in bit_]
Ξ_ = [STLSQ(Θ_[s], Y_[s], 0.01) for s in 1:nsubsys]

a1 = plot(DATA.V, DATA.I, alpha = 0.1, legend = :best, label = "Trajectory", xlabel = L"V", ylabel = L"I")
scatter!(a1, 
    DATA.V[Not(1)][.!iszero.(diff(dbs.assignments))]
  , DATA.I[Not(1)][.!iszero.(diff(dbs.assignments))]
  , color = 1, shape = :+, label = "Jumping points")

function foo(v)
    v = deepcopy(v)
    s = Int64(pop!(v))
    Θv = poly_basis(v, 2)'
    return vec([(Θv*Ξ_[s]) 0])
end

## Classification
# data = Float32.(Matrix(DATA[Not(end), [:V, :I, :t]])') # |> gpu
data = Float32.(Matrix(DATA[Not(end), [:V, :I, :dV, :dI, :t]])') # |> gpu
subs = Flux.onehotbatch(dbs.assignments[Not(1)], 1:nsubsys) # |> gpu
trng = Flux.DataLoader((data,subs), shuffle = true, batchsize = 1000)


p = size(data, 1)
SSSf = Chain( # SubSystemSelector
    Flux.Scale(5),
    Fourier(50),
    Dense( p => 50, relu),
    Dense(50 => 50, relu),
    Dense(50 => 50, relu),
    Dense(50 => nsubsys),
    softmax
) # |> gpu
Loss(x,y) = Flux.crossentropy(SSSf(x), y)
losses = [Loss(data, subs)]
optimizer = ADAM()

norm_variables = Float32.(norm.(eachrow(Matrix(data)))); norm_variables[end] = 1.0
SSSf[1].scale ./= norm_variables

if !isfile("test/SSSf.jld2")
    print("test/SSSf.jld2 not exists, ")
    
    ps = Flux.params(SSSf)
    println("Training..., Loss: ", losses[1])
    @time for epch in 1:100_000
        if epch < 10
            @time Flux.train!(Loss, ps, trng, optimizer)
        else
            Flux.train!(Loss, ps, trng, optimizer)
        end
        loss = Loss(data, subs)
        acc = sum(dbs.assignments[Not(1)] .== argmax.(eachcol(SSSf(data)))) / length(dbs.assignments[Not(1)])
        if loss < losses[end]
        # if acc > maximum(acc_)
            jldsave("C:/Temp/SSSf.jld2"; SSSf, loss, acc, losses)
            println()
            print("epoch ", lpad(length(losses), 5)
                    , ": loss = ", rpad(trunc(loss, digits = 6), 8)
                    , ", acc = ", rpad(acc, 7)
                    , " saved!")
            push!(losses, loss)
        else
            epch % 100 == 0 && print("epch- ", lpad(epch, 5), ", ")
            push!(losses, losses[end])
        end
        # push!(acc_, acc)
    end
    cp("C:/Temp/SSSf.jld2", "test/SSSf.jld2")
else
    @info "Loading SSSf..."
    fSSS = jldopen("test/SSSf.jld2")
    SSSf = fSSS["SSSf"]
    println("loss = ", fSSS["loss"])
    println("acc = ", fSSS["acc"])
    close(fSSS)
end

### Classifier testing
a1_1 = scatter(a1, 
    DATA.V[Not(1, 2)][.!iszero.(diff(argmax.(eachcol(SSSf(data)))))]
  , DATA.I[Not(1, 2)][.!iszero.(diff(argmax.(eachcol(SSSf(data)))))]
  , color = 2, shape = :+, label = "Detected Jumping points")
png(a1_1, "a1_1.png")

## ODE recovery
v_ = [[XY[1, 1:2]; dbs.assignments[1]]]
d_ = [foo(v_[end])]
dt = 0.00001
for (tk, t) in ProgressBar(enumerate(0:dt:0.025))
    push!(v_, RK4(foo, v_[end], dt))
    push!(d_, foo(v_[end]))
    # v_[end][end] = Float32[v_[end][Not(end)]; t] |> SSSf |> vec |> argmax
    v_[end][end] = Float32[v_[end][Not(end)]; d_[end][Not(end)]; t] |> SSSf |> vec |> argmax
end

### Trajectories
prdt = stack(v_)
a1_2 = plot(prdt[1,:], prdt[2,:], legend = :best, label = "Recovered system", xlabel = L"V", ylabel = L"I")
png(a1_2, "a1_2.png")

v1 = plot(prdt[1,:], ylabel = L"V")
plot!(v1, DATA.V[1:2500], color = :blue)
i1 = plot(prdt[2,:], ylabel = L"I", xlabel = L"t")
plot!(i1, DATA.I[1:2500], color = :blue)
a3 = plot(
    v1, i1
    , layout = (2,1), xformatter = x -> x*dt, legend = :none, size = (800,600)
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

v_ = [[XY[1, 1:2]; dbs.assignments[1]]]
d_ = [foo(v_[end])]
dt = 0.00001
for (tk, t) in ProgressBar(enumerate(dt:dt:0.005))
    v_[tk][end] = v_[tk][1] > DATA.Vr[tk] ? 1 : 2
    push!(v_, RK4(foo, v_[tk], dt))
    push!(d_, foo(v_[tk])) 
end
prdt = stack(v_)
dprdt = stack(d_)

tk = 21
v_[tk]
# v_[tk][1] > DATA.Vr[tk] ? 1 : 2
# DATA.Vr[tk]
# foo(v_[tk])
RK4(foo, [v_[tk][Not(end)]; 1.0], dt)
RK4(foo, [v_[tk][Not(end)]; 2.0], dt)
foo([v_[tk][Not(end)]; 2.0])
buck([v_[tk][Not(end)]; 0.0003])

a1_2 = plot(prdt[1,:], prdt[2,:], legend = :best, label = "Recovered system", xlabel = L"V", ylabel = L"I")
v1 = plot(prdt[1,:], ylabel = L"V", color = :blue)
plot!(v1, DATA.V[1:500], color = :black)
plot!(v1, DATA.Vr[1:500], color = :red)
i1 = plot(prdt[2,:], ylabel = L"I", xlabel = L"t")
plot!(i1, DATA.I[1:500], color = :blue)
a3 = plot(
    v1, i1
    , layout = (2,1), xformatter = x -> x*dt, legend = :none, size = (800,600)
)

Ξ_[1] - Ξ_[2]
Ξ_[2][1,2] - E/L

[prdt[1,1:50] prdt[2,1:50] DATA.V[1:50] dprdt[2,1:50] DATA.dI[1:50] DATA.Vr[1:50] abs.(prdt[1,1:50] - DATA.V[1:50])]


abs.(dprdt[1,:] - DATA.dV[1:length(dprdt[1,:])])
plot((prdt[1,:] .> DATA.Vr[1:length(dprdt[1,:])]) .- (DATA.V[1:length(dprdt[1,:])] .> DATA.Vr[1:length(dprdt[1,:])]))
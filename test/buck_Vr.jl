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

struct Fourier
    an
    bn
    L
end
Fourier(m::Integer) = Fourier(zeros(Float32, m), zeros(Float32, m), Float32[-5.0])
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
Y = select(DATA, [:dV, :dI]) |> Matrix .|> Float32
X = select(DATA, [ :V,  :I]) |> Matrix .|> Float32
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
data = Float32.(Matrix(DATA[Not(end), [:V, :I, :dV, :dI, :Vr]])') # |> gpu
subs = Flux.onehotbatch(dbs.assignments[Not(1)], 1:nsubsys) # |> gpu
trng = Flux.DataLoader((data,subs), shuffle = true, batchsize = 1000)


p = size(data, 1)
SSS = Chain( # SubSystemSelector
    Flux.Scale(5),
    # Fourier(50),
    Dense( p => 50, relu),
    Dense(50 => 50, relu),
    Dense(50 => 50, relu),
    Dense(50 => nsubsys),
    softmax
) # |> gpu
Loss(x,y) = Flux.crossentropy(SSS(x), y)
losses = [Loss(data, subs)]
optimizer = ADAM()

norm_variables = Float32.(norm.(eachrow(Matrix(data)))); norm_variables[end] = 1.0
SSS[1].scale ./= norm_variables

if !isfile("test/SSS.jld2")
    print("test/SSS.jld2 not exists, ")
    
    ps = Flux.params(SSS)
    println("Training..., Loss: ", losses[1])
    @time for epch in 1:100_000
        if epch < 10
            @time Flux.train!(Loss, ps, trng, optimizer)
        else
            Flux.train!(Loss, ps, trng, optimizer)
        end
        loss = Loss(data, subs)
        acc = sum(dbs.assignments[Not(1)] .== argmax.(eachcol(SSS(data)))) / length(dbs.assignments[Not(1)])
        if loss < losses[end]
        # if acc > maximum(acc_)
            jldsave("C:/Temp/SSS.jld2"; SSS, loss, acc, losses)
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
    cp("C:/Temp/SSS.jld2", "test/SSS.jld2")
else
    @info "Loading SSS..."
    fSSS = jldopen("test/SSS.jld2")
    SSS = fSSS["SSS"]
    println("loss = ", fSSS["loss"])
    println("acc = ", fSSS["acc"])
    close(fSSS)
end

### Classifier testing
a1_1 = scatter(a1, 
    DATA.V[Not(1, 2)][.!iszero.(diff(argmax.(eachcol(SSS(data)))))]
  , DATA.I[Not(1, 2)][.!iszero.(diff(argmax.(eachcol(SSS(data)))))]
  , color = 2, shape = :+, label = "Detected Jumping points")
png(a1_1, "a1_1.png")

## ODE recovery
v_ = [[XY[1, 1:2]; dbs.assignments[1]]]
d_ = [foo(v_[end])]
dt = 10^(-7)
for (tk, t) in ProgressBar(enumerate(0:dt:0.025))
    push!(v_, RK4(foo, v_[end], dt))
    push!(d_, foo(v_[end]))
    # v_[end][end] = Float32[v_[end][Not(end)]; t] |> SSS |> vec |> argmax
    v_[end][end] = Float32[v_[end][Not(end)]; d_[end][Not(end)]; DATA.Vr[tk]] |> SSS |> vec |> argmax
end

### Trajectories
prdt = stack(v_)
a1_2 = plot(a1_1, prdt[1,:], prdt[2,:], legend = :best, label = "Recovered system", xlabel = L"V", ylabel = L"I")
png(a1_2, "a1_2.png")

a3 = plot(
    plot(prdt[1,:], ylabel = L"V"), plot(prdt[2,:], ylabel = L"I", xlabel = L"t")
    , layout = (2,1), xformatter = x -> x*dt, legend = :none
)
png(a3, "a3.png")
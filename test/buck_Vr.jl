using ProgressBars
packages = [:DataFrames, :CSV, :LinearAlgebra, :Plots, :Flux, :Clustering, :JLD2, :LaTeXStrings]
@time for package in ProgressBar(packages)
    @eval using $(package)
end
include("../src/ML.jl")
include("../src/DDM.jl")
include("../src/nonsmooth.jl")
include("../src/ODEdata.jl")
default(size = (600,600), color = :black, legend = :none)

DATA = CSV.read("data/buck.csv", DataFrame)

Y = select(DATA, [:dV, :dI]) |> Matrix .|> Float32
X = select(DATA, [ :V,  :I]) |> Matrix .|> Float32
XY = [X Y]

dbs = dbscan(col_normalize(Y)', 0.01); nsubsys = length(dbs.clusters); println(nsubsys, " clusters found!")
plot(DATA.V, DATA.I, color = dbs.assignments, alpha = 0.5)

## Clustering
bit_ = [dbs.assignments .== k for k in 1:nsubsys]
Θ1 = poly_basis(X[bit_[1],:], 2)
Θ2 = poly_basis(X[bit_[2],:], 2)

Θ = [Θ1 zeros(size(Θ1))
     zeros(size(Θ2)) Θ2]
Ξ = STLSQ(Θ, Y[vcat(findall.(bit_)...),:], 0.01)

function foo(v)
    v = deepcopy(v)
    s = Int64(pop!(v))
    w = poly_basis(v, 2)
    ph = zeros(6, 2)
    ph[:,s] = w
    return [vec(reshape(ph, 1, :) * Ξ); 0]
end

## Classification
data = Float32.(Matrix(DATA[Not(end), Not(:t)])')
subs = Flux.onehotbatch(dbs.assignments[Not(1)], 1:nsubsys) # |> gpu
trng = Flux.DataLoader((data,subs) , batchsize = 1000)

norm_variables = Float32.(norm.(eachrow(data)))
p = size(data, 1)
SSS = Chain( # SubSystemSelector
    x -> x./norm_variables,
    Dense( p => 50),    relu,
    Dense(50 => 50),    relu,
    Dense(50 => 50),    relu,
    Dense(50 => nsubsys),
    softmax
) # |> gpu
Loss(x,y) = Flux.crossentropy(SSS(x), y)
optimizer = ADAM()

if !isfile("test/SSS.jld2")
    losses = [Loss(data, subs)]
    acc_ = [0.0]
    ps = Flux.params(SSS)
    println("Training...")
    @time for epch in 1:100_000
        epch % 100 == 0 && println("epoch ", lpad(epch, 5))
        Flux.train!(Loss, ps, trng, optimizer)
        loss = Loss(data, subs)
        acc = sum(dbs.assignments[Not(1)] .== argmax.(SSS.(eachcol(data)))) / length(dbs.assignments[Not(1)])
        if loss < minimum(losses)
        # if acc > maximum(acc_)
            jldsave("C:/Temp/SSS.jld2"; SSS, loss, acc)
            println("epoch ", lpad(epch, 5)
                    , ": loss = ", rpad(trunc(loss, digits = 6), 8)
                    , ", acc = ", rpad(acc, 7)
                    , " saved!")
        end
        push!(losses, loss)
        push!(acc_, acc)
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

a1 = plot(DATA.V, DATA.I, alpha = 0.1, legend = :best, label = "Trajectory", xlabel = L"V", ylabel = L"I")
scatter!(a1, 
    DATA.V[Not(1)][.!iszero.(diff(dbs.assignments))]
  , DATA.I[Not(1)][.!iszero.(diff(dbs.assignments))]
  , color = 1, shape = :+, label = "Jumping points")
scatter!(a1, 
    DATA.V[Not(1, 2)][.!iszero.(diff(argmax.(SSS.(eachcol(data)))))]
  , DATA.I[Not(1, 2)][.!iszero.(diff(argmax.(SSS.(eachcol(data)))))]
  , color = 2, shape = :+, label = "Detected Jumping points")
png(a1, "a1.png")

DATA


## ODE recovery
v_ = [[XY[1, 1:2]; dbs.assignments[1]]]
d_ = [foo(v_[end])]
dt = 0.00001
for (tk, t) in ProgressBar(enumerate(dt:dt:0.025))
    push!(v_, RK4(foo, v_[end], dt))
    push!(d_, foo(v_[end]))
    v_[end][end] = argmax(SSS(Float32[v_[end][Not(end)]; d_[end][Not(end)]; DATA.Vr[tk]]))
    # if argmax(SSS(Float32[d_[end][Not(end)]; v_[end][Not(end)]])) == 2
    #     v_[end][end] = mod(v_[end][end], 2) + 1
    # end
end

prdt = stack(v_)
a2 = plot(prdt[1,:], prdt[2,:], legend = :best, label = "Recovered system", xlabel = L"V", ylabel = L"I")
png(a2, "a2.png")

a3 = plot(
    plot(prdt[1,:], ylabel = L"V"), plot(prdt[2,:], ylabel = L"I", xlabel = L"t")
    , layout = (2,1), xformatter = x -> x*dt, legend = :none
)
png(a3, "a3.png")
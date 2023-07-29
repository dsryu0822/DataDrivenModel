@time using DataFrames, CSV, Plots
using LinearAlgebra
using Plots: mm, cm
using LaTeXStrings
default(size = (600,600))

DATA = CSV.read("data/buck.csv", DataFrame)
# plot(plot.(eachcol(DATA[1:250,:]))...)
# plot(DATA.V, DATA.I)

Y = select(DATA, [:dV, :dI]) |> Matrix
X = select(DATA, [ :V,  :I]) |> Matrix
YX = [Y X]


using Clustering, Distances
a1 = plot(eachcol(YX[:,3:4])..., color = :black, title = "Phase plane (normalized)", label = :none, xlabel = L"V", ylabel = L"I", alpha = 0.5)

dbs = dbscan(col_normalize([Y X])', 0.01); nclusters = length(dbs.clusters); println(nclusters, " clusters found!")
a1 = plot(eachcol(YX[:,3:4])..., color = dbs.assignments, title = "Phase plane", label = :none, xlabel = L"V", ylabel = L"I", alpha = 0.5)
scatter(eachcol(YX[:,1:2])..., color = dbs.assignments, title = "Phase plane of dV, dI (normalized)", label = :none, xlabel = L"dV", ylabel = L"dI")

include("../src/DDM.jl")
include("../src/ODEdata.jl")
bit_ = [dbs.assignments .== k for k in 1:nclusters]
Θ1 = poly_basis(YX[bit_[1], 3:4], 2)
Θ2 = poly_basis(YX[bit_[2], 3:4], 2)

using SparseArrays

Θ = [Θ1 zeros(size(Θ1))
     zeros(size(Θ2)) Θ2]
Ξ = STLSQ(Θ, YX[vcat(findall.(bit_)...),1:2], 0.01)

col_norm = Float32.(norm.(eachcol(YX)))
function mynormlize(x)
    return x ./ col_norm
end

using Flux, CUDA, JLD2
CUDA.functional()
FFNN = Chain(
    mynormlize,
    Dense( 4 => 50),    relu,
    Dense(50 => 50),    relu,
    Dense(50 => 50),    relu,
    Dense(50 =>  2), softmax
) # |> gpu
L(x,y) = Flux.crossentropy(FFNN(x), y)
optimizer = ADAM()

yx = Float32.(YX')[:, Not(end)] # |> gpu
# ss = Flux.onehotbatch(dbs.assignments[Not(1)], 1:2) # |> gpu
ss = Flux.onehotbatch(abs.(diff(dbs.assignments)) .+ 1, 1:2)
trng = Flux.DataLoader((yx,ss) , batchsize = 1000)

if !isfile("test/FFNN.jld2")
    losses = [L(yx, ss)]
    acc_ = [0.0]
    ps = Flux.params(FFNN)
    println("Training...")
    @time for epch in 1:100000
        epch % 100 == 0 && println("epoch ", lpad(epch, 5))
        Flux.train!(L, ps, trng, optimizer)
        loss = L(yx, ss)
        acc = sum(dbs.assignments[Not(1)] .== argmax.(FFNN.(eachcol(yx)))) / length(dbs.assignments[Not(1)])
        if loss < minimum(losses)
        # if acc > maximum(acc_)
            jldsave("C:/Temp/FFNN.jld2"; FFNN, loss, acc)
            println("epoch ", lpad(epch, 5)
                    , ": loss = ", rpad(trunc(loss, digits = 6), 8)
                    , ", acc = ", rpad(acc, 7)
                    , " saved!")
        end
        push!(losses, loss)
        push!(acc_, acc)
    end
    cp("C:/Temp/FFNN.jld2", "test/FFNN.jld2")
else
    @info "Loading FFNN..."
    fFFNN = jldopen("test/FFNN.jld2")
    FFNN = fFFNN["FFNN"]
    println("loss = ", fFFNN["loss"])
    println("acc = ", fFFNN["acc"])
    close(fFFNN)
end

bit_TP = dbs.assignments[Not(1)] .== argmax.(FFNN.(eachcol(yx)))
sum(bit_TP) / length(dbs.assignments[Not(1)])
a2 = scatter(a1, 
    eachcol(YX[Not(1),3:4][.!bit_TP,:])...,
    color = :black, shape = :x, label = "Error in classification")

bit_trans = argmax.(FFNN.(eachcol(yx))) .== 2
a3 = scatter(a1, 
    eachcol(YX[Not(1),3:4][bit_trans,:])...,
    color = :black, shape = :+, label = "Detected jumping points")


function foo(v)
    v = deepcopy(v)
    s = Int64(pop!(v))
    w = poly_basis(v, 2)
    ph = zeros(6, 2)
    ph[:,s] = w
    return [vec(reshape(ph, 1, :) * Ξ); 0]
end

using ProgressBars
v_ = [[YX[1, 3:4]; dbs.assignments[1]]]
d_ = [foo(v_[end])]
dt = 10^(-7)
for t in ProgressBar(dt:dt:0.025)
    push!(v_, RK4(foo, v_[end], dt))
    push!(d_, foo(v_[end]))
    # v_[end][end] = argmax(FFNN(Float32[d_[end][Not(end)]; v_[end][Not(end)]]))
    if argmax(FFNN(Float32[d_[end][Not(end)]; v_[end][Not(end)]])) == 2
        v_[end][end] = mod(v_[end][end], 2) + 1
    end
end
prdt = stack(v_)
plot(prdt[1,:], prdt[2,:], label = "Recovered system", xlabel = L"V", ylabel = L"I")
default(color = :black)
plot(
    plot(prdt[1,:], ylabel = L"V"), plot(prdt[2,:], ylabel = L"I", xlabel = L"t")
    , layout = (2,1), xformatter = x -> x*dt, legend = :none
)

plot(DATA.V[1:2500], DATA.I[1:2500], label = "Buck converter", xlabel = L"V", ylabel = L"I")
plot(
    plot(DATA.V[1:2500], ylabel = L"V"), plot(DATA.I[1:2500], ylabel = L"I", xlabel = L"t")
    , layout = (2,1), xformatter = x -> x*dt, legend = :none
)

sum(last.(d_))
sum(Int64.(last.(v_)) .== 2)
@time using DataFrames, CSV, Plots
using LinearAlgebra
using Plots: mm, cm
using LaTeXStrings

DATA = CSV.read("data/buck.csv", DataFrame)
# plot(plot.(eachcol(DATA[1:250,:]))...)
# plot(DATA.V, DATA.I)

Y = select(DATA, [:dV, :dI]) |> Matrix
X = select(DATA, [ :V,  :I]) |> Matrix

YX = [Y X]
YX = YX ./ norm.(eachcol(YX))'

using Clustering, Distances

dbs = dbscan(YX', 0.01); println(length(dbs.clusters), " clusters found!")
# dbs = dbscan(YX', 0.01, metric = Cityblock()); println(length(dbs.clusters), " clusters found!")
a1 = plot(eachcol(YX[:,3:4])..., color = dbs.assignments, title = "Phase plane (normalized)", label = :none, xlabel = L"V", ylabel = L"I", alpha = 0.5)
plot(eachcol(YX[:,1:2])..., color = dbs.assignments, title = "Phase plane of dV, dI (normalized)", label = :none, xlabel = L"dV", ylabel = L"dI")

using Flux, CUDA, JLD2

FFNN = Chain(
    Dense(4,50, relu),
    Dense(50,50, relu),
    Dense(50,2),
    softmax
)
L(x,y) = Flux.crossentropy(FFNN(x), y)
optimizer = ADAM()

yx = Float32.(YX[Not(end),:])'
ss = Flux.onehotbatch(dbs.assignments[Not(1)], 1:2)
trng = Flux.DataLoader((yx,ss), batchsize = 256)

if !isfile("test/FFNN.jld2")
    losses = [L(yx, ss)]
    for epch in 1:10000
        Flux.train!(L, Flux.params(FFNN), trng, optimizer)
        loss = L(yx, ss)
        if loss < minimum(losses)
            jldsave("C:/FFNN.jld2"; FFNN)
            println("epoch ", lpad(epch, 5), ": loss = ", loss, " saved!")
        end
        push!(losses, loss)
    end
    cp("C:/FFNN.jld2", "test/FFNN.jld2")
else
    @info "Loading FFNN..."
    fFFNN = jldopen("test/FFNN.jld2")
    FFNN = fFFNN["FFNN"]
    println("loss = ", Flux.crossentropy(FFNN(yx), ss))
    close(fFFNN)
end

sum(dbs.assignments[Not(1)] .== argmax.(FFNN.(eachcol(yx)))) / length(dbs.assignments[Not(1)])

a2 = scatter(a1, 
    eachcol(YX[Not(1),3:4][dbs.assignments[Not(1)] .!== argmax.(FFNN.(eachcol(yx))),:])...,
    color = :black, shape = :x, label = "predicted switching points")
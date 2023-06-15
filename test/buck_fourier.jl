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


struct Fourier
    an
    bn
end
Fourier(m::Integer) = Fourier(zeros(m), zeros(m))
function (layer::Fourier)(x::AbstractVector)
    approximated = [x[1:(end-1)];
    sum([ (layer.an .* cos.((1:length(layer.an)) .* x[end]))
    ; (layer.bn .* sin.((1:length(layer.bn)) .* x[end]))
    ])]
    return approximated .|> Float32
end
(layer::Fourier)(x::AbstractMatrix) = hcat(layer.(eachcol(x))...)
Flux.@functor Fourier

function Base.show(io::IO, l::Fourier)
    print(io, "Fourier(", length(l.an), ")")
end

DATA = CSV.read("data/buck.csv", DataFrame)

Y = select(DATA, [:dV, :dI]) |> Matrix .|> Float32
X = select(DATA, [ :V,  :I]) |> Matrix .|> Float32
XY = [X Y]

dbs = dbscan(col_normalize(Y)', 0.01); nsubsys = length(dbs.clusters); println(nsubsys, " clusters found!")
plot(DATA.V, DATA.I, color = dbs.assignments, alpha = 0.5)

## Clustering
bit_ = [dbs.assignments .== k for k in 1:nsubsys]
Θ_ = [poly_basis(X[s,:], 2)     for s in bit_]
Y_ = [Y[vcat(findall(s)...),:] for s in bit_]
Ξ_ = [STLSQ(Θ_[s], Y_[s], 0.01) for s in 1:nsubsys]

function foo(v)
    v = deepcopy(v)
    s = Int64(pop!(v))
    Θv = poly_basis(v, 2)
    return vec([(Θv*Ξ_[s]) 0])
end

## Classification
data = Float32.(Matrix(DATA[Not(end), [:V, :I, :dV, :dI, :t]])') # |> gpu
subs = Flux.onehotbatch(dbs.assignments[Not(1)], 1:nsubsys) # |> gpu
trng = Flux.DataLoader((data,subs) , batchsize = 100) # |> gpu

norm_variables = Float32.(norm.(eachrow(data))); norm_variables[end] = 1
p = size(data, 1)
SSSf = Chain( # SubSystemSelector
    x -> x ./ norm_variables,
    Fourier(50),
    Dense( p => 50),    relu,
    Dense(50 => 50),    relu,
    Dense(50 => 50),    relu,
    Dense(50 => nsubsys),
    softmax
) # |> gpu
Loss(x,y) = Flux.crossentropy(SSSf(x), y)
optimizer = ADAM()

if !isfile("test/SSSf.jld2")
    losses = [Loss(data, subs)]
    acc_ = [0.0]
    ps = Flux.params(SSSf)
    println("Training...")
    @time for epch in 1:100_000
        epch % 100 == 0 && println("epoch ", lpad(epch, 5))
        Flux.train!(Loss, ps, trng, optimizer)
        loss = Loss(data, subs)
        acc = sum(dbs.assignments[Not(1)] .== argmax.(SSSf.(eachcol(data)))) / length(dbs.assignments[Not(1)])
        if loss < minimum(losses)
        # if acc > maximum(acc_)
            jldsave("C:/Temp/SSSf.jld2"; SSSf, loss, acc)
            println("epoch ", lpad(epch, 5)
                    , ": loss = ", rpad(trunc(loss, digits = 6), 8)
                    , ", acc = ", rpad(acc, 7)
                    , " saved!")
        end
        push!(losses, loss)
        push!(acc_, acc)
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

a1 = plot(DATA.V, DATA.I, alpha = 0.1, legend = :best, label = "Trajectory", xlabel = L"V", ylabel = L"I")
scatter!(a1, 
    DATA.V[Not(1)][.!iszero.(diff(dbs.assignments))]
  , DATA.I[Not(1)][.!iszero.(diff(dbs.assignments))]
  , color = 1, shape = :+, label = "Jumping points")
scatter!(a1, 
    DATA.V[Not(1, 2)][.!iszero.(diff(argmax.(SSSf.(eachcol(data)))))]
  , DATA.I[Not(1, 2)][.!iszero.(diff(argmax.(SSSf.(eachcol(data)))))]
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
    v_[end][end] = argmax(SSSf(Float32[v_[end][Not(end)]; d_[end][Not(end)]; DATA.Vr[tk]]))
    # if argmax(SSSf(Float32[d_[end][Not(end)]; v_[end][Not(end)]])) == 2
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
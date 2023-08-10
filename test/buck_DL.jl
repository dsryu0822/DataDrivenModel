@time include("buck_.jl")

subsystem = (DATA.dI .> 0) .+ 1

_DATA = deepcopy(DATA)
_DATA[!, :now] = subsystem
__DATA = _DATA[1:1:end,:]
__next = __DATA.now[Not(1)]
__DATA = __DATA[Not(end), :]
__DATA[!, :next] = __next
_DATA = __DATA
__DATA = Nothing

## Classification
const label = deepcopy(_DATA.now)
gx = Float32.(Matrix(select(_DATA, [:V, :I, :Vr]))') |> gpu
gy = Flux.onehotbatch(label, 1:nsubsys) |> gpu
gdata = Flux.DataLoader((gx,gy), shuffle = true, batchsize = 500_000)
# data = gdata

ANN = init_ANN(gdata)
ANN = save_ANN(ANN, gdata; lastepch = 100_000)
# ANN, loss, acry, epch = load_ANN("C:/Temp/SSS_000023.jld2")
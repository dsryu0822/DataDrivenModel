@time include("buck_.jl")

Stalinize(arr, m::Real, M::Real) = arr[m .≤ arr .≤ M]
Stalinize(arr, mM) = arr[first(mM) .≤ arr .≤ last(mM)]

subsystem = (DATA.dI .> 0) .+ 1

_DATA = deepcopy(DATA)
_DATA[!, :now] = subsystem
__DATA = _DATA[1:1000:end,:]
__next = __DATA.now[Not(1)]
__DATA = __DATA[Not(end), :]
__DATA[!, :next] = __next

# indx_switch = findall(__DATA.now .!= __DATA.next)
# indx_switch = sort(Stalinize(vcat([indx_switch .+ k for k in -1:1]...), 1, 5000000))
# __DATA = __DATA[indx_switch, :]

_DATA = __DATA
__DATA = Nothing

## Classification
const label = deepcopy(_DATA.now)
gx = Float32.(Matrix(select(_DATA, [:V, :Vr]))') |> gpu
gy = Flux.onehotbatch(label, 1:nsubsys) |> gpu
gdata = Flux.DataLoader((gx,gy), shuffle = true, batchsize = min(size(gx, 2), 500_000))
# data = gdata

# ANN = load_ANN("C:/Temp/SSS.jld2")
ANN = init_ANN(gdata)
ANN = save_ANN(ANN, gdata; lastepch = 100_000)
# ANN, loss, acry, epch = load_ANN("C:/Temp/SSS_000023.jld2")


using DecisionTree

labels = _DATA.now
features = Matrix(_DATA[:, [:V, :Vr]])
model = build_forest(labels, features, 2, 1, 0.9, 2)
n_folds=5; n_subfeatures=2
nfoldCV_forest(labels, features, n_folds, n_subfeatures)

using XGBoost

bst = xgboost((features, labels), num_round = 10)

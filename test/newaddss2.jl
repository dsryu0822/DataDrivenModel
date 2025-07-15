include("../core/header.jl")


import Base.rand
rand(df::AbstractDataFrame; n = 1) = df[rand(1:nrow(df), n), :]
doublerange(n) = Base.product(1:n, 1:n)
centerpick(n, m) = round.(Int64, range(1, n, m+2))[2:end-1]
has_neighbor(arr) = !isempty(intersect(arr, [arr .+ 1; arr .- 1]))

##########################################################################
#                                                                        #
#                            Soft impact model                           #
#                                                                        #
##########################################################################
vrbl = [:dt, :du, :dv], [:t, :u, :v]
# cnfg = (; f_ = [cospi], λ = 1e-1)
dt = 1e-5; tspan = [0, 20]; #θ = 1e-16;
if isfile("hyperparameter/0 data.csv")
    data = CSV.read("hyperparameter/0 data.csv", DataFrame)
    test = CSV.read("hyperparameter/0 test.csv", DataFrame)
else
    # data = factory_soft(DataFrame, 0.1, ic = [0, .0446272, -0.119564]; tspan, dt)
    # test = factory_soft(DataFrame, 0.1, ic = [0, .000129715, 0.301469]; tspan, dt)
    temp = factory_soft(DataFrame, 0.1, ic = [0, .0446272, -0.119564]; tspan, dt)
    data = temp[1:(nrow(temp) ÷ 2), :]
    test = temp[(nrow(temp) ÷ 2):end, :]
    CSV.write("hyperparameter/0 data.csv", data, bom = true)
    CSV.write("hyperparameter/0 test.csv", test, bom = true)
end

jumpt = detect_jump(data, vrbl)
sets = set_divider(jumpt)
datasets = [data[set, :] for set in sets]
ground_truth = [ifelse(rand(ds.u) > 0.05, 3, ifelse(rand(ds.u) < -0.05, 2, 1)) for ds in datasets]
ground_truth_ = [findall(ground_truth .== k) for k in eachindex(1:3)]
ground_truth__ =[any([[i,j] ⊆ gt for gt in ground_truth_]) for (i,j) in doublerange(length(datasets))]
sampled = [ds[centerpick(nrow(ds), 50), :] for ds in datasets]

θ1 = 1e-0
B_ = vec([(; N, M, λ = 0) for (N, M) in Base.product(0:3, 0:3)])
f__ = [[SINDy(smp, vrbl...; B...) for smp in sampled] for B in B_]
ranks__ = [[rank(Θ(smp[:, last(vrbl)]; B...)) for smp in sampled] for B in B_]
fullrank = first.([[f.sparse_matrix.m for f in f_] for f_ in f__])
f__[6][5] |> print

i = 6
j = 7
sum.(abs2, residual.(f__[i], Ref(sampled[j])))

ts_mse(f, ds) = norm(solve(f, [ds[1, last(vrbl)]...], 0:dt:500dt) - Matrix(ds[1:501, last(vrbl)]))
ts_mse(f__[3][6], datasets[7])


SINDy(datasets[1], vrbl...; B_[6]..., λ = 1e-2) |> print
SINDy(datasets[24], vrbl...; B_[6]..., λ = 1e-2) |> print
Θ(datasets[1][:, last(vrbl)]; B_[6]...)
plot(datasets[1].u, legend = :none)
plot(cospi.(datasets[1].u), legend = :none)
DataFrame(
    실제_데이터_길이 = nrow.(datasets),
    랭크 = [rank(Θ(ds[:, last(vrbl)]; B_[6]...)) for ds in datasets],
)


_S = [[i == j ? Inf : SINDy([sampled[i]; sampled[j]], vrbl...; B_[6]...).MSE for i in eachindex(sampled)] for j in eachindex(sampled)]
S_ = vcat.(sampled, sampled[argmin.(_S)])
DataFrame(
    샘플_데이터_길이 = "50 + 50",
    랭크 = [rank(Θ(S[:, last(vrbl)]; B_[6]...)) for S in S_],
)
SINDy(S_[24], vrbl...; B_[6]..., λ = 1e-2) |> print

function add_fold!(data::AbstractDataFrame, k = 5)
    data[!, :fold] .= 1 .+ mod.(shuffle(1:nrow(data)), k)
    return data
end
add_fold!
folder = [rand(datasets[1], n = 20) for _ in 1:10]
f_ = [SINDy(fold, vrbl...; B_[6]...) for fold in folder]





n0 = 176
_data = data[1:n0, :]
cnfg = cook(last(vrbl), poly = 0:2)
M = Θ(_data[:, last(vrbl)], cnfg)
R = rankAnalysis2(M)
corM = rankAnalysis3(M)
itr = 0
for _ in 1:nrow(cnfg)
    itr += 1
    println("="^60, "\nitration: $itr")
    R = rankAnalysis2(M)
    println("rank / term = $(1 + maximum(R)) / $(nrow(cnfg))")
    if (1 + maximum(R)) == nrow(cnfg)
        @info "all terms are linearly independent"
        break
    end
    suspect = findall(R .== maximum(R))
    del1, del2 = argmax(corM).I
    println("suspect: $(cnfg.term[del1]) and $(cnfg.term[del2]) with correlation $(corM[del1, del2])")
    f1_ = [SINDy(_data[_data.fold .!= k, :], vrbl, cnfg[Not(del1), :]) for k in 1:5]
    f2_ = [SINDy(_data[_data.fold .!= k, :], vrbl, cnfg[Not(del2), :]) for k in 1:5]
    mmsw1 = sum.(abs2, [residual(f1_[k], _data[_data.fold .== k, :]) for k in 1:5])
    mmse2 = sum.(abs2, [residual(f2_[k], _data[_data.fold .== k, :]) for k in 1:5])
    del = mmsw1 < mmse2 ? del1 : del2
    println("delete suspect: $(cnfg.term[del])")
    cnfg = cnfg[Not(del), :]
    M = M[:, Not(del)]
    corM = corM[Not(del), Not(del)]
end
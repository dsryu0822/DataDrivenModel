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

N_ = 0:3
M_ = 0:3
schedules = DataFrame(ID = Int64[], N = Int64[], M = Int64[], worked = String[], ss = Int64[], E = Float64[])
for (N, M) in Base.product(N_, M_) push!(schedules, (0, N, M, "", 0, 0.0, )) end
schedules.ID = 1:size(schedules, 1)

# schedules.cut = 

# 이거 내가 직접 눈으로 보고 넣어야함

manualθ = [
    1 1e+0
    2 1e-0
    3 1e-0
    4 1e-5
    5 1e-0
    6 1e-20
    7 1e-20
    8 1e-10
    9 1e-7
   10 1e-20
   11 1e-20
   12 1e-10
   13 1e-10
   14 1e-20
   15 1e-20
   16 1e-20
]

schedules.θ = manualθ[:, 2]
for dr = eachrow(schedules)
    ID = dr.ID; N = dr.N; M = dr.M; θ = dr.θ
    cnfg = (; N, M, λ = 0)

    n = length(datasets)
    min_m = 2SINDy(datasets[1][1:10, :], vrbl...; cnfg...).sparse_matrix.m
    _sampled = [ds[centerpick(nrow(ds), min_m), :] for ds in datasets]
    f__ = fill(SINDy(datasets[1][1:10, :], vrbl...; cnfg...), n, n)
    mse__ = fill(Inf, n, n)
    for (i,j) = doublerange(length(datasets))
        if i ≤ j continue end
        f__[i, j] = SINDy([_sampled[i]; _sampled[j]], vrbl...; cnfg...)
        f__[j, i] = f__[i, j]
        mse__[i, j] = f__[i, j].MSE
        mse__[j, i] = mse__[i, j]
    end
    sampled = [[_sampled[k]; _sampled[argmin(mse__[k,:])]] for k in eachindex(_sampled)]
    mse__ = fill(Inf, n, n)
    for (i,j) = doublerange(length(datasets))
        if i ≤ j continue end
        f__[i, j] = SINDy([sampled[i]; sampled[j]], vrbl...; cnfg...)
        f__[j, i] = f__[i, j]
        mse__[i, j] = f__[i, j].MSE
        mse__[j, i] = mse__[i, j]
    end

    subsys = connected_components(SimpleGraph(mse__ .< θ))
    f_ = [SINDy([sampled[ss]...;], vrbl...; cnfg...) for ss in subsys]; # print.(f_)
    open("hyperparameter/$ID subsys.txt", "w") do io
        println(io, "N = $N, M = $M")
        for k in eachindex(f_)
            f = SINDy([sampled[subsys[k]]...;], vrbl...; cnfg..., λ = 1e-3)
            println(io, pretty_table(f))
        end
    end

    feasibility = !any((length.(subsys) .== 1) .|| has_neighbor.(subsys))
    B_ = [SINDy(sampled[k], vrbl...; cnfg...) for k in eachindex(sampled)]
    E = sum(sum.([[norm(f_[i].sparse_matrix - B_[j].sparse_matrix) for j in subsys[i]] for i in eachindex(f_)]))

    dr .= [ID, N, M, ifelse(feasibility, "✅", "❌"), length(subsys), E, θ]
    
    scatter(vec(mse__[0 .< mse__ .< Inf]), yscale = :log10, legend = :none,
            title = "N = $N, M = $M", ylabel = "MSE", xlabel = "Pairwise sampled data",
            size = [600, 1200], left_margin = 10mm, color = [:white, :black][vec(ground_truth__[0 .< mse__ .< Inf]) .+ 1])
    hline!([θ], color = :red)
    png("hyperparameter/$ID mse n=$(N)_m=$(M).png")
end
schedules

default()
performance = plot()
for m = 0:3
    plot!(performance, schedules[schedules.M .== m, :].N, schedules[schedules.M .== m, :].E, 
            label = "M = $m", xlabel = "N", ylabel = "E",
            lw = 2, size = [400, 400])
end
plot(performance, yscale = :log10, legend = :topleft)

##########################################################################
#                                                                        #
#                             Gear system                                #
#                                                                        #
##########################################################################
dt = 1-2; tspan = [0, 200]
data = factory_gear(DataFrame, -0.2; ic = [0, .1, .1, .1], tspan)
vrbl = [:dx, :dv, :dΩ, :dθ], [:x, :v, :Ω, :θ]

jumpt = detect_jump(data, vrbl; dos = 1)
sets = set_divider(jumpt)
datasets = [data[set, :] for set in sets]
ground_truth = [rand(ds.x) < 1 ? rand(ds.x) > -1 ? 1 : 2 : 3 for ds in datasets]

N_ = 0:2
M_ = 0:2
C_ = 1:2
schedules = DataFrame(ID = Int64[], N = Int64[], M = Int64[], C = Int64[], worked = String[], ss = Int64[], E = Float64[])
for (N, M, C) in Base.product(N_, M_, C_) push!(schedules, (0, N, M, C, "", 0, 0.0, )) end
schedules.ID = 1:size(schedules, 1)

# manualθ = [
#     1 1e-20
#     2 1e-3
#     3 1e-10
#     4 1e-10
#     5 1e-20
#     6 1e-25
#     7 1e-26
#     8 1e-30
#     9 1e-12
#    10 1e-2
#    11 1e-2
#    12 1e-10
#    13 1e-10
#    14 1e-15
#    15 1e-15
#    16 1e-15
#    17 1e-15
#    18 1e-15
# ]
manualθ = [
    1 1e-2
    2 1e-3
    3 1e-5
    4 1e-5
    5 1e-20
    6 1e-25
    7 1e-10
    8 1e-25
    9 1e-25
   10 1e-2
   11 1e-3
   12 1e-5
   13 1e-5
   14 1e-20
   15 1e-28
   16 1e-10
   17 1e-15
   18 1e+0
]
schedules.θ = manualθ[:, 2]

schedules = schedules[1:16, :]
for dr = eachrow(schedules)
    ID = dr.ID; N = dr.N; M = dr.M; C = dr.C; θ = dr.θ
    sin2(x) = sin(2x); sin3(x) = sin(3x); sin4(x) = sin(4x);
    cos2(x) = cos(2x); cos3(x) = sin(3x); cos4(x) = cos(4x);
    sincos = [sin, cos, sin2, cos2, sin3, cos3, sin4, cos4]
    cnfg = (; N, f_ = sincos[1:(2M)], C, λ = 0)

    # cnfg = (; N = 1, M = 1, λ = 0)
    # # min_m = minimum(nrow.(datasets))
    n = length(datasets)
    min_m = SINDy(datasets[1][1:10, :], vrbl...; cnfg...).sparse_matrix.m
    # sampled = rand.(datasets; n = min_m)
    sampled = [ds[centerpick(nrow(ds), min_m), :] for ds in datasets]
    f__ = fill(SINDy(datasets[1][1:10, :], vrbl...; cnfg...), n, n)
    mse__ = fill(Inf, n, n)
    dist__ = zeros(n, n)
    for (i,j) = doublerange(length(datasets))
        if i ≤ j continue end
        f__[i, j] = SINDy([sampled[i]; sampled[j]], vrbl...; cnfg...)
        f__[j, i] = f__[i, j]
        mse__[i, j] = f__[i, j].MSE
        if abs(i - j) ≤ 1 mse__[i, j] = Inf end
        mse__[j, i] = mse__[i, j]
    end
    # B_ = getproperty.(f__[argmin(mse__, dims = 2)], :sparse_matrix)
    # for (i,j) = doublerange(length(datasets))
    #     if i ≤ j continue end
    #     mse__[i, j] = norm(B_[i] - B_[j])
    #     if abs(i - j) ≤ 1 mse__[i, j] = Inf end
    #     mse__[j, i] = mse__[i, j]
    # end

    # subsys = connected_components(SimpleGraph(mse__ .< θ))
    # f_ = [SINDy([sampled[ss]...;], vrbl...; cnfg...) for ss in subsys]
    # # f_ = [SINDy([sampled[ss]...;], vrbl...; cnfg..., λ = 1e-3) for ss in subsys]; print.(f_)

    # feasibility = !any((length.(subsys) .== 1) .|| has_neighbor.(subsys))
    # B_ = [SINDy([sampled[k]; sampled[argmin(mse__[k,:])]], vrbl...; cnfg...) for k in eachindex(sampled)]
    # E = mean(mean.([[norm(f_[i].sparse_matrix - B_[j].sparse_matrix) for j in subsys[i]] for i in eachindex(f_)]))

    # dr .= [ID, N, M, C, ifelse(feasibility, "✅", "❌"), length(subsys), E, θ]    

    scatter(filter(!iszero, mse__), yscale = :log10,
            title = "N = $N, M = $M, C = $C", ylabel = "MSE", xlabel = "Pairwise sampled data",
            size = [600, 1200], left_margin = 10mm)
    hline!([θ], color = :red)
    png("hyperparameter/$ID mse n=$(N)_m=$(M)_c=$(C).png")
end
schedules
f__[end] |> print
cnfg = (; N = 1, f_ = [sin, cos], C = 2)
# cnfg = (; N = 1, M = 1)
SINDy([datasets[ground_truth .== 3]...;], vrbl...; cnfg..., λ = 1e-2) |> print
[findall(ground_truth .== 1),
findall(ground_truth .== 2),
findall(ground_truth .== 3)]


##########################################################################
#                                                                        #
#                           Hindmarsh-Rose model                         #
#                                                                        #
##########################################################################
data = factory_hrnm(DataFrame, 0.1)
vrbl = [:dt, :dx, :dy, :dz], [:t, :x, :y, :z]

jumpt = detect_jump(data, vrbl)
sets = set_divider(jumpt)
datasets = [data[set, :] for set in sets]
ground_truth = [ifelse(rand(ds.z) > 1, 3, ifelse(rand(ds.z) < -1, 2, 1)) for ds in datasets]
ground_truth_ = [findall(ground_truth .== k) for k in eachindex(1:3)]
ground_truth__ =[any([[i,j] ⊆ gt for gt in ground_truth_]) for (i,j) in doublerange(length(datasets))]

N_ = 0:3
M_ = 0:3
schedules = DataFrame(ID = Int64[], N = Int64[], M = Int64[], worked = String[], ss = Int64[], E = Float64[])
for (N, M) in Base.product(N_, M_) push!(schedules, (0, N, M, "", 0, 0.0, )) end
schedules.ID = 1:size(schedules, 1)

manualθ = [
    1 1e-0
    2 1e-0
    3 1e-5
    4 1e-10
    5 1e-3
    6 1e-5
    7 1e-10
    8 1e-20
    9 1e-3
   10 1e-5
   11 1e-10
   12 1e-25
   13 1e-5
   14 1e-5
   15 1e-10
   16 1e-25
]

schedules.θ = manualθ[:, 2]
for dr = eachrow(schedules)
    ID = dr.ID; N = dr.N; M = dr.M; θ = dr.θ
    sin2(x) = sin(2x); sin3(x) = sin(3x); sin4(x) = sin(4x);
    cos2(x) = cos(2x); cos3(x) = sin(3x); cos4(x) = cos(4x);
    sincos = [sin, cos, sin2, cos2, sin3, cos3, sin4, cos4]
    cnfg = (; N, f_ = sincos[1:(2M)], λ = 0)

    n = length(datasets)
    min_m = 2SINDy(datasets[1][1:10, :], vrbl...; cnfg...).sparse_matrix.m
    _sampled = [ds[centerpick(nrow(ds), min_m), :] for ds in datasets]
    f__ = fill(SINDy(datasets[1][1:10, :], vrbl...; cnfg...), n, n)
    mse__ = fill(Inf, n, n)
    for (i,j) = doublerange(length(datasets))
        if i ≤ j continue end
        f__[i, j] = SINDy([_sampled[i]; _sampled[j]], vrbl...; cnfg...)
        f__[j, i] = f__[i, j]
        mse__[i, j] = f__[i, j].MSE
        mse__[j, i] = mse__[i, j]
    end
    sampled = [[_sampled[k]; _sampled[argmin(mse__[k,:])]] for k in eachindex(_sampled)]
    mse__ = fill(Inf, n, n)
    for (i,j) = doublerange(length(datasets))
        if i ≤ j continue end
        f__[i, j] = SINDy([sampled[i]; sampled[j]], vrbl...; cnfg...)
        f__[j, i] = f__[i, j]
        mse__[i, j] = f__[i, j].MSE
        mse__[j, i] = mse__[i, j]
    end

    subsys = connected_components(SimpleGraph(mse__ .< θ))
    f_ = [SINDy([sampled[ss]...;], vrbl...; cnfg...) for ss in subsys]; # print.(f_)
    open("hyperparameter/$ID subsys.txt", "w") do io
        println(io, "N = $N, M = $M")
        for k in eachindex(f_)
            f = SINDy([sampled[subsys[k]]...;], vrbl...; cnfg..., λ = 1e-3)
            println(io, pretty_table(f))
        end
    end

    feasibility = !any((length.(subsys) .== 1) .|| has_neighbor.(subsys))
    B_ = [SINDy(sampled[k], vrbl...; cnfg...) for k in eachindex(sampled)]
    E = sum(sum.([[norm(f_[i].sparse_matrix - B_[j].sparse_matrix) for j in subsys[i]] for i in eachindex(f_)]))

    dr .= [ID, N, M, ifelse(feasibility, "✅", "❌"), length(subsys), E, θ]
    
    scatter(vec(mse__[0 .< mse__ .< Inf]), yscale = :log10, legend = :none,
            title = "N = $N, M = $M", ylabel = "MSE", xlabel = "Pairwise sampled data",
            size = [600, 1200], left_margin = 10mm, color = [:white, :black][vec(ground_truth__[0 .< mse__ .< Inf]) .+ 1])
    hline!([θ], color = :red)
    png("hyperparameter/$ID mse n=$(N)_m=$(M).png")
end
sort(schedules, :E)

cnfg = (; N = 3, f_ = sincos[1:2], λ = 1e-3)
SINDy([datasets[ground_truth .== 1]...;], vrbl...; cnfg..., λ = 1e-3) |> print

SINDy([sampled[1]; sampled[2]], vrbl...; cnfg...)
subsys = connected_components(SimpleGraph(mse__ .< 5e-30))

SINDy(data, vrbl...; cnfg..., λ = 1e-2)

scatter(data.z, data.dz, shape = :pixel)
scatter(test.u[1:100:end], test.v[1:100:end], shape = :pixel, xlims = [-.1, .1])

mse__ = fill(Inf, n, n)
for (i,j) = doublerange(length(datasets))
    if i ≤ j continue end
    f__[i, j] = SINDy([sampled[i]; sampled[j]], vrbl...; cnfg...)
    f__[j, i] = f__[i, j]
    if abs(i - j) == 1
        mse__[i, j] = Inf
        mse__[j, i] = Inf
        continue
    end
    mse__[i, j] = f__[i, j].MSE
    mse__[j, i] = mse__[i, j]
end


spshow(A) = DataFrame(replace(A, -1 => ""), string.(1:n))

M = deepcopy(mse__)
A = [abs(i - j) == 1 ? 0 : -1 for (i,j) in doublerange(length(datasets))]
spshow(A)

@assert count(A .== -1) > n
# i, j = argmax(replace(M, Inf => 0)).I
# A[i,j] = 0; M[i,j] = Inf;
# A[j,i] = 0; M[j,i] = Inf;

i, j = argmin(M).I
A[i,j] = 1; M[i,j] = Inf;
A[j,i] = 1; M[j,i] = Inf;
for bit = 0:1
    for k in findall(A[i, :] .== bit)
        if j == k continue end
        A[k, j] = bit; M[k, j] = Inf
        A[j, k] = bit; M[j, k] = Inf
    end
    for k in findall(A[j, :] .== bit)
        if i == k continue end
        A[k, i] = bit; M[k, i] = Inf
        A[i, k] = bit; M[i, k] = Inf
    end
end
print(spshow(A))
scatter(vec(mse__[0 .< M .< Inf]), yscale = :log10,
title = "N = 3, M = 1", ylabel = "MSE",
size = [600, 1200], left_margin = 10mm, color = [:white, :black][vec(ground_truth__[0 .< M .< Inf]) .+ 1])
@info ""

ground_truth_
connected_components(SimpleGraph(A))

bag = [[k] for k in eachindex(ground_truth)]
out = [[k-1, k+1] for k in eachindex(ground_truth)]

M = deepcopy(mse__)

i,j = argmin(M).I
M[i,j] = Inf; M[j,i] = Inf;
_i = findfirst(i .∈ bag)
_j = findfirst(j .∈ bag)
_i, _j = sort([_i, _j])
if _i != _j
    if isempty(bag[_i] ∩ out[_j]) && isempty(bag[_j] ∩ out[_i])
        push!(bag, unique(sort([bag[_i]; bag[_j]])))
        push!(out, unique(sort([out[_i]; out[_j]])))
        deleteat!(bag, _j); deleteat!(bag, _i);
        deleteat!(out, _j); deleteat!(out, _i);
    else
        @warn "!!!"
    end
end
bag .=> out
scatter(vec(mse__[0 .< M .< Inf]), yscale = :log10,
title = "N = 3, M = 1", ylabel = "MSE",
size = [600, 1200], left_margin = 10mm, color = [:white, :black][vec(ground_truth__[0 .< M .< Inf]) .+ 1])

ground_truth_
plot(datasets[12].z)
plot(datasets[13].z)
vec(mse__)
vec(mse__)

[1:19 nrow.(datasets)]
f__[6, 12]
f__[1, 16]
plot(datasets[16].z, size = [400, 400], legend = :none)
plot(data.z, size = [800, 300], legend = :none)

plot(data.x, data.y, data.z, xlabel = "x", ylabel = "y", zlabel = "z",
     lw = 2, legend = :none, color = :black, size = [400, 400])
scatter(data.z[1:10:end], data.dz[1:10:end], xlabel = "z", ylabel = "dz",
     lw = 2, legend = :none, color = :black, size = [400, 400], shape = :pixel)
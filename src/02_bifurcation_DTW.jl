function dtw(ts1::AbstractVector, ts2::AbstractVector; dist = (a, b) -> abs(a - b))
    n, m = length(ts1), length(ts2)

    # 1) 누적 비용 행렬
    D = fill(Inf, n + 1, m + 1)
    D[1, 1] = 0.0
    for i in 1:n, j in 1:m
        cost = dist(ts1[i], ts2[j])
        D[i+1, j+1] = cost + min(D[i, j+1], D[i+1, j], D[i, j])
    end

    # 2) 역추적으로 정렬 경로(path) 복원
    i, j = n, m
    path = Tuple{Int,Int}[]
    while i > 0 && j > 0
        push!(path, (i, j))
        choices = (D[i, j], D[i, j+1], D[i+1, j])
        _, idx = findmin(choices)
        if idx == 1
            i -= 1; j -= 1
        elseif idx == 2
            i -= 1
        else
            j -= 1
        end
    end
    reverse!(path)
    L = length(path)

    # 3) 경로 상의 위치(k번째 대응)를 공통 [0,1] 정렬 시간축으로 사용
    p = L > 1 ? [(k - 1) / (L - 1) for k in 1:L] : [0.0]

    # 4) 각 인덱스 i(또는 j)가 경로에서 등장한 위치들의 평균을 그 인덱스의 [0,1] 좌표로 사용
    sum1, cnt1 = zeros(n), zeros(Int, n)
    sum2, cnt2 = zeros(m), zeros(Int, m)
    for (k, (ii, jj)) in enumerate(path)
        sum1[ii] += p[k]; cnt1[ii] += 1
        sum2[jj] += p[k]; cnt2[jj] += 1
    end

    t1 = sum1 ./ cnt1   # length(t1) == n == length(ts1)
    t2 = sum2 ./ cnt2   # length(t2) == m == length(ts2)

    return t1, t2
end

function w1(a::AbstractVector, b::AbstractVector)
    a = filter(isfinite, a)
    b = filter(isfinite, b)
    if isempty(a) || isempty(b)
        return Inf
    end
    q = range(0.0, 1.0, length = 512)
    return mean(abs.(quantile(a, q) .- quantile(b, q)))
end



bfcnA = callbfcn("G:/BF/foodchain/bfcnA.jld2")
lineA = sort(DataFrame(k = [keys(bfcnA)...], v = maximum.([values(bfcnA)...])), :k)
pltB_ = []
pltC_ = []
pltBbif_ = []
pltCbif_ = []
for (P0, P1) = [[93, 94], [93, 95], [94, 95], [93, 96], [94, 96], [95, 96]]
    bfcnB = callbfcn("G:/BF/foodchain/bfcnB_$(P0)$(P1).jld2")
    bfcnC = callbfcn("G:/BF/foodchain/bfcnC_$(P0)$(P1).jld2")
    # scatter(dict2bifurcation(bfcnA)..., msw = 0, ms = 1, color = :black)
    # scatter(dict2bifurcation(bfcnC)..., msw = 0, ms = 1, color = :blue)
    lineB = sort(DataFrame(k = [keys(bfcnB)...], v = maximum.([values(bfcnB)...])), :k)
    lineC = sort(DataFrame(k = [keys(bfcnC)...], v = maximum.([values(bfcnC)...])), :k)
    # plot(lineA.k, lineA.v, color = :black)
    # plot(lineC.k, lineC.v, color = :blue)

    tB1, tB2 = dtw(lineA.v, lineB.v)
    yA_ = [[bfcnA[lineA.k[t]] for t in 1:nrow(lineA)]...;]
    xA_ = [[fill(tB1[t], length(bfcnA[lineA.k[t]])) for t in 1:nrow(lineA)]...;]
    yB_ = [[bfcnB[lineB.k[t]] for t in 1:nrow(lineB)]...;]
    xB_ = [[fill(tB2[t], length(bfcnB[lineB.k[t]])) for t in 1:nrow(lineB)]...;]
    pltB = scatter(xA_, yA_, msw = 0, ms = 1, color = :black)
    scatter!(xB_, yB_, msw = 0, ms = 1, color = :red)
    push!(pltB_, pltB)
    
    tC1, tC2 = dtw(lineA.v, lineC.v)
    yA_ = [[bfcnA[lineA.k[t]] for t in 1:nrow(lineA)]...;]
    xA_ = [[fill(tC1[t], length(bfcnA[lineA.k[t]])) for t in 1:nrow(lineA)]...;]    
    yC_ = [[bfcnC[lineC.k[t]] for t in 1:nrow(lineC)]...;]
    xC_ = [[fill(tC2[t], length(bfcnC[lineC.k[t]])) for t in 1:nrow(lineC)]...;]
    pltC = scatter(xA_, yA_, msw = 0, ms = 1, color = :black)
    scatter!(xC_, yC_, msw = 0, ms = 1, color = :red)
    push!(pltC_, pltC)
    
    push!(pltBbif_, plot([w1(
        [(bfcnA[p] for p in lineA.k[t .≤ tB1 .< t + 1e-2])...;],
        [(bfcnB[p] for p in lineB.k[t .≤ tB2 .< t + 1e-2])...;]
    ) for t = 0:1e-2:0.99], color = :red, ylims = [0, 0.05]))
    push!(pltCbif_, plot([w1(
        [(bfcnA[p] for p in lineA.k[t .≤ tC1 .< t + 1e-2])...;],
        [(bfcnC[p] for p in lineC.k[t .≤ tC2 .< t + 1e-2])...;]
    ) for t = 0:1e-2:0.99], color = :blue, ylims = [0, 0.05]))
end
plot([pltB_; pltC_]..., ylims = [0.55, 0.8], ms = 0.2, msw = 0, layout = (4, 3), size = (600, 600)); png("temp1")
plot([pltBbif_; pltCbif_]..., ylims = [0, 0.05], layout = (4, 3), size = (600, 600)); png("temp2")


plot([w1(
        [(bfcnA[p] for p in lineA.k[t .≤ tB1 .< t + 1e-2])...;],
        [(bfcnB[p] for p in lineB.k[t .≤ tB2 .< t + 1e-2])...;]
    ) for t = 0:1e-2:0.99], color = :red, ylims = [0, 0.05], xlims = [70, 80])


pltB = scatter(xA_, yA_, msw = 0, ms = 1, color = :black, xlims = [0.7, 0.8])
scatter!(xB_, yB_, msw = 0, ms = 1, color = :red)


[w1(
        [(bfcnA[p] for p in lineA.k[t .≤ tB1 .< t + 1e-2])...;],
        [(bfcnB[p] for p in lineB.k[t .≤ tB2 .< t + 1e-2])...;]
    ) for t = 0:1e-2:0.99][76]

t=0.78

plot(lineA.k, lineB.k)




function dtw(X::AbstractVector, R::AbstractVector, Q::AbstractVector; dist = (a, b) -> abs(a - b))
    n, m = length(R), length(Q)
    @assert length(X) == n "length(X)는 length(R)과 같아야 합니다."

    # 1) 누적 비용 행렬
    D = fill(Inf, n + 1, m + 1)
    D[1, 1] = 0.0
    for i in 1:n, j in 1:m
        cost = dist(R[i], Q[j])
        D[i+1, j+1] = cost + min(D[i, j+1], D[i+1, j], D[i, j])
    end

    # 2) 역추적으로 정렬 경로(path) 복원
    i, j = n, m
    path = Tuple{Int,Int}[]
    while i > 0 && j > 0
        push!(path, (i, j))
        choices = (D[i, j], D[i, j+1], D[i+1, j])
        _, idx = findmin(choices)
        if idx == 1
            i -= 1; j -= 1
        elseif idx == 2
            i -= 1
        else
            j -= 1
        end
    end
    reverse!(path)

    # 3) Q의 각 인덱스 j가 경로에서 대응된 R의 인덱스 i들을 X 좌표로 변환해 평균
    sumX = zeros(m)
    cnt = zeros(Int, m)
    for (i, j) in path
        sumX[j] += X[i]
        cnt[j] += 1
    end
    _X = sumX ./ cnt   # length(_X) == length(Q) == m

    return _X
end


bfcnA = callbfcn("G:/BF/foodchain/bfcnA.jld2")
lineA = sort(DataFrame(k = [keys(bfcnA)...], v = maximum.([values(bfcnA)...])), :k)

bfcnB = callbfcn("G:/BF/foodchain/bfcnB_$(P0)$(P1).jld2")
bfcnC = callbfcn("G:/BF/foodchain/bfcnC_$(P0)$(P1).jld2")
# scatter(dict2bifurcation(bfcnA)..., msw = 0, ms = 1, color = :black)
# scatter(dict2bifurcation(bfcnC)..., msw = 0, ms = 1, color = :blue)
lineB = sort(DataFrame(k = [keys(bfcnB)...], v = maximum.([values(bfcnB)...])), :k)
lineC = sort(DataFrame(k = [keys(bfcnC)...], v = maximum.([values(bfcnC)...])), :k)

_XB = dtw(lineA.k, lineA.v, lineB.v)
scatter(dict2bifurcation(bfcnA)..., msw = 0, ms = 1, color = :black)
scatter!(
    [[fill(_XB[t], length(bfcnB[lineB.k[t]])) for t in eachindex(_XB)]...;],
    [[bfcnB[k] for k in lineB.k]...;],
    msw = 0, ms = 1, color = :red
)


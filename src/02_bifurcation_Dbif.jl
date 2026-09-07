include.("../core/" .* readdir("core")[[1,2,3,4,6]])

using MATLAB

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
DbifB_ = []
DbifC_ = []
for (P0, P1) = [[93, 94], [93, 95], [94, 95], [93, 96], [94, 96], [95, 96]]
    p0, p1 = P0/100, P1/100
    bfcnB = callbfcn("G:/BF/foodchain/bfcnB_$(P0)$(P1).jld2")
    bfcnC = callbfcn("G:/BF/foodchain/bfcnC_$(P0)$(P1).jld2")

    pm = 0.88; pM = 1.00;
    βm = (pm - p0) / (p1 - p0)
    βM = (pM - p0) / (p1 - p0)

    p_ = range(pm, pM, length = 2001)
    β_ = range(βm, βM, length = 2001)
    DbifB = fill(Inf, length(p_))
    DbifC = fill(Inf, length(p_))
    for k in eachindex(β_)
        DbifB[k] = try w1(bfcnA[p_[k]], bfcnB[β_[k]]) catch e Inf end
        DbifC[k] = try w1(bfcnA[p_[k]], bfcnC[β_[k]]) catch e Inf end
    end
    push!(DbifB_, DbifB)
    push!(DbifC_, DbifC)
end

mat"""
fig = tiledlayout(4, 3, 'Padding', 'compact');

nexttile
plot(linspace(-5, 7, 2001), $(DbifB_[1]), 'r')
xlim([-5, 7]); xticks([-5, 0, 1, 7]); ylim([0, 0.05]);

nexttile
plot(linspace(-2.5, 3.5, 2001), $(DbifB_[2]), 'r')
xlim([-2.5, 3.5]); xticks([-2.5, 0, 1, 3.5]); ylim([0, 0.05]);

nexttile
plot(linspace(-6, 6, 2001), $(DbifB_[3]), 'r')
xlim([-6, 6]); xticks([-6, 0, 1, 6]); ylim([0, 0.05]);

nexttile
plot(linspace(-1.7, 2.3, 2001), $(DbifB_[4]), 'r')
xlim([-1.7, 2.3]); xticks([-1.7, 0, 1, 2.3]); ylim([0, 0.05]);

nexttile
plot(linspace(-3, 3, 2001), $(DbifB_[5]), 'r')
xlim([-3, 3]); xticks([-3, 0, 1, 3]); ylim([0, 0.05]);

nexttile
plot(linspace(-7, 5, 2001), $(DbifB_[6]), 'r')
xlim([-7, 5]); xticks([-7, 0, 1, 5]); ylim([0, 0.05]);

nexttile
plot(linspace(-5, 7, 2001), $(DbifC_[1]), 'b')
xlim([-5, 7]); xticks([-5, 0, 1, 7]); ylim([0, 0.05]);

nexttile
plot(linspace(-2.5, 3.5, 2001), $(DbifC_[2]), 'b')
xlim([-2.5, 3.5]); xticks([-2.5, 0, 1, 3.5]); ylim([0, 0.05]);

nexttile
plot(linspace(-6, 6, 2001), $(DbifC_[3]), 'b')
xlim([-6, 6]); xticks([-6, 0, 1, 6]); ylim([0, 0.05]);

nexttile
plot(linspace(-1.7, 2.3, 2001), $(DbifC_[4]), 'b')
xlim([-1.7, 2.3]); xticks([-1.7, 0, 1, 2.3]); ylim([0, 0.05]);

nexttile
plot(linspace(-3, 3, 2001), $(DbifC_[5]), 'b')
xlim([-3, 3]); xticks([-3, 0, 1, 3]); ylim([0, 0.05]);

nexttile
plot(linspace(-7, 5, 2001), $(DbifC_[6]), 'b')
xlim([-7, 5]); xticks([-7, 0, 1, 5]); ylim([0, 0.05]);
"""
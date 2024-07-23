using ProgressMeter
packages = [:CSV, :DataFrames, :Plots]
@showprogress for package in packages
    @eval using $(package)
end

cd("G:/DDM/lyapunov")



@time files = CSV.read.(readdir("soft", join = true), DataFrame)
findall(nrow.(files) .< 11)
count(nrow.(files) .< 11)
scatter(nrow.(files))
histogram(nrow.(files))

lpnv1, lpnv2, lpnv3 = [], [], []
for k = 1:2001
    _t = files[k].t[end]
    temp = [sum(files[k].λ1) ./ _t, 
            sum(files[k].λ2) ./ _t, 
            sum(files[k].λ3) ./ _t]
    sort!(temp)
    push!(lpnv1, temp[1])
    push!(lpnv2, temp[2])
    push!(lpnv3, temp[3])
end
plot(size = [600, 300], legend = :none)
plot!(.1:.0001:.3, lpnv1, color = 1, ms = 1)
plot!(.1:.0001:.3, lpnv2, color = 2, ms = 1)
plot!(.1:.0001:.3, lpnv3, color = 3, ms = 1)
png("temp")
findall(lpnv1 .> -0.5)
# completed = eachrow.(files[nrow.(files) .== 6]);
completed = eachrow.(files);

temp = stack(getindex.(completed, nrow.(completed)))'
temp = stack(getindex.(completed, 6))'
temp[:, (end-2):end] ./= temp[:, 3]
for i in axes(temp, 1)
    temp[i, (end-2):end] .= sort(temp[i, (end-2):end])
end

scatter(temp[:, end])
scatter(temp[:, end], xlims = [1545, 1550], xticks = 1545:1550)

plot(size = [600, 300], legend = :none)
plot!(temp[:, end-2], color = 1, ms = 1)
plot!(temp[:, end-1], color = 2, ms = 1)
plot!(temp[:, end-0], color = 3, ms = 1)
plot!(xlims = [1560, 1570], xticks = 1560:1570)
png("temp")



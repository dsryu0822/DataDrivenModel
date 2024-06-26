using ProgressMeter
packages = [:CSV, :DataFrames, :Plots]
@showprogress for package in packages
    @eval using $(package)
end

cd("G:/DDM/lyapunov")



@time files = CSV.read.(readdir("soft", join = true), DataFrame)

scatter(nrow.(files))
histogram(nrow.(files))


# completed = eachrow.(files[nrow.(files) .== 6]);
completed = eachrow.(files);
temp = stack(getindex.(completed, nrow.(completed)))'
temp = stack(getindex.(completed, 6))'
temp[:, 6:8] ./= (temp[:, 3] .- 50)
for i in axes(temp, 1)
    temp[i, 6:8] .= sort(temp[i, 6:8])
end

scatter(temp[:, end])
scatter(temp[:, end], xlims = [1545, 1550], xticks = 1545:1550)

plot(size = [600, 300], legend = :none)
plot!(temp[:, 6], color = 1, ms = 1)
plot!(temp[:, 7], color = 2, ms = 1)
plot!(temp[:, 8], color = 3, ms = 1)
plot!(xlims = [1560, 1570], xticks = 1560:1570)
png("temp")
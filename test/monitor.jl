using ProgressMeter
packages = [:CSV, :DataFrames, :Plots]
@showprogress for package in packages
    @eval using $(package)
end

cd("G:/DDM/lyapunov")



@time files = CSV.read.(readdir("soft", join = true), DataFrame)

scatter(nrow.(files))
histogram(nrow.(files))


completed = eachrow.(files[nrow.(files) .== 52]);
temp = stack(getindex.(completed, 52))'
temp[:, 6:8] ./= 51
for i in axes(temp, 1)
    temp[i, 6:8] .= sort(temp[i, 6:8])
end

plot(size = [600, 300])
scatter!(temp[:, 2], temp[:, 6], color = 1, ms = 1)
scatter!(temp[:, 2], temp[:, 7], color = 2, ms = 1)
scatter!(temp[:, 2], temp[:, 8], color = 3, ms = 1)
png("temp")
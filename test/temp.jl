using Plots, CSV, DataFrames

lyapunov_table = CSV.read("G:/DDM/lyapunov/hrnm.csv", DataFrame)

plot(size = [600, 300], legend = :none)
plot!(lyapunov_table.f, lyapunov_table.λ1, color = 1, ms = 1)
plot!(lyapunov_table.f, lyapunov_table.λ2, color = 2, ms = 1)
plot!(lyapunov_table.f, lyapunov_table.λ3, color = 3, ms = 1)
plot!(lyapunov_table.f, lyapunov_table.λ4, color = 4, ms = 1)
png("temp")
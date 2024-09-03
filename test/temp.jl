include("../core/header.jl")

schedules = CSV.read("bifurcation/soft_schedules.csv", DataFrame)
schedules[!, :λ1] .= .0; schedules[!, :λ2] .= .0; schedules[!, :λ3] .= .0;

for dr = eachrow(schedules)
    filename = "lyapunov/soft_rcvd/$(lpad(dr.idx, 5, '0')).csv"
    data = CSV.read(filename, DataFrame)
    # dr[[:λ1, :λ2, :λ3]] .= sort(collect(data[end, (end-2):end]) / 100, rev=true)
    dr[[:λ1, :λ2, :λ3]] .= sort(collect(data[end, (end-2):end]) / 10, rev=true)
end

# CSV.write("lyapunov/soft.csv", schedules, bom = true)
CSV.write("lyapunov/soft_rcvd.csv", schedules, bom = true)
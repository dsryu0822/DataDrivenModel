include("../src/visual.jl")
include("../src/DDM.jl")
include("../src/ODEdata.jl")

# sunspot = CSV.read("data/sunspot_data.csv", DataFrame)[29951:end,:]
sunspot = CSV.read("data/sunspot_data.csv", DataFrame)[61000:end,:]

y = sunspot[:, 7]
plot(y)

df = DataFrame([Matrix(sunspot[Not(1), [1, 7]]) diff(Matrix(sunspot[:, [1, 7]]), dims = 1)], :auto)
STLSQ(df, [4], [1,2]; M = 3)

scatter(diff(y), ms = 0.4)
scatter(log10.(abs.((circshift(y, 1) ./ y) .- 1)), ms = 0.5)

plot(power.Voltage, power.Global_intensity, alpha = 0.1)

DATA = CSV.read("G:/earthquakes.csv", DataFrame)
earthquakes = DATA[:, [1, 7]]
earthquakes = DATA[occursin.("Japan", DATA.place), [1, 7]]
using StatsBase
earthquakes.time = earthquakes.time .÷ 10000000
eq = combine(groupby(earthquakes, :time), :time => length => "freq")
eq[!, "mgnt"] = combine(groupby(earthquakes, :time), :magnitudo => mean => "mgnt").mgnt
# plot(eq.freq)
# plot(eq.freq, xlims = (50000, 50000+100))
# plot(eq.mgnt)
# plot(log10.(eq.freq))
scatter(eq.freq, eq.mgnt, ms = 0.5)
plot(eq.freq, eq.mgnt, ms = 0.5)

scatter(eq.mgnt)
scatter(abs.(diff(eq.mgnt)))
scatter(abs.(diff(eq.freq)))

DATA = CSV.read("G:/emotions.csv", DataFrame)[1:end, :]
DATA = DATA[DATA.Activity .== 4, :]
plot(DATA.alx, DATA.aly, DATA.alz)


using EDF

DATA = EDF.read("G:/eeg/Subject00_1.edf")
DATA.io
DATA.signals[1] |> propertynames
plot(DATA.signals[1].samples)

DATA.signals[end].records
DATA.signals[end].samples_per_record
DATA.signals[end-1].header
EEG = stack(getproperty.(DATA.signals[1:(end-1)], :samples))[1:2000, :]

 plot(EEG[:, 1], color = 1, legend = :none)
plot!(EEG[:, 2], color = 2, legend = :none)
plot!(EEG[:, 4], color = 4, legend = :none)
plot!(EEG[:, 8], color = 8, legend = :none)
plot()

dt = 1/500
df = FDM1(EEG, dt)
sampled = 1:50
STLSQ_ = [STLSQ(df[sampled, :], 1:1, 22:22; M = 50)]


stranger = Int64[]
error_ = Float64[]
for k in ProgressBar(1:1999)
    error = sum(abs2, Matrix(df[[k], 1:21]) - (Θ(df[k, 22:42], M = 100) * first(STLSQ_).matrix))
    push!(error_, error)
    if error > 1e-5
        push!(stranger, k)
    end
end
scatter(log10.(error_))


function factory_STLSQ(STLSQed)
    function f(s, x)
        # return [1; vec(Θ(x, f_ = [cospi, sign, abs]) * STLSQed[s].matrix)]
        return vec(Θ(x, K = 2) * STLSQed[s].matrix)
    end
    return f
end
g = factory_STLSQ(STLSQ_)


x_ = [collect(df[1, 22:42])]
x = x_[end]

@time for t in ProgressBar(0:dt:4)
    s = 1
    x, dx = RK4(g, s, x_[end], dt)
    push!(x_, x)
end
_x_ = stack(x_)'


 plot(EEG[1:2000, 2])
plot!(_x_[1:2000, 2])



using Plots, CSV, DataFrames, LaTeXStrings, FFTW

# DATA = CSV.read("data/lorenz.csv", DataFrame)[1000:10:20000,:]

# plot(DATA.x1, DATA.x2, DATA.x3, color = :black, lw = 2, alpha = 0.5, formatter = (_...) -> "")
# scatter(DATA.x1, DATA.x2, DATA.x3, color = :black, lw = 2, marker = :x, alpha = 0.5, formatter = (_...) -> "")
# quiver(DATA.x1, DATA.x3, quiver = (diff(DATA.x1)/2, diff(DATA.x3)/2), size = (600,600), color = :black, xlabel = L"x(t)", ylabel = L"z(t)"); png("lorenz")

DATA = CSV.read("G:/buck/buck_000006.csv", DataFrame)

DATA[1:100000,:]
f = DATA.Vr[1:100000]
a1 = plot(abs.(f), label = "given data")

m = 10^(-3)
μ = 10^(-6)
Fs = 1/(400μ) #진동수
T = 0.1μ #샘플링 간격
L = length(f) #신호의 길이

ℱf = fft(f) # 푸리에 변환
ξ = Fs*[i for i in 0:L/2-1]/L #주파수 도메인(절반)

plot(
plot(ξ, abs.(ℱf[1:Int(L/2)])*2/L, title=L"Fourier transform of $f$")
, plot(ξ, log.(abs.(ℱf[1:Int(L/2)])*2/L), title=L"log Fourier transform of $f$")
, plot(ξ, abs.(ℱf[1:Int(L/2)])*2/L, title=L"Fourier transform of $f$", xlims = (0, 5))
, plot(ξ, log.(abs.(ℱf[1:Int(L/2)])*2/L), title=L"log Fourier transform of $f$", xlims = (0, 5))
, size = (1600, 900)
); png("a2")

plot(ξ, log.(abs.(ℱf[1:Int(L/2)])*2/L), title=L"log Fourier transform of $f$", xlims = (0, 100))

findall(Bool.(diff(DATA.dI .> 0)))


drop = -30; plot(abs.(ifft(ℱf[log10.(abs.(ℱf)) .> drop])), title = "drop: $drop")


f̂ = abs.(ifft(ℱf))
plot!(a1, f̂, label = "reconstructed")

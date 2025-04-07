include("../core/header.jl")

bfcn = vcat(CSV.read.(filter(x -> occursin("soft_bfcn", x), readdir()), DataFrame)...)
lpnv = vcat(CSV.read.(filter(x -> occursin("soft_lyapunov", x), readdir()), DataFrame)...)
sort!(bfcn, :idcs); sort!(lpnv, :idx)
lpnv = lpnv[.!iszero.(lpnv[:,3]), :]
done = unique(bfcn.idcs)

CSV.write("soft_bfcn_rcvd.csv", bfcn, bom = true)
CSV.write("soft_lyapunov_rcvd.csv", lpnv, bom = true)

plot()
plot!(lpnv.bp, lpnv.λ1, color = :red)
plot!(lpnv.bp, lpnv.λ2, color = :red)
plot!(lpnv.bp, lpnv.λ3, color = :red)
png("temp1")

scatter(bfcn.hrzn, bfcn.vrtc, color = :red,
msw = 0, ms = .5, ma = .5);
png("temp2")


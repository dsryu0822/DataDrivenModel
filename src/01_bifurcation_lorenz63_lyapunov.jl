include.("../core/" .* readdir("core")[[1,2,3,4,6]])


using Symbolics
jacobian(Matrix, f0)
J3 = jacobian(Function, f0)
J3(rand(3))

Float64.(substitute(J, Dict(rname .=> rand(3))))
foo2(rand(3))

substitute(J, Dict(rname .=> rand(3))) |> propertynames

J_func = build_function(J, rname, expression=Val{false})[1]
J_func(rand(3))  # 바로 Float64 행렬 반환

using ChaosTools

function lorenz_rule(u, p, t)
    σ = p[1]; ρ = p[2]; β = p[3]
    du1 = σ*(u[2]-u[1])
    du2 = u[1]*(ρ-u[3]) - u[2]
    du3 = u[1]*u[2] - β*u[3]
    return SVector{3}(du1, du2, du3)
end

lor = CoupledODEs(lorenz_rule, fill(10.0, 3), [10, 32, 8/3])
@time λλ = lyapunovspectrum(lor, 10000; Δt = 0.1)

function lorenz_rule(du, u, p, t)
    σ, ρ, β = p
    du[1] = σ*(u[2]-u[1])
    du[2] = u[1]*(ρ-u[3]) - u[2]
    du[3] = u[1]*u[2] - β*u[3]
end

u0 = fill(10.0, 3)   # 일반 Vector 사용 가능
lor = CoupledODEs(lorenz_rule, u0, [10, 32, 8/3])
@time λλ = lyapunovspectrum(lor, 10000; Δt = 0.1)

 
define(String, f0) |> print
foo = define(Function, f0)
lor = CoupledODEs(foo, [-38.8057, 5.32695, 192.705], [0.])
@time λλ = lyapunovspectrum(lor, 100000)



g0 |> print
s = g0

n = length(s.lname)+1
sz = nrow(s.recipe) ÷ n
Θx = s.recipe[1:sz, :tex]
matrices = [s.matrix[(i-1)*sz+1 : i*sz, :] for i in 1:n]
s_matrix = matrices[1]
ds_matrix = sum(matrices[2:end])
nmrt = Θx * s_matrix
dnmr = 1 .- (Θx * ds_matrix)
vec(nmrt ./ dnmr)

col = s_matrix[:, 1]
for j in eachindex(s.lname)
    s_col = s_matrix[:, j]
    ds_col = -ds_matrix[:, j]
    bit_snz = .!iszero.(s_col)
    bit_dsnz = .!iszero.(ds_col)
    @info join(string.(s_col[bit_snz].nzval) .* Θx[bit_snz], " + ")
    @info join(string.(ds_col[bit_dsnz].nzval) .* Θx[bit_dsnz], " + ")
end
s |> print
col[bit_nzterm] |> propertynames

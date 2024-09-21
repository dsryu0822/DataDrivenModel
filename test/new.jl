include("../core/header.jl")
vrbl = [:dt, :dx, :dy, :dz], [:t, :x, :y, :z]
cnfg = (; N = 3, f_ = [cos])
dt = 1e-4; θ1 = 1e-2; θ2 = 1e-27; θ3 = 1e-1; min_rank = 32;

data = CSV.read("lyapunov/hrnm_traj/00001.csv", DataFrame)

f_ = [SINDy(df, vrbl...; cnfg...) for df in groupby(data, :subsystem)]

typeof(s)
function functionalizer(s::STLSQresult)
    rname = eval(Meta.parse("@variables $(join(string.(s.rname), " "))"))
    fnexp = vec(sum(Θ(rname, N = s.N, M = s.M, f_ = s.f_, C = s.C)' .* s.matrix, dims = 1))
    return v -> substitute(fnexp, Dict(rname .=> v))
end
g_ = functionalizer.(f_)

using BenchmarkTools

x = collect(data[5, 1:4])
@btime f_[1]($x)
@btime g_[1]($x)

function Θ(X::AbstractMatrix;
    N = 1, M = 0, f_ = Function[], C = 1, λ = 0, sparse_rows = [])
    # λ is just for dummy argument for add_subsystem! function
    dim = size(X, 2)
    ansatz = []
    i = 0

    for k in 0:N
        for case = collect(multiexponents(dim, k))
            i += 1; i ∈ sparse_rows && continue
            push!(ansatz, prod(X .^ case', dims = 2))
        end
    end
    ΘX = hcat(ansatz...)
    for f in f_
        i += 1; i ∈ sparse_rows && continue
        ΘX = [ΘX f.(X)]
    end
    for m in 1:M
        i += 1; i ∈ sparse_rows && continue
        ΘX = [ΘX cospi.(m*X)]
    end
    for m in 1:M
        i += 1; i ∈ sparse_rows && continue
        ΘX = [ΘX sinpi.(m*X)]
    end

    dim = size(ΘX, 2)
    for c in 2:C
        for (j1, j2) in combinations(2:dim, c)
            i += 1; i ∈ sparse_rows && continue
            ΘX = [ΘX (ΘX[:, j1] .* ΘX[:, j2])]
        end
    end

    return ΘX
end

function slower(s, x)
    return vec(Θ(x; N = s.N, M = s.M, f_ = s.f_, C = s.C) * s.matrix)
end
function faster(s, x)
    return vec(Θ(x; N = s.N, M = s.M, f_ = s.f_, C = s.C, sparse_rows) * s_matrix)
end

s_matrix = f_[1].matrix[.!all.(map(x -> iszero.(x), eachrow(f_[1].matrix))),:]
sparse_rows = findall(all.(map(x -> iszero.(x), eachrow(f_[1].matrix))))

f_[1].matrix[sparse_rows, :]
findall(Nothing)


@btime slower(f_[1], $x)
@btime faster(f_[1], $x)


@btime 0.0 * 321545634184
@btime 97.0 * 321545634184
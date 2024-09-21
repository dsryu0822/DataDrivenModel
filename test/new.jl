include("../core/header.jl")
vrbl = [:dt, :dx, :dy, :dz], [:t, :x, :y, :z]
cnfg = (; N = 3, f_ = [cos])
dt = 1e-4; θ1 = 1e-2; θ2 = 1e-27; θ3 = 1e-1; min_rank = 32;

data = CSV.read("lyapunov/hrnm_traj/00001.csv", DataFrame)

@time f_ = [SINDy(df, vrbl...; cnfg...) for df in groupby(data, :subsystem)]

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


data = CSV.read("lyapunov/soft_traj/00001.csv", DataFrame)

vrbl = [:dt, :du, :dv], [:t, :u, :v]
cnfg = (; f_ = [cospi, sign], λ = 1e-2)
dt = 1e-6; θ1 = 1e-8; θ2 = 1e-12; θ3 = 1e-5; min_rank = 21;

f_ = [SINDy(df, vrbl...; cnfg...) for df in groupby(data, :subsystem)]
Dtree = dryad(data, last(vrbl)); # print_tree(Dtree)
J_ = jacobian.(Function, f_)

function lyapunov_exponent3(_data::DataFrame, J_, DT::Root{Float64, Int64}, bf_param;
    U = I(ncol(_data)), T = (last(_data.t) - first(_data.t)))
    last_vrbl = Symbol.(names(_data))

    λ = zeros(size(U, 1))
    for k = 1:nrow(_data)
        s = apply_tree(DT, collect(_data[k, :]))
        String_Dict = "substitute($(J_[s]), Dict($(join(["$vb => $(_data[k, vb])" for vb in last_vrbl], ", "))))"
        J = Float64.(eval(Meta.parse(String_Dict)))
        U, V = gram_schmidt(U)
        λ += V |> eachcol .|> norm .|> log
        U = RK4(J, U, dt)
    end
    return sort(λ / T, rev=true)
end
function lyapunov_exponent2(_data::DataFrame, J_, DT::Root{Float64, Int64}, bf_param;
    U = I(ncol(_data)), T = (last(_data.t) - first(_data.t)))

    λ = zeros(size(U, 1))
    for k = 1:nrow(_data)
        s = apply_tree(DT, collect(_data[k, :]))
        J = J_[s](collect(_data[k, :]))
        U, V = gram_schmidt(U)
        λ += V |> eachcol .|> norm .|> log
        U = RK4(J, U, dt)
    end
    return sort(λ / T, rev=true)
end
RK4(f_[1], [eps(), .05853, .47898], dt)

f_[1] |> propertynames
f_[1].dense_matrix
f_[1].sparse_rows
f_[1]([eps(), .05853, .47898])

s = f_[1]
Θ([eps(), .05853, .47898]; N = s.N, M = s.M, f_ = s.f_, C = s.C, sparse_rows = s.sparse_rows)
print(f_[1])
@time data = DataFrame(solve(f_, [eps(), .05853, .47898], dt, 0:dt:0.1, Dtree), last(vrbl))

@time λ = lyapunov_exponent1(data[:, last(vrbl)], jacobian.(Matrix, f_), Dtree, 0.1)
@time λ = lyapunov_exponent2(data[:, last(vrbl)], jacobian.(Function, f_), Dtree, 0.1)

J_ = jacobian.(Function, f_)
J_ = jacobian.(Matrix, f_)

function jacobian(T::Type, s::STLSQresult)
    rname = eval(Meta.parse("@variables $(join(string.(s.rname), " "))"))
    fnexp = vec(sum(Θ(rname, N = s.N, M = s.M, f_ = s.f_, C = s.C)' .* s.sparse_matrix, dims = 1))

    J = Symbolics.jacobian(fnexp, rname)
    if T == Matrix
        return J
    elseif T == Function
        return x -> Float64.(substitute(J, Dict(rname .=> x)))
    else
        error("Type not supported: Only `Function` or `Matrix` are supported.")
    end
end
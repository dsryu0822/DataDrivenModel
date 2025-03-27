include("../core/header.jl")

function lyapunov_exponent(_data::DataFrame, J_, DT, bf_param;
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
wsum(a, b, s = 0) = (1 - s) * a + s * b

function dt2df(dtree)
    io = IOBuffer()
    print_tree(io, dtree)
    str_dtree = replace(String(take!(io)), r"\w.\d+/\d+\w." => "",
                        "Feature " => "", "?" => "", " : " => "", "    "=>"\t", " "=>"")
    lines = split(str_dtree, "\n")[1:(end-1)]
    
    virtualdata = DataFrame([Float64[] for _ in last(vrbl)], last(vrbl))
    virtualdata.subsystem = Int64[]
    
    depth = (length.(findall.("\t", lines)) + length.(findall.("<", lines)))
    depth[1] = 0
    for dpth = sort(unique(depth), rev = true)
        bit_dpth = (depth .== dpth)
        leaf = strip.(lines[bit_dpth], Ref(['\t', '├', '└', '─']))
        featarg = parse(Int64, split(leaf[1], "<")[1])
        featval = parse(Float64, split(leaf[1], "<")[2])
        push!(virtualdata, [[ifelse(featarg == k, featval - 128eps(), 0) for k in eachindex(last(vrbl))]; parse(Int64, leaf[2])])
        push!(virtualdata, [[ifelse(featarg == k, featval + 128eps(), 0) for k in eachindex(last(vrbl))]; parse(Int64, leaf[3])])
        depth[bit_dpth] .= -1
        depth[findlast(bit_dpth)] = dpth - 1
    end
    sort!(virtualdata, last(vrbl))
    return virtualdata
end
function bifurcation_(sysname::AbstractString, data, d)
    if sysname == "soft"
        idx_sampled = diff(abs.(data.u) .> (d/2)) .< 0
        return data[Not(1), :v][idx_sampled]    
    elseif sysname == "gear"
        idx_sampled = diff([0; mod.(data.Ω, 2π)]) .< 0
        return data.v[idx_sampled]
    elseif sysname == "hrnm"
        idx_sampled = abs.(diff(diff(data.z) / 1e-3)) .> 0.1
        return data[Not([1, end]), :x][idx_sampled]
    end
end

##########################################################################
#                                                                        #
#                            Soft impact model                           #
#                                                                        #
##########################################################################
import DecisionTree.apply_tree
apply_tree(tree::Function, features::AbstractVector) = tree(features)
kozimo(bp) = (u -> ifelse(u[2] < -bp, 3, ifelse(u[2] > bp, 2, 1)))

sysname = "soft" # d ∈ [0.1, 0.3]
vrbl = [:dt, :du, :dv], [:t, :u, :v]
# cnfg = (; f_ = [cospi], λ = 1e-1) # λ = 5e-1 → 1e-2 → 1e-3
cnfg = (; N = 2, M = 2, f_ = [cospi], λ = 1e-1) # λ = 5e-1 → 1e-2 → 1e-3
dt = 1e-5; tspan = [30, 50]; θ = 1e-6;

# done = unique(vcat(CSV.read.(filter(x -> occursin("bfcn", x), readdir()), DataFrame)...).idcs)

_idx = 1:2001
tasks = Dict(device => _idx, "chaos1" => 1:650, "chaos2" => 651:1300, "chaos3" => 1301:length(_idx))
schedules = DataFrame(idx = _idx, bp = LinRange(0.1, 0.3, length(_idx)))[tasks[device], :]
for k in eachindex(last(vrbl)) schedules[!, "λ$k"] .= 0.0 end
bp1, bp2 = sort([0.100, 0.101]);

data_ = [factory_(sysname)(DataFrame, bp1; tspan, dt),
         factory_(sysname)(DataFrame, bp2; tspan, dt)]
for data in data_ add_subsystem!(data, vrbl, cnfg; θ) end
replace!(data_[2].subsystem, 2 => 3, 3 => 2)

f__ = [[SINDy(df, vrbl...; cnfg...) for df in groupby(data, :subsystem)] for data in data_]
_cnfg_ = [getproperty.(Ref(f_), f_ |> propertynames) for f_ in f__[1]];
Dtree_ = [dryad(data, last(vrbl)) for data in data_] # print_tree.(Dtree_)
df_Dtree_ = dt2df.(Dtree_);
for df_Dtree in df_Dtree_ df_Dtree.subsystem .= [3, 1, 1, 2] end

# jacobian.(Matrix, f__[1])[3]
J_ = [
    vv -> [0 0 0; 0 0 1; -π*sinpi(vv[1]) 0 0],
    vv -> [0 0 0; 0 0 1; -π*sinpi(vv[1]) -160000 -172.363],
    vv -> [0 0 0; 0 0 1; -π*sinpi(vv[1]) -160000 -172.363],
]
M1 = Matrix(df_Dtree_[1][:, 1:(end-1)])
M2 = Matrix(df_Dtree_[2][:, 1:(end-1)])

@info "$(now()): Preprocess for $(sysname) done!"

idcs = [Int64[] for _ in _idx]
hrzn = [Float64[] for _ in _idx]
vrtc = [Float64[] for _ in _idx]
@showprogress @threads for dr = eachrow(schedules)
# try
    pin = (dr.bp - bp1) / (bp2 - bp1)

    M0 = DataFrame(wsum(M1, M2, pin), last(vrbl))
    # M0.subsystem = df_Dtree_[2].subsystem
    # Dtree = dryad(M0, last(vrbl)) # print_tree(Dtree)
    Dtree = kozimo(M0.u[end])

    J_ = [
        vv -> [0 0 0 0; (-dr.bp*sin(vv[1])) (6vv[2] + 0.9vv[4] - 3*(vv[2]^2)) 1 (0.9vv[2]); 0 (-10vv[2]) -1 0; 0 0.8 0 -0.1],
        vv -> [0 0 0 0; (-dr.bp*sin(vv[1])) (6vv[2] + 0.9vv[4] - 3*(vv[2]^2)) 1 (0.9vv[2]); 0 (-10vv[2]) -1 0; 0 0.8 0 -0.1],
        vv -> [0 0 0 0; (-dr.bp*sin(vv[1])) (6vv[2] + 0.9vv[4] - 3*(vv[2]^2)) 1 (0.9vv[2]); 0 (-10vv[2]) -1 0; 0 0.8 0 -0.1],
    ]
    
    f_ = []
    for (f_1, f_2, _cnfg) = collect(zip(f__[1], f__[2], _cnfg_))
        __cnfg = collect(_cnfg)
        __cnfg[5] = wsum(f_1.sparse_matrix, f_2.sparse_matrix, pin)
        __cnfg[7] = wsum(f_1.dense_matrix, f_2.dense_matrix, pin)
        push!(f_, STLSQresult(__cnfg...))
    end
    # J_ = []; while true try J_ = jacobian.(Function, f_); break; catch; print("."); end end
    data = DataFrame(solve(f_, [eps(), .05853, .47898], dt, 0:dt:150, Dtree), last(vrbl));
    data = data[(nrow(data) ÷ 5):end, :]

    λ = lyapunov_exponent(data[:, last(vrbl)], J_, Dtree, dr.bp)
    dr[names(schedules)[3:end]] .= λ
    
    vrtc[dr.idx] = bifurcation_(sysname, data, dr.bp)
    hrzn[dr.idx] = repeat([dr.bp], length(vrtc[dr.idx]))
    idcs[dr.idx] = repeat([dr.idx], length(vrtc[dr.idx]))

    try
        CSV.write("...$(device)ing $(sysname)_lyapunov_rcvd.csv", schedules, bom = true)
        bfcn = DataFrame(idcs = vcat(idcs...), hrzn = vcat(hrzn...), vrtc = vcat(vrtc...))
        CSV.write("...$(device)ing $(sysname)_bfcn_rcvd.csv", bfcn, bom = true)
    catch e
        @error "\n$(now()): $e in $(dr.idx)"
    end
end
CSV.write("!$(device) $(sysname)_lyapunov_rcvd.csv", schedules, bom = true)
bfcn = DataFrame(idcs = vcat(idcs...), hrzn = vcat(hrzn...), vrtc = vcat(vrtc...))
CSV.write("!$(device) $(sysname)_bfcn_rcvd.csv", bfcn, bom = true)

# ##########################################################################
# #                                                                        #
# #                             Gear system                                #
# #                                                                        #
# ##########################################################################
# sysname = "gear" # Fe ∈ [-0.2, 0.2]
# vrbl = [:dx, :dv, :dΩ, :dθ], [:x, :v, :Ω, :θ]
# cnfg = (; N = 1, f_ = [cos], C = 2,  λ = 1e-4)
# dt = 1e-2; tspan = [0, 1000]; θ = 1e-20; dos = 1

# _idx = 1:801
# tasks = Dict(device => _idx, "chaos1" => 1:650, "chaos2" => 651:1300, "chaos3" => 1301:length(_idx))
# schedules = DataFrame(idx = _idx, bp = LinRange(-0.2, 0.2, length(_idx)))[tasks[device], :]
# for k in eachindex(last(vrbl)) schedules[!, "λ$k"] .= 0.0 end
# bp1, bp2 = sort([-0.2, -0.19]);

# data_ = [factory_(sysname)(DataFrame, bp1; tspan, dt),
#          factory_(sysname)(DataFrame, bp2; tspan, dt)]
# for data in data_ add_subsystem!(data, vrbl, cnfg; θ, dos) end

# f__ = [[SINDy(df, vrbl...; cnfg...) for df in groupby(data, :subsystem)] for data in data_]
# _cnfg_ = [getproperty.(Ref(f_), f_ |> propertynames) for f_ in f__[1]];
# Dtree_ = [dryad(data, last(vrbl)) for data in data_] # print_tree.(Dtree_)
# # df_Dtree_ = dt2df.(Dtree_);

# # M0 = Matrix(_df_Dtree[:, 1:(end-1)])
# @info "$(now()): Preprocess for $(sysname) done!"

# idcs = [Int64[] for _ in _idx]
# hrzn = [Float64[] for _ in _idx]
# vrtc = [Float64[] for _ in _idx]
# @showprogress @threads for dr = eachrow(schedules)
#     pin = (dr.bp - bp1) / (bp2 - bp1)

#     Dtree = deepcopy(Dtree_[1])
#     J_ = [
#         vv -> [0 1 0 0;                     0 -0.12                      -dr.bp*sin(vv[3]) -dr.bp*sin(vv[4]); 0 0 0 0; 0 0 0 0],
#         vv -> [0 1 0 0; -(1 + 0.06cos(vv[3])) -0.12 (dr.bp + 0.06 + 0.06*vv[1])*sin(vv[3]) -dr.bp*sin(vv[4]); 0 0 0 0; 0 0 0 0],
#         vv -> [0 1 0 0; -(1 + 0.06cos(vv[3])) -0.12 (dr.bp - 0.06 + 0.06*vv[1])*sin(vv[3]) -dr.bp*sin(vv[4]); 0 0 0 0; 0 0 0 0],
#     ]

#     f_ = []
#     for (f_1, f_2, _cnfg) = collect(zip(f__[1], f__[2], _cnfg_))
#         __cnfg = collect(_cnfg)
#         __cnfg[5] = wsum(f_1.sparse_matrix, f_2.sparse_matrix, pin)
#         __cnfg[7] = wsum(f_1.dense_matrix, f_2.dense_matrix, pin)
#         push!(f_, STLSQresult(__cnfg...))
#     end

#     data = DataFrame(solve(f_, [0.1, 0.1, 0.1, eps()], dt, 0:dt:last(1500), Dtree), last(vrbl))
#     data = data[(nrow(data) ÷ 3):end, :]

#     λ = lyapunov_exponent(data[:, last(vrbl)], J_, Dtree, dr.bp, T = last(tspan))
#     dr[names(schedules)[3:end]] .= λ
    
#     vrtc[dr.idx] = bifurcation_(sysname, data, dr.bp)
#     hrzn[dr.idx] = repeat([dr.bp], length(vrtc[dr.idx]))
#     idcs[dr.idx] = repeat([dr.idx], length(vrtc[dr.idx]))

#     try
#         CSV.write("...$(device)ing $(sysname)_lyapunov_rcvd.csv", schedules, bom = true)
#         bfcn = DataFrame(idcs = vcat(idcs...), hrzn = vcat(hrzn...), vrtc = vcat(vrtc...))
#         CSV.write("...$(device)ing $(sysname)_bfcn_rcvd.csv", bfcn, bom = true)
#     catch e
#         @error "\n$(now()): $e in $(dr.idx)"
#     end
# end
# CSV.write("!$(device) $(sysname)_lyapunov_rcvd.csv", schedules, bom = true)
# bfcn = DataFrame(idcs = vcat(idcs...), hrzn = vcat(hrzn...), vrtc = vcat(vrtc...))
# CSV.write("!$(device) $(sysname)_bfcn_rcvd.csv", bfcn, bom = true)

##########################################################################
#                                                                        #
#                           Hindmarsh-Rose model                         #
#                                                                        #
##########################################################################
sysname = "hrnm" # l ∈ [0, 1]
vrbl = [:dt, :dx, :dy, :dz], [:t, :x, :y, :z]
cnfg = (; N = 3, f_ = [cos])
dt = 1e-3; θ = 1e-10; tspan = [0, 1000]

_idx = 1:1001
tasks = Dict(device => _idx, "chaos1" => 1:650, "chaos2" => 651:1300, "chaos3" => 1301:length(_idx))
schedules = DataFrame(idx = _idx, bp = LinRange(0, 1, length(_idx)))[tasks[device], :]
for k in eachindex(last(vrbl)) schedules[!, "λ$k"] .= 0.0 end
bp1, bp2 = sort([0.01, 0.02]);

data_ = [factory_(sysname)(DataFrame, bp1; tspan, dt),
         factory_(sysname)(DataFrame, bp2; tspan, dt)]
for data in data_ add_subsystem!(data, vrbl, cnfg; θ) end

f__ = [[SINDy(df, vrbl...; cnfg...) for df in groupby(data, :subsystem)] for data in data_]
_cnfg_ = [getproperty.(Ref(f_), f_ |> propertynames) for f_ in f__[1]];
Dtree_ = [dryad(data, last(vrbl)) for data in data_] # print_tree.(Dtree_)

# M0 = Matrix(_df_Dtree[:, 1:(end-1)])
@info "$(now()): Preprocess for $(sysname) done!"

idcs = [Int64[] for _ in _idx]
hrzn = [Float64[] for _ in _idx]
vrtc = [Float64[] for _ in _idx]
@showprogress @threads for dr = eachrow(schedules)
    pin = (dr.bp - bp1) / (bp2 - bp1)

    Dtree = deepcopy(Dtree_[1])
    J_ = [
        vv -> [0 0 0 0; (-dr.bp*sin(vv[1])) (6vv[2] + 0.9vv[4] - 3*(vv[2]^2)) 1 (0.9vv[2]); 0 (-10vv[2]) -1 0; 0 0.8 0 -0.1],
        vv -> [0 0 0 0; (-dr.bp*sin(vv[1])) (6vv[2] + 0.9vv[4] - 3*(vv[2]^2)) 1 (0.9vv[2]); 0 (-10vv[2]) -1 0; 0 0.8 0 -0.1],
        vv -> [0 0 0 0; (-dr.bp*sin(vv[1])) (6vv[2] + 0.9vv[4] - 3*(vv[2]^2)) 1 (0.9vv[2]); 0 (-10vv[2]) -1 0; 0 0.8 0 -0.1],
    ]

    f_ = []
    for (f_1, f_2, _cnfg) = collect(zip(f__[1], f__[2], _cnfg_))
        __cnfg = collect(_cnfg)
        __cnfg[5] = wsum(f_1.sparse_matrix, f_2.sparse_matrix, pin)
        __cnfg[7] = wsum(f_1.dense_matrix, f_2.dense_matrix, pin)
        push!(f_, STLSQresult(__cnfg...))
    end

    data = DataFrame(solve(f_, [eps(), eps(), eps(), 0.1], dt, 0:dt:last(tspan), Dtree), last(vrbl));
    data = data[(nrow(data) ÷ 5):end, :]

    λ = lyapunov_exponent(data[:, last(vrbl)], J_, Dtree, dr.bp)
    dr[names(schedules)[3:end]] .= λ
    
    vrtc[dr.idx] = bifurcation_(sysname, data, dr.bp)
    hrzn[dr.idx] = repeat([dr.bp], length(vrtc[dr.idx]))
    idcs[dr.idx] = repeat([dr.idx], length(vrtc[dr.idx]))

    try
        CSV.write("...$(device)ing $(sysname)_lyapunov_rcvd.csv", schedules, bom = true)
        bfcn = DataFrame(idcs = vcat(idcs...), hrzn = vcat(hrzn...), vrtc = vcat(vrtc...))
        CSV.write("...$(device)ing $(sysname)_bfcn_rcvd.csv", bfcn, bom = true)
    catch e
        @error "\n$(now()): $e in $(dr.idx)"
    end
end
CSV.write("!$(device) $(sysname)_lyapunov_rcvd.csv", schedules, bom = true)
bfcn = DataFrame(idcs = vcat(idcs...), hrzn = vcat(hrzn...), vrtc = vcat(vrtc...))
CSV.write("!$(device) $(sysname)_bfcn_rcvd.csv", bfcn, bom = true)


# plot(data.u[1:100:end])

# f__[2][1].sparse_matrix.nzval |> collect

# Dtree_[2].node.right.featval
# Dtree_[2].node.featval
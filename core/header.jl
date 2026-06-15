using Base.Threads: @threads, nthreads # Base.Threads.nthreads()
import Pkg
packages = [:Combinatorics, :LinearAlgebra, :SparseArrays, :DataFrames,
            :PrettyTables, :Symbolics, :CSV, :DecisionTree, :Random,
            :StatsBase, :Dates, :ProgressMeter, :Plots, :LaTeXStrings,
            :Colors, :ColorSchemes, :JLD2, :Graphs, :DifferentialEquations
            ] .|> string
try
    @time "All packages load" eval(Meta.parse("using $(join(packages, ", "))"))
catch e
    # required = setdiff(packages, getproperty.(values(Pkg.dependencies()), :name))
    required = setdiff(packages, keys(Pkg.installed()))
    print("\r\x1b[F")
    if required |> !isempty
        @info "Installing required packages: $(required)"
        Pkg.add.(required)
    end
    @time "All packages load" eval(Meta.parse("using $(join(packages, ", "))"))
end
default(legend = :none, framestyle = :box, dpi = 180)
mm = Plots.mm
cm = Plots.cm
# RK4 = DifferentialEquations.RK4

device = gethostname() # In julia v1.11, it could be replaced by `Sys.username()`
# if occursin(lowercase(pwd()), lowercase(@__DIR__))
#     if Sys.iswindows()
#         if device != "Sickbook"
#             cd("G:/DDM")
#         end
#         # device = ENV["COMPUTERNAME"];
#     elseif Sys.islinux()
#         cd("/home/$(ENV["LOGNAME"])/g/DDM")
#         # device = ENV["LOGNAME"];
#     end
# end

include("../core/utils.jl")
include("../core/DDM.jl")

@info "$(now()) - @$device $(nthreads()) threads"
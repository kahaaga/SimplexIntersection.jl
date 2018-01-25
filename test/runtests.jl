using SimplexIntersection
using Base.Test
using BenchmarkTools

include("speed_test.jl")
@time speed_test(3,1)
#@benchmark speed_test(3)
@time include("test_NDintersection.jl")

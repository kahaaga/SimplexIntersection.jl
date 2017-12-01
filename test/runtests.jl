using SimplexIntersection
using Base.Test

# write your own tests here
tic()
println("Test 1")
@time include("test_simplexintersection.jl")
#println("Test 2")
#@time @test include("test2.jl")
toc()

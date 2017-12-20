using SimplexIntersection
using Base.Test

# write your own tests here
tic()
@time include("test_simplexintersection.jl")
#@time include("test_NDintersection.jl")
toc()

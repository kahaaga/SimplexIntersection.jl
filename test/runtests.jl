using SimplexIntersection
using Base.Test

#!/usr/bin/env julia
# Run tests

println("Testing simplex intersection against MATLAB results")
@time @test include("test_simplexintersection.jl")

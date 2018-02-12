using SimplexIntersection
using Base.Test

include("speed_test.jl")
speed_test_nontrivial(3, 1)
speed_test_nontrivial(4, 1)
speed_test_nontrivial(5, 1)

speed_test_sharing(3, 1)
speed_test_sharing(4, 1)
speed_test_sharing(5, 1)
#
nreps = 50
t3D_nontrivial = @elapsed speed_test_nontrivial(3, nreps)
println("Speed test:", nreps, " nontrivial intersections in 3D ... ", t3D_nontrivial/nreps, " s per intersection")
#
# t4D_nontrivial = @elapsed speed_test_nontrivial(4, nreps)
#
# println("Speed test:", nreps, " nontrivial intersections in 4D ... ", t4D_nontrivial/nreps, " s per intersection")

# t5D_nontrivial = @elapsed speed_test_nontrivial(5, nreps)
# println("Speed test:", nreps, " nontrivial intersections in 5D ... ", t5D_nontrivial/nreps, " s per intersection")
# #
# t3D_shared = @elapsed speed_test_sharing(3, nreps)
# println("Speed test:", nreps, " shared intersections in 3D ... ", t3D_shared/nreps, " s per intersection")
#
# t4D_shared = @elapsed speed_test_sharing(4, nreps)
# println("Speed test:", nreps, " shared intersections in 4D ... ", t4D_shared/nreps, " s per intersection")
#
# t5D_shared = @elapsed speed_test_sharing(5, nreps)
# println("Speed test:", nreps, " shared intersections in 5D ... ", t5D_shared/nreps, " s per intersection")

# #
# # include("test_delaunay.jl")
# # delaunaytest(10, 3, 100)
# # delaunay_benchmark_old(10, 3, 1)
# # delaunay_benchmark_new(10, 3, 1)
# #
# # println("\nDELAUNAY TRIANGULATION ...")
# # println("Old delaunay (1000 reps, 10 points in 3D):")
# # @time delaunay_benchmark_old(10, 3, 1000)
# # println("New delaunay (1000 reps, 10 points in 3D):")
# # @time delaunay_benchmark_new(10, 3, 1000)
# #
# # println("Old delaunay (1000 reps, 10 points in 4D):")
# # @time delaunay_benchmark_old(10, 4, 1000)
# # println("New delaunay (1000 reps, 10 points in 4D):")
# # @time delaunay_benchmark_new(10, 4, 1000)

#@time include("test_NDintersection.jl")

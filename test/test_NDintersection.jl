include("nd_test_function.jl")
include("intersectiontest_funs.jl")

reps = 1
k = 2

function init()
    for E = 3:3
        nd_Test(k, E, reps)
    end
end

reps = 10

@testset "nD discrepancy: strictly contained" begin
    init()
    println("nD discrepancy test with ", reps, " reps:")
    #Es = 2:6
    @testset "E = $E" for E in 3:5
        # Trigger once for precompilation

        t1 = time_ns()/10^9
        discrepancies = nd_Test(k, E, reps)
        t2 = time_ns()/10^9
        elapsed =  t2 - t1
        println("E = ", E, " | ", (k^E), " possible intersections.\t| Max discrepancy: ", maximum(discrepancies), " | Per intersection: ", elapsed/reps/(k^E)*10^3, " ms | Total time for one triangulation: ", elapsed/reps*10^3, " ms | total time: ", elapsed, "seconds")
        @test maximum(discrepancies) < 1/10^10
    end
end

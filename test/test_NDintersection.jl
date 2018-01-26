include("nd_test_function.jl")
include("intersectiontest_funs.jl")

reps = 50
function init()
    for i = 2:6
        nd_Test(2,i,1)
    end
end

k = 2

@testset "nD discrepancy: strictly contained" begin
    init()
    println("nD discrepancy test with ", reps, " reps:")
    Es = 2:6
    @testset "E = $i" for i in 1:length(Es)
        # Trigger once for precompilation

        t1 = time_ns()/10^9
        discrepancies = nd_Test(k, Es[i], reps)
        t2 = time_ns()/10^9
        elapsed =  t2 - t1
        println("E = ", Es[i], " | ", (k^Es[i]), " possible intersections.\t| Max discrepancy: ", maximum(discrepancies), " | Per intersection: ", elapsed/reps/(k^Es[i])*10^3, " ms | Total time for one triangulation: ", elapsed/reps*10^3, " ms | total time: ", elapsed, "seconds")
        @test maximum(discrepancies) < 1/10^10
    end
end

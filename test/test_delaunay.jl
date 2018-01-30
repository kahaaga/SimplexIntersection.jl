using SimplexIntersection

function delaunaytest(n::Int, dim::Int, reps::Int)
    for i = 1:reps
        pts = rand(n, dim)
        @test SimplexIntersection.QHull.delaunay_tesselation(pts) ==
        SimplexIntersection.QHull.delaunay_tesselation(pts)
    end
end

function delaunay_benchmark_old(n::Int, dim::Int, reps::Int)
    for i = 1:reps
        pts = rand(n, dim)
        SimplexIntersection.QHull.delaunay_tesselation(pts)
    end
end

function delaunay_benchmark_new(n::Int, dim::Int, reps::Int)
    for i = 1:reps
        pts = rand(n, dim)
        SimplexIntersection.QHull.delaunayn(pts)
    end
end

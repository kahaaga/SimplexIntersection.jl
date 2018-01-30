using SimplexIntersection

function speed_test_nontrivial(dim,N)
    for i = 1:N
        S1,S2 = SimplexIntersection.nontrivially_intersecting_simplices(dim)
        SimplexIntersection.simplexintersection(S1.', S2.')
    end
end


function speed_test_sharing(dim,N)
    for i = 1:N
        S1,S2 = SimplexIntersection.simplices_sharing_vertices(dim)
        SimplexIntersection.simplexintersection(S1.', S2.')
    end
end

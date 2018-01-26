function speed_test(dim,N)

    T = time_ns()/10^9
    for i =1:N
        S1,S2 = SimplexIntersection.nontrivially_intersecting_simplices(dim)

         SimplexIntersection.simplexintersection(S1.', S2.')

    end
    (time_ns()/10^9 - T)/N 
end

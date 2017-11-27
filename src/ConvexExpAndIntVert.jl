function ConvexExpAndIntVert(IntVert,NewPoint,alpha,lambda,ReferenceBoundary,TargetBoundary,Switch,ConvexExpIntVert,n)
    ##
    #OUTCOME:
    #
    # MinimalBoundIntVert:
    # Array with as many rows as IntVert has, and 2n+2 clumns. Each row contains the convex expansion coefficinets of the corresponding intersecting point
    # in terms of the original vertices of the simplices in the order set by SimplexVert1 and SimplexVert2, respectively
    # First the faces of Simplex1 and then those of Simplex2.
    # Recall that IndexVert contains the elements 1:n+1 but ordered such that
    # the corresponding vertices of the simplex contained in the other one,
    # appear first
    ##


    Aux = zeros(1, (2*n)+2)
    Aux[vcat(ReferenceBoundary, TargetBoundary + n + 1)] = hcat(alpha, lambda)
    #println("Aux " , Aux)
    if Switch == 1
        Aux = zeros(1, 2*n + 2)
        Aux[vcat(TargetBoundary, ReferenceBoundary + n + 1)] = hcat(lambda, alpha)
    end

    return Update(IntVert, NewPoint), Update(ConvexExpIntVert, Aux)
end

function ConvexExp(Vert1Inside2, Vert1Outside2)

    # Outcome:
    # beta: is a n x n matrix in which each column contains the beta parameters of the
    # convex expansion of the vertex of Simplex1 with the same index as the column in terms of the vertices of Simplex2

    beta = copy(Vert1Outside2)

    if size(Vert1Inside2, 1) > 1
        beta = [Vert1Inside2 Vert1Outside2]
    end
    beta = transpose(sortrows(transpose(beta)))
    return beta[2:end, 1:end]
end

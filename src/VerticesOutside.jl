
"""
VerticesOutside()

Given the vertices generating the target simplex, the information about which vertices are shared by both simplces, and which vertices of the target simplex are contained in the source simplex, compute the convex expansion of the remaining vertices (the vertices of the target simplex not contained in the source simplex). This has already been done for the vertices of the target simplex lying inside the circumsphere of the source simplex, but not the vertices lying outside the circumsphere.

Input arguments
---------------
SourceSimplex::Array{Float64, 2}            An N-by-N+1 array. Each column is a vertex of the source simplex.
TargetSimplex::Array{Float64, 2}            An N-by-N+1 array. Each column is a vertex of the target simplex.
InfoVerticesTargetInSource::Array{Int-Float64, 2}   The first row contains integers, remaining rows are floats. The first rows contains a permutation of the elements in 'IndexOfNoncommonVerticesOfSourceInsideCircumphereOfTarget' such that the indices correspond to contained vertices appear first (so that the first 'NumVertInNotComm' columns correspond to vertices that are contained in the source simplex). Below each index, we put the convex expansion of the vertex in the target simplex with the same index in terms of the vertices of the source simplex. The sum of the coefficients should be one. For the vertices that are contained, all coefficients in the convex expansion are non-negative. For the remaining vertices that are not contained, at least one coefficient in the convex expansion is negative.      CHANGE THIS TO RETURN INDICES AND THE CONVEX EXPANSIONS SEPARATELY
Ncomm::Int                                  Number of vertices that are shared by both simplices.
n::Int                                      The dimension of the space. Redundant. Loose it!
IndexCommVert::Array{Int, 1}                Row vector. The indices of the vertices of the target simplex that are shared by both simplices.

Returns
-------
Vert1Outside2::Array{Int-Float64, 2}         A (n+2)-by-M matrix, where M = n + 1 - 'numofTargetinSourceNotCom'. First row contains integers corresponding to the indices of the vertices of the target simplex lying outside the source simplex. The remaining matrix contain the corresponding convex expansions in terms of the vertices of the source simplex. In each column, at least one of the coefficients in the convex expansion must be negative (otherwise the point would be contained).

"""

function VerticesOutside(Y,X,InfoVertices1in2,numof1in2NotCom,All1in2,O2,tolerance,n,Ncomm,IndexCommVert)
    # Y vertices of the holding simplex
    X = copy(X)
    Y = copy(Y)

    indices1in2 =  round.(Int64, InfoVertices1in2[1, :]) # First row of infovertices1in2

    if (Ncomm == 0)
        if (numof1in2NotCom > 0)
            # Vertices of S1 outside S2
            if (All1in2 == 1)
                # All the checked vertices belong to S2
                # There must be outside the circumsphere
                IndexVert1Outside2 = complementary(indices1in2, n + 1) # CONVERT TO INT
                Vert1Outside2 = InfoVerticesOutside(Y, X, IndexVert1Outside2, O2, tolerance)
            else
                #Some of the checked vertices do not belong to S2
                if size(InfoVertices1in2, 2) == n + 1
                    # No vertex outside the circumsphere
                    Vert1Outside2 = InfoVertices1in2[:, (numof1in2NotCom + 1):(n + 1)]
                else
                    # Some vertices outside the circumsphere
                    IndexVert1Outside2 = complementary(indices1in2, n + 1)
                    Vert1Outside2 = InfoVerticesOutside(Y, X, IndexVert1Outside2, O2, tolerance)
                    Vert1Outside2 = [Vert1Outside2 InfoVertices1in2[:, (numof1in2NotCom + 1):size(InfoVertices1in2, 2)]]
                end
            end
        else # numof1in2NotCom=0
            # Either all checked vertices lie outside the simplex or there was no vertex to check
            if (size(InfoVertices1in2, 1) == 1)
                # There was no vertex to check in the first place (no vertices in the circumsphere)
                Vert1Outside2 = InfoVerticesOutside(Y, X, collect(1:(n + 1)), O2, tolerance)
            else
                # Vertices of S1 outside S2 (all of the checked vertices + ...)
                Vert1Outside2 = InfoVertices1in2[:, :]
                if (size(InfoVertices1in2, 2) < n + 1)
                    # Some vertices outside the circumsphere
                    Aux = complementary(indices1in2, n + 1)
                    Aux = InfoVerticesOutside(Y, X, Aux, O2, tolerance)
                    Vert1Outside2 = [Vert1Outside2 Aux]
                end
            end

        end
    else # Ncomm>0
        if (numof1in2NotCom > 0)
            if (All1in2 == 1)
                # All of the checked vertices belong to S2
                # There must be outside the circumsphere
                Vert1Outside2 = complementary([IndexCommVert; indices1in2], n + 1)
                Vert1Outside2 = InfoVerticesOutside(Y, X, Vert1Outside2, O2, tolerance)
            else
                # Some of the checked vertices do not belong to S2
                if size(InfoVertices1in2, 2) + Ncomm == n + 1
                    # No vertex outside the circumsphere
                    Vert1Outside2 = InfoVertices1in2[:, (numof1in2NotCom + 1):size(InfoVertices1in2, 2)]
                else
                    # Some vertices outside the circumsphere
                    IndexVert1Outside2 = complementary([IndexCommVert indices1in2'], n + 1)
                    Vert1Outside2 = InfoVerticesOutside(Y, X, IndexVert1Outside2, O2, tolerance)
                    Vert1Outside2 = [Vert1Outside2 InfoVertices1in2[:, (numof1in2NotCom + 1):size(InfoVertices1in2, 2)]]
                end
            end
        else #numof1in2NotCom=0
            # Either all checked vertices lie outside the simplex or there was no vertex to check
            # Vertices of S1 outside S2
            if size(InfoVertices1in2, 1) == 1
                # There was no vertex to check in the first place (no vertices in the circumsphere)
                IndexVert1Outside2 = complementary(IndexCommVert, n + 1)
                Vert1Outside2 = InfoVerticesOutside(Y, X, IndexVert1Outside2, O2, tolerance)
            else
                if Ncomm + size(InfoVertices1in2, 2) == n + 1
                    # No vertex outside the circumsphere
                    Vert1Outside2 = InfoVertices1in2[:, :]
                else
                    # Some vertices outside the circumsphere
                    IndexVert1Outside2 = complementary([IndexCommVert; indices1in2], n + 1)
                    Vert1Outside2 = InfoVerticesOutside(Y, X, IndexVert1Outside2, O2, tolerance)
                    Vert1Outside2 = [Vert1Outside2 InfoVertices1in2[:, :]]
                end
            end
        end
    end
    return Vert1Outside2
end

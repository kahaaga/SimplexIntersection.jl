
"""
VerticesInside()

Given vertices of the target simplex lying in the circumsphere of the source simplex, select those that actually lie inside the source simplex.

Input arguments
---------------
sorted_vertex_indices::Array{Int, 1}  Indices of target simplex vertices not shared with the source simplex, but lying inside the circumsphere of the source simplex.

convex_expansion::Array{Float64, 2}   Convex expansions of target simplex vertices in terms of source simplex vertices.
Columns are ordered according to 'sorted_vertex_indices'. Each column contains the
convex expansion of a target simplex vertex.  Coefficients should sum to 1 in each row.
For the vertices that are contained, all coefficients in the convex expansion are
non-negative. Remaining target simplex vertices (not contained in the source simplex)
have at least one negative coefficient in its convex expansion.

n_targvertices_in_src_circum::Int     The number of target simplex vertices that lie inside the circumsphere of the source simplex,
but are not vertices of both simplices (non-shared vertices).

numofTargetinSourceNotCom::Int        The number of target simplex vertices that lie inside the circumsphere of the source simplex,
but are not vertices of both simplices (non-shared vertices).

inds_sharedvertices_target::Array{Int, 1}   Row vector. Contains the indices of vertices in the target simplex that are shared between both simplices.
inds_sharedvertices_source:Array{Int, 1}    Row vector. Contains the indices of vertices in the source simplex that are shared between both simplices.
n_sharedvertices::Int                       Number of vertices that are shared by both simplices.
dim::Int                                    The dimension of the space. Redundant. Loose it!

------
InfoVerticesTargetInSource::Array{Int-Float64}   The first row contains integers, remaining rows are floats. The first rows contains a permutation of the elements in 'IndexOfNoncommonVerticesOfSourceInsideCircumphereOfTarget' such that the indices correspond to contained vertices appear first (so that the first 'NumVertInNotComm' columns correspond to vertices that are contained in the source simplex). Below each index, we put the convex expansion of the vertex in the target simplex with the same index in terms of the vertices of the source simplex. The sum of the coefficients should be one. For the vertices that are contained, all coefficients in the convex expansion are non-negative. For the remaining vertices that are not contained, at least one coefficient in the convex expansion is negative.      CHANGE THIS TO RETURN INDICES AND THE CONVEX EXPANSIONS SEPARATELY
IndexComVertTarget::Array{Int, 1}           Row vector. Contains the indices of vertices in the target simplex that are shared between both simplices.
IndexComVertSource::Array{Int, 1}           Row vector. Contains the indices of vertices in the source simplex that are shared between both simplices.
Ncomm::Int                                  Number of vertices that are shared by both simplices.
n::Int                                      The dimension of the space. Redundant. Loose it!
-----

Returns
-------
sorted_vertex_indices::Array{Int, 1}  Row vector containing a permutation of the elements in 'IndexOfNoncommonVerticesOfSourceInsideCircumphereOfTarget' such that the indices correspond to vertices contained in the source simplex appear first (so that the first 'NumVertInNotComm' columns correspond to vertices that are contained in the source simplex).
convex_expansion::Array{Float64, 2}   Convex expansion of the target simplex vertices in terms of the vertices of the source simplex. Each column contains the convex expansion for a given point in the target simplex, and the columns are ordered according to 'sorted_vertex_indices'.  Coefficients should sum to 1 in each row. For the vertices that are contained, all coefficients in the convex expansion are non-negative. Remaining target simplex vertices (not contained in the source simplex) have at least one negative coefficient in its convex expansion.      CHANGE THIS TO RETURN INDICES AND THE CONVEX EXPANSIONS SEPARATELY


----
Vert1Inside2::Array{Int-Float64, 2}         A (n+2)-by-M matrix. M depends on how many vertices lie outside or not. First row is integers corresponding to the indices of the vertices of the target simplex lying inside the source simplex. The remaining matrix contain the corresponding convex expansions in terms of the vertices of the source simplex.
----
"""
function VerticesInside(InfoVertices1in2,numof1in2NotCom,IndexComVert1,IndexComVert2,n,Ncomm)

    # Y vertices of the holding simplex

    # X vertices of the simplex which location we want to determine

    # InfoVertices1in2  information on the location of the vertices inside the circumsphere

    # IndexCommVert1   indices in (X) of the shared vertices

    # IndexCommVert2   indices in (Y) of the shared vertices

    # O2  orientation of the holding simplex

    # n  dimension

    # Ncomm  number of shared vertices

    # numof1in2 number of vertices of S1 in the circumsphere of S2

    # numof1in2NotCom  number of vertices of S1 in the circumsphere of S2 that are not common

    # OUTCOME

    # VerticesInside is a (n+2)xM matrix
    # the first row are the indices of
    # the contained vertices starting with the shared ones (if any), in the orther specified by IndexComVert1
    # below each index the (n+1)-column vector
    # contains the coefficients of the convex expansion of the corresponding vertex in the vertices of the other simplex


    # If the result is empty, the program just returns 0

    Vert1Inside2 = [0]

    if Ncomm == 0
        if numof1in2NotCom > 0
            Vert1Inside2 = InfoVertices1in2[:, 1:numof1in2NotCom]
        end
    else # Ncomm>0
        Vert1Inside2 = zeros(n + 2, Ncomm)
        Vert1Inside2[1, :] = IndexComVert1
        for i = 1:Ncomm
            Vert1Inside2[IndexComVert2[i]+1, i] = 1
        end
        if numof1in2NotCom > 0
            # There are other indices in the circumsphere of S2 not common
            # and yet contained
            Vert1Inside2 = [Vert1Inside2 InfoVertices1in2[:, 1:numof1in2NotCom]]
        end
    end
    return Vert1Inside2
end

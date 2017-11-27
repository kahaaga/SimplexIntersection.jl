"""
ContainedVertices()

Given a source and a target simplex, and the indices of the source simplex vertices that are not common to both simplices, but still lie inside the circumsphere of the target simplex.
The function finds the vertices actually belonging to the source simplex.


Input arguments
---------------
SourceSimplex::Array{Float64, 2}      An N-by-N+1 matrix containing vertices of the source simplex.
TargetSimplex::Array{Float64, 2}      An N-by-N+1 matrix containing vertices of the target simplex.
IndexOfNoncommonVerticesOfSourceInsideCircumphereOfTarget::Array{Int, 1}    Row vector containing the indices of the source simplex vertices that are not shared between vertices, but lie inside the circumsphere of the target simplex.
SimplexOrientation::Float64           The orientation of the source simplex (which is the determinant of the matrix of simplex vertices, but with a row of ones appended to the top).
tolerance::Float64                    Any convex expansion coefficient less than 'tolerance' in absolute value is set to zero. Used to decide if a point belongs to a simplex.

Returns
-------
NumVertInNotComm::Int                 The number of vertices generating the target simplex that are not shared between both simplces, but lie inside the source simplex.
AllVerticesAreContained::Int          1 if all the checked non-shared vertices of the target simplex are inside the source simplex 0 otherwise.
InfoVertices::Array{Float64, 2}       The first row contains integers, remaining rows are floats. The first rows contains a permutation of the elements in 'IndexOfNoncommonVerticesOfSourceInsideCircumphereOfTarget' such that the indices correspond to contained vertices appear first (so that the first 'NumVertInNotComm' columns correspond to vertices that are contained in the source simplex). Below each index, we put the convex expansion of the vertex in the target simplex with the same index in terms of the vertices of the source simplex. The sum of the coefficients should be one. For the vertices that are contained, all coefficients in the convex expansion are non-negative. For the remaining vertices that are not contained, at least one coefficient in the convex expansion is negative.      CHANGE THIS TO RETURN INDICES AND THE CONVEX EXPANSIONS SEPARATELY

sorted_vertex_indices::Array{Int, 1}  Row vector containing a permutation of the elements in 'IndexOfNoncommonVerticesOfSourceInsideCircumphereOfTarget' such that the indices correspond to vertices contained in the source simplex appear first (so that the first 'NumVertInNotComm' columns correspond to vertices that are contained in the source simplex).
convex_expansion::Array{Float64, 2}   Convex expansion of the target simplex vertices in terms of the vertices of the source simplex. Each column contains the convex expansion for a given point in the target simplex, and the columns are ordered according to 'sorted_vertex_indices'.  Coefficients should sum to 1 in each row. For the vertices that are contained, all coefficients in the convex expansion are non-negative. Remaining target simplex vertices (not contained in the source simplex) have at least one negative coefficient in its convex expansion.      CHANGE THIS TO RETURN INDICES AND THE CONVEX EXPANSIONS SEPARATELY
"""

function ContainedVertices(
  SimplexVert,
  SimplexVertCheck,
  PreIndexContainedVertices,
  SimplexOrientation,
  tolerance)

  n = size(SimplexVert, 1)

  NumVertInNotComm =  0
  All = 0
  InfoVertices= [0]

  if PreIndexContainedVertices[1] > 0
    Points = copy(SimplexVertCheck[:, PreIndexContainedVertices])

    np = size(Points, 2)

    InfoVertices = zeros(n + 2, np)

    InfoVertices[1, :] = PreIndexContainedVertices

    Contained = collect(1:np)


    for i = 1:np
      firstentry=1

      for j = 1:(n + 1)
        Aux = copy(SimplexVert)
        Aux[:, j] = Points[:, i]
        beta = det([ones(1, n + 1); Aux]) / SimplexOrientation

        if abs(beta) <= tolerance
          beta = 0
        end

        InfoVertices[j + 1, i] = beta

        if beta < 0
          # First time a negative beta is found
          if firstentry == 1
            Contained[i] = 0
            firstentry = -1
          end
        end
      end
    end

    oldnp = np
    Cont = find(Contained) # Contained indices
    np = size(Cont, 1)

    # all vertices are contained
    if np == oldnp
      NumVertInNotComm = np
      All = 1

    # some of the vertices are not contained
    elseif np > 0
      NotCont =  complementary(Cont, oldnp) #Not contained indices
      InfoVertices = InfoVertices[:, vcat(Cont, NotCont)]
      NumVertInNotComm = np

      All = 0
    end
  end

  return NumVertInNotComm::Int64, All::Int64, InfoVertices
end

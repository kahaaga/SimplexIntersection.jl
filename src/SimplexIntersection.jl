module SimplexIntersection
export simplexintersection

using SimplexSplitting

include("CommonVertices.jl")
include("NonCommonVertices.jl")
include("ContainedVertices.jl")
include("ShareFace_nD.jl")
include("geometry.jl")
include("heaviside.jl")
include("VerticesInside.jl")
include("VerticesOutside.jl")
include("Circumsphere.jl")
include("InsideCircum.jl")
include("complementary.jl")
include("InfoVerticesOutside.jl")
include("IntersectingBoundaries.jl")
include("ConvexExpansion.jl")
include("Binary.jl")
include("BoundaryxBoundary.jl")
include("QR.jl")
include("qhull.jl")
include("UpdateNonZeroSearchingIndex.jl")
include("MinimalBoundaries.jl")
include("ConvexExpAndIntVert.jl")
include("Update.jl")
include("IntersectingOfPolytopeVertices.jl")
include("VolumeBT.jl")
include("TriangulationPolytopeFaces.jl")
include("NullSpace.jl")
include("TriangulationNonSimplicialFaces.jl")
include("SimplexChecks.jl")

"""
    SimplexIntersection()

Computes the volume of intersection between two n-dimensional simplices
by boundary triangulation.

# Input arguments
`S1::Array{Float64, 2}` Simplex 1 represented a matrix of dimension nx(n+1), where each
  column vector is a vertex.
`S2::Array{Float64, 2}` Simplex 2 represented a matrix of dimension nx(n+1), where each
  column vector is a vertex.
`tolerance::Float64` A threshold to determine whether the convex expansion parameters
  of two simplices are different enough that vertices are distinct from each other.
`what::Bool` What to return. "volume" return only the volume (a Float64), "vertices"
  returns an Array{Float64, 2} where each column is a vertex of the intersecting polytope,
  and "both" returns a tuple (volume, vertices).
"""
function simplexintersection(S1::Array{Float64, 2}, S2::Array{Float64, 2};
                              tolerance::Float64 = 1/10^10, what = "volume")

  # Dimension
  n = size(S1, 1)

  # Centroid and radii
  c1 = Circumsphere(S1)[2:n+1]
  c2 = Circumsphere(S2)[2:n+1]
  r1 = Circumsphere(S1)[1]
  r2 = Circumsphere(S2)[1]

  # Orientation of simplices
  orientation_S1 = det([ones(1, n + 1); S1])
  orientation_S2 = det([ones(1, n + 1); S2])

  # Set volume to zero and intersection vertices to an empty array initially
  # Change only if there is intersection
  IntVol = 0.0

  # -------------------------------------
  # Simplices intersect in some way
  # -------------------------------------

  # If the (distance between centroids)^2-(sum of radii)^2 < 0,
  # then the simplices intersect in some way.
  dist_difference = (transpose(c1 - c2) * (c1 - c2) - (r1 + r2)^2)[1]

  if (dist_difference < 0)
    # Find the number of points of each simplex contained within the
    # circumsphere of the other simplex
    Index1in2, numof1in2 = InsideCircum(S1, r2, c2, n)
    Index2in1, numof2in1 = InsideCircum(S2, r1, c1, n)

    # ------------------------------------------------------------------------
    # At least one circumsphere contains vertices of the other simplex
    # ------------------------------------------------------------------------
    if numof1in2 * numof2in1 > 0
      commonvertices = CommonVertices(copy(S1), copy(S2), Index1in2, Index2in1,
                                      numof1in2, numof2in1, n)

      Ncomm = commonvertices[1]             # Number of common vertices
      InternalComIndex1 = commonvertices[2] # Indices in the vector Index1in2 of the indices corresponding to the contained vertices
      InternalComIndex2 = commonvertices[3] # Indices in the vector Index2in1 of the indices cooresponding to the contained vertices
      IndexComVert1 = commonvertices[4]     # Indices of common vertices for simplex S1
      IndexComVert2 = commonvertices[5]     # Indices of common vertices for simplex S1

      # -------------------------------------
      # The simplices coincide
      # -------------------------------------
      if Ncomm == n + 1
        IntVert = S1 # intersecting vertices are just the vertices of one of the simplices
        IntVol = abs(orientation_S1)

      # -------------------------------------
      # The simplices share a face
      # -------------------------------------
      elseif Ncomm == n
        IntVol = ShareFace_nD(S1, S2, IndexComVert1, IndexComVert2,
                              orientation_S1, orientation_S2, tolerance)

      # -------------------------------------
      # The simplices might have nontrivial intersection
      # -------------------------------------
      else
        # Indices of the non common vertices of S1 in the circumsphere of S2
        NonComVert1in2 = NonCommonVertices(InternalComIndex1, Index1in2, numof1in2, Ncomm)

        # Indices of the non common vertices of S2 in the circumsphere of S1
        NonComVert2in1 = NonCommonVertices(InternalComIndex2, Index2in1, numof2in1, Ncomm)

        numof1in2NotCom, All1in2, InfoVertices1in2 = ContainedVertices(S2, S1,
                                                NonComVert1in2, orientation_S2, tolerance)
        numof2in1NotCom, All2in1, InfoVertices2in1 = ContainedVertices(S1, S2,
                                                NonComVert2in1, orientation_S1, tolerance)

        # ----------------------
        # S1 is contained in S2
        # ----------------------
        if (numof1in2NotCom + Ncomm == n + 1)
          IntVert = S1
          IntVol = abs(orientation_S1)

        # ----------------------
        # S2 is contained in S1
        # ----------------------
        elseif (numof2in1NotCom + Ncomm == n + 1)
          IntVert = S2
          IntVol = abs(orientation_S2)

        # -----------------------------
        # Intersection is more complex
        # -----------------------------
        else
          Vert1Inside2 = VerticesInside(InfoVertices1in2, numof1in2NotCom, IndexComVert1,
                                        IndexComVert2, n, Ncomm)
          Vert2Inside1 = VerticesInside(InfoVertices2in1, numof2in1NotCom, IndexComVert2,
                                        IndexComVert1, n, Ncomm)
          Vert1Outside2 = VerticesOutside(S2, S1, InfoVertices1in2, numof1in2NotCom,
                              All1in2, orientation_S2, tolerance, n, Ncomm, IndexComVert1)
          Vert2Outside1 = VerticesOutside(S1, S2, InfoVertices2in1, numof2in1NotCom,
                              All2in1, orientation_S1, tolerance, n, Ncomm, IndexComVert2)

          IntVert, ConvexExpIntVert = IntersectingBoundaries(S1, S2,
                                                              Vert1Inside2, Vert2Inside1,
                                                              Vert1Outside2, Vert2Outside1,
                                                              Ncomm, n, tolerance)
          dim = size(IntVert, 2)

          # Are there enough points to form a polytope?
          if dim > 1
            IntVert, ConvexExpIntVert = IntersectionPolytopeVertices(S1, S2,
                                          IntVert, ConvexExpIntVert,
                                          Vert1Inside2, Vert2Inside1,
                                          numof1in2NotCom, numof2in1NotCom,
                                          Ncomm, n)
            IntVol = VolumeBT(IntVert, ConvexExpIntVert, n)
          end
        end
      end
    end
  end

  if what == "volume"
    return IntVol
  elseif what == "vertices"
    return IntVert
  elseif what == "both"
    return IntVol, IntVert
  end
end

end # module

module SimplexIntersection

using SimplexSplitting

include("Circumsphere.jl")
include("barycentric-coordinates.jl")
include("ShareFace_nD.jl")
include("geometry.jl")
include("heaviside.jl")
include("complementary.jl")
include("InfoVerticesOutside.jl")
include("intersection-of-boundaries.jl")
include("some-vertex-in-circumsphere.jl")
include("ConvexExpansion.jl")
include("Binary.jl")
include("BoundaryxBoundary.jl")
include("sharing-a-face.jl")
include("QR.jl")
include("UpdateNonZeroSearchingIndex.jl")
include("MinimalBoundaries.jl")
include("ConvexExpAndIntVert.jl")
include("Update.jl")
include("IntersectingOfPolytopeVertices.jl")
include("polytope-generating-vertices.jl")
include("VolumeBT.jl")
include("TriangulationPolytopeFaces.jl")
include("NullSpace.jl")
include("TriangulationNonSimplicialFaces.jl")
include("SimplexChecks.jl")
include("volume-computation.jl")
include("shared-vertices.jl")
include("simplexoperations.jl")
include("nd_test.jl")
"""
    SimplexIntersection()

Computes the volume of intersection between two n-dimensional simplices
by boundary triangulation.

Input arguments
---------------
S1::Array{Float64, 2} Simplex 1 represented a matrix of dimension nx(n+1), where each column vector is a vertex.
S2::Array{Float64, 2} Simplex 2 represented a matrix of dimension nx(n+1), where each column vector is a vertex.

Returns
-------
"""

function simplexintersection(S1::Array{Float64, 2}, S2::Array{Float64, 2}; tolerance::Float64 = 1/10^10, what = "volume")

  # Dimension
  const n = size(S1, 1)

  # Centroid and radii
  const c1 = Circumsphere(S1)[2:n+1]
  const c2 = Circumsphere(S2)[2:n+1]
  const r1 = Circumsphere(S1)[1]
  const r2 = Circumsphere(S2)[1]

  # Orientation of simplices
  const orientation_S1 = det([ones(1, n + 1); S1])
  const orientation_S2 = det([ones(1, n + 1); S2])

  if abs(orientation_S1) < tolerance || abs(orientation_S2) < tolerance
    return 0
  end

  # Set volume to zero initially. Change only if there is intersection
  IntVol = 0
  IntVert = []

# -------------------------------------
# Simplices intersect in some way
# -------------------------------------

  # If the (distance between centroids)^2-(sum of radii)^2 < 0,
  # then the simplices intersect in some way.
  dist_difference::Float64 = ((c1 - c2).' * (c1 - c2) - (r1 + r2)^2)[1]

  if dist_difference < 0
    # Find the number of points of each simplex contained within the
    # circumsphere of the other simplex

    vertices1InCircum2 = SomeVertexInCircumsphere(S1, r2, c2)
    vertices2InCircum1 = SomeVertexInCircumsphere(S2, r1, c1)

    # At least one circumsphere contains vertices of the other simplex

    if vertices1InCircum2 + vertices2InCircum1 >= 1
      ConvexExp1in2, ConvexExp2in1, ordered_vertices1, ordered_vertices2, numof1in2, numof2in1 = BarycentricCoordinates(S1,S2,orientation_S1,orientation_S2,tolerance)

      # Trivial intersections
      TriviallyContained = heaviside0([numof1in2 numof2in1] - (n+1))
      IsSomeContained = sum(TriviallyContained, 2)[1]

      if IsSomeContained == 2.0 # The simplices coincide
        IntVol = abs(orientation_S1)
      elseif IsSomeContained == 1.0 # One simplex is contained in the other
        if TriviallyContained[1] == 1.0 # Simplex1 is contained in Simplex2
          IntVol = abs(orientation_S1)
        else # Simplex2 is contained in Simplex1
          IntVol = abs(orientation_S2)
        end
      else # No simplex contains the other

        Ncomm, ordered_vertices1, ordered_vertices2 = SharedVertices(ConvexExp1in2,ordered_vertices1,ordered_vertices2,numof1in2,numof2in1)

        # Is there any shared face?
        if Ncomm == n
          #print("The simplices share a face")
          IntVol = SharingAFace(S2, ConvexExp1in2, ordered_vertices1, ordered_vertices2)

        else # The simplices do not share a face.
          IntVert, ConvexExpIntVert  = IntersectionOfBoundaries(S1,S2,ConvexExp1in2,ConvexExp2in1,ordered_vertices1,ordered_vertices2,numof1in2,numof2in1,Ncomm,tolerance);
          dim = size(IntVert, 2)

          if dim > 1
            IntVert,ConvexExpIntVert = PolytopeGeneratingVertices(S1,S2,IntVert,ConvexExpIntVert,ConvexExp1in2,ConvexExp2in1,ordered_vertices1,ordered_vertices2,numof1in2,numof2in1,Ncomm);
            IntVol = VolumeComputation(IntVert, ConvexExpIntVert)
          end
        end
      end
    else
      #println("No circumsphere of either simplex contains vertices of the other simplex")
    end
  end

  #if what == "volume"
  #  return IntVol[1]
  #elseif what == "vertices"
  #  return IntVert
  #elseif what == "both"
  #  return IntVol[1], IntVert
  #end
  return IntVol
end

export simplexintersection,
        childsimplex,
        outsidepoints, insidepoints,
        issingular,
        radius, centroid, orientation, volume,
            Circumsphere, nontrivially_intersecting_simplices, simplices_sharing_vertices, intersecting_simplices, nd_Test
end #module

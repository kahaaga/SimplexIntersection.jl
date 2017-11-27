function VolumeBT(polytope, convex_expansion_polytopevertices, n)
  #println("\nVolumeBT()\n")
  IntVol = 0
  dim = size(polytope, 1)
  #println("\n")
  # The simplices intersect in a simplex
  if dim == n + 1
    #println("The simplices intersect in a simplex\n")
    IntVol = abs(det([ones(1, n + 1); transpose(polytope)]))

    #The simplices intersect in a polytope
  elseif dim > n + 1
    #println("The simplices intersect in a polytope\n")

    # triangulation_polfaces of all polytope faces
    triangulation_polfaces, n_nonsingular_faces = TriangulationPolytopeFaces(convex_expansion_polytopevertices, dim, n)
    #@show triangulation_polfaces
    # If the number of true faces (nonsingular = object of dimension n - 1) is less than
    # n + 1, then the faces don't enclose a volume.
    if (n_nonsingular_faces >= n + 1)
      #println("numofNonSimplicialPolytopefaces >= n + 1")

      # Centroid of the polotype
      polytope_centroid = ones(1, dim) * polytope/dim
      D = size(triangulation_polfaces, 1)
      IntVol = zeros(D, 1)
      for a = 1:D
        #Ind contains the indices of the intersecting points (polytope)
        #that furnish each simplex in the triangulation_polfaces of the faces
        #of the polytope
        Ind = round.(Int64, triangulation_polfaces[a, :])
        #@show size(polytope)
        #@show Ind, size(Ind)
        #@show polytope_centroid, size(polytope_centroid)
        #@show polytope
        Vertices = [polytope[vec(Ind), :]; polytope_centroid]
        IntVol[a] = abs(det([ones(n + 1, 1) Vertices]))
      end
      IntVol = ones(1, D) * IntVol
    end
  else
    IntVol = [0]
  end

  return IntVol
end

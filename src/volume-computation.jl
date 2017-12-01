function VolumeComputation(polytope, convex_expansion_of_polytope)

    IntVol = 0
    num_polytope_generators = size(polytope, 1)
    dim = size(polytope, 2)

    if num_polytope_generators == dim + 1
        IntVol = abs(det([ones(1, dim + 1); polytope.']))
    elseif num_polytope_generators > dim + 1
        # The simplices intersect in a polytope
        #return convex_expansion_of_polytope
        triangulation_polfaces, n_nonsingular_faces = TriangulationPolytopeFaces(convex_expansion_of_polytope, num_polytope_generators, dim)

        if n_nonsingular_faces >= dim + 1
            polytope_centroid = ones(1, num_polytope_generators) * polytope / num_polytope_generators
            D = size(triangulation_polfaces, 1)
            IntVol = zeros(D, 1)
            for a = 1:D
                Ind = round.(Int64, triangulation_polfaces[a, :])
                Vertices = [polytope[vec(Ind), :]; polytope_centroid]
                IntVol[a] = abs(det([ones(dim + 1, 1) Vertices]))
            end
            IntVol = ones(1, D) * IntVol
        end
    end
    return IntVol[1]
end

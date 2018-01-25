using SimplexIntersection
using SimplexSplitting

# include("../src/Circumsphere.jl")
# include("../src/barycentric-coordinates.jl")
# include("../src/ShareFace_nD.jl")
# include("../src/geometry.jl")
# include("../src/heaviside.jl")
# include("../src/complementary.jl")
# include("../src/InfoVerticesOutside.jl")
# include("../src/intersection-of-boundaries.jl")
# include("../src/some-vertex-in-circumsphere.jl")
# include("../src/ConvexExpansion.jl")
# include("../src/Binary.jl")
# include("../src/BoundaryxBoundary.jl")
# include("sharing-a-face.jl")
# include("QR.jl")
# include("UpdateNonZeroSearchingIndex.jl")
# include("MinimalBoundaries.jl")
# include("ConvexExpAndIntVert.jl")
# include("Update.jl")
# include("IntersectingOfPolytopeVertices.jl")
# include("polytope-generating-vertices.jl")
# include("VolumeBT.jl")
# include("TriangulationPolytopeFaces.jl")
# include("NullSpace.jl")
# include("TriangulationNonSimplicialFaces.jl")
# include("SimplexChecks.jl")
# include("volume-computation.jl")
# include("shared-vertices.jl")
# include("simplexoperations.jl")

# Splits a canonical simplex of dimension E with size division factor k. Inside the simplex,
# it generates a random subsimplex. Compute its volume by computing the intersection of
# that simplex with the simplices forming the splitting of the original simplex. Compare
# the volume obtained by the simplex intersection routine to the analytical volume.
# N is the number of times to perform the test, and tolerance is the usual tolerance in the
# simplex intersection routine.
function nD_timing(ks, dims, N; tolerance = 1/10^12, plot = false)


    # Define vertices of canonical simplex
    Times = zeros(size(dims,1),size(ks,1))
    num_int_vert = zeros(size(dims,1),size(ks,1))
    max_discrepancy = 0

    for ind = 1:size(dims,1)
        E = dims[ind]
        for jnd = 1:size(ks,1)
            #counter = 0
            T = 0
            k = ks[jnd]
            canonical_simplex_vertices = zeros(E + 1, E)
            canonical_simplex_vertices[2:(E+1), :] = eye(E)
            simplex_indices = zeros(Int, 1, E + 1)
            simplex_indices[1, :] = round.(Int, collect(1:E+1))

            refined = refine_triangulation(canonical_simplex_vertices, simplex_indices, [1], k)
            triang_vertices, triang_simplex_indices = refined[1], refined[2]

            #differences = Vector{Float64}(N)

            # Repeat the test N times
            for i = 1:N
                # Build convex expansion coefficients for the random simplex. We need these in a
                # form that guarantees that the random simplex lies within the canonical simplex.
                beta = rand(E + 1, E + 1)
                beta = beta .* repmat(1 ./ sum(beta, 2), 1, E + 1)
                # Ensure that we have convex expansions

                # Linear combination of the original vertices and the convex expansion coefficients.
                # Creates a matrix containing the vertices of the random simplex (which now is
                # guaranteed to lie _within_ the original simplex.
                random_simplex = beta * canonical_simplex_vertices
                random_simplex_orientation = det(hcat(ones(E + 1), random_simplex))
                random_simplex_volume = abs(random_simplex_orientation)
                random_simplex = random_simplex.'

                # Intersection between each of the subsimplices with the random simplex
                intersecting_volumes = Vector{Float64}(size(triang_simplex_indices, 1))

                for j = 1:size(triang_simplex_indices, 1)
                    # Get the subsimplex vertices
                    subsimplex = triang_vertices[triang_simplex_indices[j, :], :].'
                    tic()
                    intvol, num_p = intersection(random_simplex, subsimplex, tolerance = tolerance)
                    timef = toc();
                    if (timef > T)
                        T = timef
                    end

                    intersecting_volumes[j] = intvol
                    if (num_p > num_int_vert[ind,jnd])
                        num_int_vert[ind,jnd] = num_p
                    end
                end

                # Compute the discrepancies
                numeric_volume = sum(intersecting_volumes)
                analytic_volume = random_simplex_volume
                tmp = abs((numeric_volume - analytic_volume)/analytic_volume)
                if (tmp>max_discrepancy)
                    max_discrepancy = tmp
                end

            end

            Times[ind,jnd] = T
        end
    end

    return Times,num_int_vert,max_discrepancy
end




function intersection(S1::Array{Float64, 2}, S2::Array{Float64, 2}; tolerance::Float64 = 1/10^10)

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
    num_p = 0

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
                num_p = n +1
                IntVol = abs(orientation_S1)
            elseif IsSomeContained == 1.0 # One simplex is contained in the other
                num_p = n+1
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
                    num_p = n+1
                    IntVol = SharingAFace(S2, ConvexExp1in2, ordered_vertices1, ordered_vertices2)

                else # The simplices do not share a face.
                    IntVert, ConvexExpIntVert  = IntersectionOfBoundaries(S1,S2,ConvexExp1in2,ConvexExp2in1,ordered_vertices1,ordered_vertices2,numof1in2,numof2in1,Ncomm,tolerance);
                    dim = size(IntVert, 2)
                    num_p = size(IntVert,1)

                    if dim > 1

                        IntVert,ConvexExpIntVert = PolytopeGeneratingVertices(S1,S2,IntVert,ConvexExpIntVert,ConvexExp1in2,ConvexExp2in1,ordered_vertices1,ordered_vertices2,numof1in2,numof2in1,Ncomm);
                        num_p = size(IntVert,1)
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
    return IntVol, num_p
end

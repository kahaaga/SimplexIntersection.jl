function TriangulationPolytopeFaces(ConvexExpIntVert, dim, n)
    TriangNSFaces = 0
    TriangAllPolytopeFaces = 0

    FacesContainingIntVert = heaviside0(-ConvexExpIntVert)
    # transforms any zero convex parameter into 1 and any positive convex
    # parameter into 0.

    NumofIntPointsInEachFace = round.(Int64, ones(1, dim) * FacesContainingIntVert)
    # array of dimension 1 x 2*n+2
    PolytopeFacesIndices = transpose(heaviside0(NumofIntPointsInEachFace - n)) .* (1:2*n+2)
    PolytopeFacesIndices = vec(find(PolytopeFacesIndices))

    # contains the elements in 1:2*n+2 corresponding to simplex faces containing n or more points from IntVert.
    # These are the indices of the potential faces of the polytope
    # At this stage different indices in PolytopeFaces might correspond to
    # the same actual face

    # The number of faces that the intersection region polytope has
    numofPolFaces = length(PolytopeFacesIndices)

    # The number vertices furnishing each polytope face. Array of dimension 1 x numofPolFaces
    NumofIntPointsInEachFace = round.(Int64, NumofIntPointsInEachFace[PolytopeFacesIndices])

    # The indices of the intersecting points (in IntVert) furnishing each polytope face.
    # Array of dimension dim x numofPolFaces. Each column represents a potential polytope face.
    IntPointsInEachPolFace = FacesContainingIntVert[:, PolytopeFacesIndices] .* repmat(collect(1:dim), 1, numofPolFaces)


    # Go through each faces and decide whether it is a true face or a boundary of a face.
    # Disregard stuff if it is not a true face. 'NonSingularIndices' starts out with only
    # zeros. If the potential polytope face is a true face, then set the value for that
    # face to 1.
    NonSingularIndices = zeros(size(PolytopeFacesIndices))

    for a = 1:numofPolFaces
        Indices = find(IntPointsInEachPolFace[:, a])
        Aux = heaviside0(-ones(1, size(Indices, 1)) * ConvexExpIntVert[Indices, :])
        Multiplicity = [Aux * ones(2*n + 2, 1);
                        Aux * [ones(n + 1, 1);
                        2 * ones(n + 1, 1)]]
        #Multiplicity(1): number of times that the face with index PolytopeFacesIndices(a) appears
        #Multiplicity(2): (number of faces of simplex1 containing the face) + 2*(number of faces of simplex2 containing the face)
        if (Multiplicity[1] == 1 || (Multiplicity[1] == 2 &&
                                    Multiplicity[2] == 3 &&
                                    PolytopeFacesIndices[a] <= n + 1))
            NonSingularIndices[a] = 1
        end
    end

    # Get back the corresponding indices of the true polytope faces (only nonzeros make
    # sense, hence find())

    NonSingularIndices = find(NonSingularIndices .* collect(1:numofPolFaces))
    numofNSPolFaces = size(NonSingularIndices, 1)

    if numofNSPolFaces >= n + 1
        #println("# Non-simplicial faces >= n + 1\n")
        # all these arrays contain the corresponding information but only for the
        # non singular faces. And therefore, the second dimension goes over 1:numofNSPolFaces
        NonSingularPolytopeFaces = PolytopeFacesIndices[NonSingularIndices]
        NumOfVertInEachFace = NumofIntPointsInEachFace[NonSingularIndices]
        NonSingularIntPointsInEachPolFace = IntPointsInEachPolFace[1:end, vec(NonSingularIndices)]
        # THe number of polytope faces that are actually simplices
        Simplicial = vec(find(heaviside0(n - NumOfVertInEachFace) .* collect(1:numofNSPolFaces)))
        NonSimplicialFaces = 0

        if length(Simplicial) == 0
            #println("No polytope faces are simplices")
            NonSimplicialFaces = NonSingularPolytopeFaces
            NonSimplicial = 1:numofNSPolFaces
            VerticesSFaces = 0
        elseif length(Simplicial) == numofNSPolFaces
            #println("All polytope faces are simplices")

            VerticesSFaces = NonSingularIntPointsInEachPolFace
            inner = reshape(VerticesSFaces, size(Simplicial, 1) * dim, 1)
            inner_nonzeros = inner[find(inner)]
            inner_reshaped_transposed = transpose(reshape(inner_nonzeros, n, size(Simplicial, 1)))
            VerticesSFaces = inner_reshaped_transposed
        else
            #println("Some polytope faces are simplices")

            VerticesSFaces = NonSingularIntPointsInEachPolFace[:, Simplicial]
            inner = reshape(VerticesSFaces, size(Simplicial, 1)*dim, 1)

            inner_nonzeros = inner[find(inner)]
            inner_reshaped_transposed = transpose(reshape(inner_nonzeros, n, size(Simplicial, 1)))
            VerticesSFaces = inner_reshaped_transposed
            NonSimplicial = complementary(Simplicial, numofNSPolFaces)
            NonSimplicialFaces = NonSingularPolytopeFaces[NonSimplicial]
        end


        if (NonSimplicialFaces[1] > 0)
            #println("One or more nonsimplical faces.")

            VerticesNSFaces = NonSingularIntPointsInEachPolFace[:, NonSimplicial]
            SimplexIndexNS = round.(Int64, ceil.(NonSimplicialFaces / (n + 1)))
            SimplexFaceIndexNS = round.(Int64, NonSimplicialFaces - (SimplexIndexNS - 1) * (n + 1))

            TriangNSFaces = TriangulationNonSimplicialFaces(VerticesNSFaces, SimplexIndexNS, SimplexFaceIndexNS, ConvexExpIntVert, n)
        else
            #println("No nonsimplical faces.")
        end
        #VerticesSFaces
        TriangAllPolytopeFaces = Update(VerticesSFaces, TriangNSFaces)
    else
        #println("# Non-simplicial faces < n + 1\n")

    end
    return TriangAllPolytopeFaces, numofNSPolFaces
end

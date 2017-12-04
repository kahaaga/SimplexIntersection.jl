using .QHull.delaunay_tesselation


#TODO:: Needs revision. Doesn't work.
function TriangulationNonSimplicialFaces(VerticesNSFaces,SimplexIndexNS,SimplexFaceIndexNS,ConvexExpIntVert,n)
    #println("\n\n@TriangulationNonSimplicialFaces...\n\n")
    TriangNSFaces = [0]
    numofNSFaces = size(VerticesNSFaces, 2)

    for a = 1:numofNSFaces
        # Contains the indices of the elements in IntVert generating the corresponding
        # polytope face.
        Indices = find(VerticesNSFaces[:, a])

        FaceVert = complementary(SimplexFaceIndexNS[a], n + 1) + (SimplexIndexNS[a] - 1) * (n + 1)
        BarCoordinates = ConvexExpIntVert[Indices, FaceVert]

        # The rows of this array contain the (strictly positive) concvex exp. coefficients of
        # the IntPoints generating the polytope face, in terms of the vertices
        # generating the simplex sharing the boundary with the polytope face
        i = size(BarCoordinates, 1)
        BarCoordinates = [ones(1, i); transpose(BarCoordinates)]

        # The columns of Null form a basis of the Null space of BarCoordinates
        # PermC is a permutation of the columns of BarCoordinates (or the
        # rows of Null) such that Null(PermC,:)=[M;eye(N)] with M some
        # matrix and N the number of columns in Null.
        Null, PermC = NullSpace(BarCoordinates, n)

        CopBarCoord = zeros(n, i)
        CopBarCoord[:, PermC[1:n]] = eye(n)
        CopBarCoord[:, PermC[(n + 1):i]] = - Null[PermC[1:n], :]

        CoplanarCoord=CopBarCoord[2:n, :]
        CoplanarCoord[:, PermC[1]] = zeros(n - 1, 1)
        C = transpose(CoplanarCoord)

        # Triangulating
        points, simplices = delaunay_tesselation(C)
        #disp('you filthy bastard')

        NewTriang = simplices

        # Notice that the order in the argument of delaunayn
        # corresponds to the order given by Indices
        # The rows of NewTriang contain the indices of the vertices generating
        # the corresponding simplex in the order given by Indices, on the other hand
        # the entries of Indices correspond to the indices of the original intersecting
        # points in IntVert
        d = size(NewTriang, 1)

        # Triangulation of the non simplicial face a. The indices of the vertices generating
        # the simplices now correspond to the order given by IntVert
        Aux = reshape(transpose(NewTriang), n * d, 1)
        NewTriang = transpose(reshape(Indices[Aux], n, d))

        if a == 1
            TriangNSFaces = NewTriang;
        else
            TriangNSFaces = [TriangNSFaces; NewTriang]
        end
    end
    return TriangNSFaces
end

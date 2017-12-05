
"""
ShareFace_nD()

If simplices share an entire face, then the intersecting volume is easily
computable, and you don't need to compute all intersecting point. This
function does exactly that. When two simplices share a face, they will either
intersect with zero volume (they only share a face) OR the intersection will
be another simplex. This is true in any dimension.

Input arguments
---------------
source_simplex::Array{Float64, 2} Source simplex. Represented a matrix of dimension nx(n+1), where each column vector is a vertex.
target_simplex::Array{Float64, 2} Target simplex. Rrepresented a matrix of dimension nx(n+1), where each column vector is a vertex.

indices_of_shared_vertices_in_source::Vec{Int}
indices_of_shared_vertices_in_::Vec{Int}
O_Reference::Float64                           Orientation of source simplex.
O_Target::Float64                              Orientation of target simplex.
tolerance::Float64


Returns
-------
intersecting_volume::Float64                   The intersection between the simplices.
"""
function ShareFace_nD(Reference_Simplex, Target_Simplex, IndexComVert_Reference,
    IndexComVert_Target, O_Reference, O_Target, tolerance)

    # IndexComvert_Reference and IndexComvert_Target are ordered such that they correspond to the same points
    n = size(Reference_Simplex, 1)

    # Array of dimensions 1x1. Index of the vertex NOT shared by the two simplices for the reference simplex (all other vertices are shared because they form the shared face)
    IndexExtra_Reference = complementary(IndexComVert_Reference, n + 1)

    # Array of dimensions 1x1.  Index of the vertex NOT shared by the two simplices for the reference simplex (all other vertices are shared because they form the shared face)
    IndexExtra_Target = complementary(IndexComVert_Target, n + 1)

    # Convex expansion of the extra index of Reference_Simplex in terms of the
    # vertices of Target_Simplex. Array of dimensions 1x(n+1)
    ConvexExpansion_Reference = zeros(1, n + 1)

    # Take the extra vertex of the reference simplex and express it as a linear convex
    # combination of the vertices of the target simplex.
    for j=1:n+1
        Aux = copy(Target_Simplex)
        ref = copy(Reference_Simplex)
        Aux[:, j] = ref[:, IndexExtra_Reference]
        beta = det([ones(1, n + 1); Aux]) / O_Target
        if abs(beta) < tolerance
            beta = 0
        end
        ConvexExpansion_Reference[j] = beta
    end

    minimum_coefficient = minimum(ConvexExpansion_Reference)

    # Maximum coefficient of the coefficients in terms of the shared vertices
    maximum_coeff_n = maximum(ConvexExpansion_Reference[IndexComVert_Target])

    # Convex coefficient appearing in front of the extra (non-shared) vertex of the target simplex
    extra_convex_coeff_= ConvexExpansion_Reference[IndexExtra_Target]


    if ConvexExpansion_Reference[IndexExtra_Target][1] >= 0
        # S1 is contained in S2

        if minimum_coefficient >= 0
            IntVol = abs(O_Reference)
            IntVert = Reference_Simplex

        # S2 is contained in S1

        elseif maximum_coeff_n <= 0
            IntVol = abs(O_Target)
            IntVert = Target_Simplex
        # S1 and S2 intersect non-trivially, such that the intersection is a new simplex.
        else
            negativeCoeff = heaviside(-ConvexExpansion_Reference[IndexComVert_Target]).' # row vector
            NegativeIndices = round.(Int64, negativeCoeff .* IndexComVert_Target.') # row vector

            tmp = IndexComVert_Target.' - NegativeIndices
            NonNegativeIndices = tmp[find(tmp)].' # row vector

            # Should be a number
            Sigma = sum(ConvexExpansion_Reference[NonNegativeIndices])

            # Column vector
            TargetCoefficients = ConvexExpansion_Reference[[NonNegativeIndices IndexExtra_Target]] / ((ConvexExpansion_Reference[IndexExtra_Target] + Sigma))[1]

            # Should be a vertex (collumn vectof) MY PRECIOUS
            IntPoint = Target_Simplex[:, vec([NonNegativeIndices IndexExtra_Target].')] * TargetCoefficients.'

            IntVert = [Target_Simplex[:, IndexComVert_Target] IntPoint].'

            IntVol = abs(det([ones(n + 1, 1) IntVert]))
        end

    # The simplices only have the shared face in common, so they have no shared volume.

    else
        IntVol = 0
        IntVert = []
    end

    return IntVol, IntVert
end

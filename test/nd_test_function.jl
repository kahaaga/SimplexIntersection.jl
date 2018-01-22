using SimplexIntersection
using SimplexSplitting

# Splits a canonical simplex of dimension E with size division factor k. Inside the simplex,
# it generates a random subsimplex. Compute its volume by computing the intersection of
# that simplex with the simplices forming the splitting of the original simplex. Compare
# the volume obtained by the simplex intersection routine to the analytical volume.
# N is the number of times to perform the test, and tolerance is the usual tolerance in the
# simplex intersection routine.
function nd_Test(k, E, N; tolerance = 1/10^12, plot = false)


    # Define vertices of canonical simplex
    canonical_simplex_vertices = zeros(E + 1, E)
    canonical_simplex_vertices[2:(E+1), :] = eye(E)
    simplex_indices = zeros(Int, 1, E + 1)
    simplex_indices[1, :] = round.(Int, collect(1:E+1))

    refined = refine_triangulation(canonical_simplex_vertices, simplex_indices, [1], k)
    triang_vertices, triang_simplex_indices = refined[1], refined[2]

    differences = Vector{Float64}(N)

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
            intvol = simplexintersection(random_simplex, subsimplex, tolerance = tolerance, what = "volume")
            intersecting_volumes[j] = intvol
        end

        # Compute the discrepancies
        numeric_volume = sum(intersecting_volumes)
        analytic_volume = random_simplex_volume
        differences[i] = abs((numeric_volume - analytic_volume)/analytic_volume)

    end
    return differences
end

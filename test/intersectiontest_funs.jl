"""
    intersectiontest(dim::Int, n_reps::Int, intersection_type::String, tolerance::Float64)

Generate a `n_reps` pairs of `dim`-dimensional simplices that intersect in the way specified by `intersecton_type`. Valid intersection types are `"nontrivial"` and `"sharingvertices"`. Checks if the total volume spanned by each pair of simplices matches that obtained when accounting for the shared volume.
"""
function intersectiontest(dim::Int, n_reps::Int, intersection_type::String; tolerance = 1/10^12)

    discrepancies = zeros(Float64, n_reps)
    elapsed = zeros(Float64, n_reps)
    for i = 1:n_reps
        s₁, s₂ = intersecting_simplices(dim = dim, intersection_type = intersection_type)
        before = time_ns()/10^9
        intvol = simplexintersection(s₁.', s₂.', tolerance = tolerance, what = "volume")
        after =  time_ns()/10^9
        elapsed[i] = after - before
        vol₁ = volume(s₁)
        vol₂ = volume(s₂)
        total_vol = vol₁ + vol₂
        discrepancies[i] = total_vol - ((vol₁ - intvol) + (vol₂ - intvol) + intvol)
        if discrepancies[i] > 1/10^5
            @show vol₁, vol₂, discrepancies[i]
        end
    end
    return discrepancies, elapsed
end

function InsideCircum(S, r, C, n)

    C = repmat(C, 1, n + 1)
    rsquared = repeat([r^2], outer = [n + 1])
    tmp1 = r^2 * ones(1, n+1)
    tmp2 = ones(1, n) * (S - C).^2
    inds = round.(Int64, heaviside0(tmp1 - tmp2) .* transpose(collect(1:n+1)))
    inds = find(vec(inds)) # Disregard zero indices
    nverts = countnz(inds)::Int # Count nonzeros

  return ifelse(nverts == 0, [0], inds), nverts
end

function inside_circumsphere(S1, S2)
    # Centroid and radii
    n = size(S1, 1)

    r2 = Circumsphere(S2)[1]
    c2 = Circumsphere(S2)[2:n+1]

    C = repmat(c2, 1, n + 1)
    rsquared = repeat([r2^2], outer = [n + 1])
    tmp1 = r2^2 * ones(1, n+1)
    tmp2 = ones(1, n) * (S1 - C).^2
    inds = heaviside0(tmp1 - tmp2) .* transpose(collect(1:n+1))
    inds = round.(Int64, inds)
    inds = find(vec(inds)) # Disregard zero indices
    nverts = countnz(inds)::Int # Count nonzeros

  return ifelse(nverts == 0, [0], inds), nverts
end

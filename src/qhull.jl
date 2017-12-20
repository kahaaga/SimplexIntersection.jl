## This is a modification of the code reference below:
##
## CHull.jl
## (c) 2013 David Al van Leeuwen
## A Julia wrapper around a python wrapper around the qhull Convex Hull library

## This code is licensed under the GNU General Public License, version 2
## See the file LICENSE in this distribution

module QHull

export ConvexHull, convexhull, DelaunayTesselation, delaunay_tesselation#, display, show

using PyCall
const spatial = PyNULL()

function __init__()
    copy!(spatial, pyimport_conda("scipy.spatial", "scipy"))
end

type ConvexHull{T<:Real}
    points::Array{T}
    vertices::Vector{Int}
end

type Triangulation{T<:Real}
    points::Array{Float64, 2}
    simplices::Array{Int64, 2}
end

function convexhull{T<:Real}(points::Array{T})
    py = spatial[:ConvexHull](points)
    points = convert(Array{T}, py["points"])
    vertices = convert(Array{Int}, py["vertices"])
    return points, vertices
    #ConvexHull(points, vertices)
end

function delaunay_tesselation{T<:Real}(points::Array{T})
    py = spatial[:Delaunay](points, furthest_site = true, qhull_options = ["Qt Qbb Qc Qx QJ"])
    points = convert(Array{T}, py["points"])

    # Get simplices forming the triangulation
    s = convert(Vector{Vector{Int}}, py["simplices"])
    simplices = convert(Vector{Vector{Int}}, py["simplices"])

    # Convert simplices to array
    n_simplices = length(simplices)
    dim = size(points)[2]

    simplices_arr = zeros(n_simplices, dim + 1)

    for i in 1:n_simplices
      simplices_arr[i, :] = s[i]
    end
    return points, ceil.(Int64, simplices_arr) + 1 # Add 1 to account for base difference in indices
end

end

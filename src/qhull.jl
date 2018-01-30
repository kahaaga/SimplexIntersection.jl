## This is a modification of the code reference below:
##
## CHull.jl
## (c) 2013 David Al van Leeuwen
## A Julia wrapper around a python wrapper around the qhull Convex Hull library

## This code is licensed under the GNU General Public License, version 2
## See the file LICENSE in this distribution

module QHull

export ConvexHull, delaunayn, delaunay_tesselation

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
end

function delaunayn(points::Array{Float64})
    py = spatial[:Delaunay](points)
    indices = Array{Int64, 2}(length(py["simplices"]), size(points, 2) + 1)
    pyarray_to_array!(py["simplices"], indices, Int)

    return indices .+ 1 # Add 1 to account for base difference in indices
end

function pyarray_to_array!(pyobject, arr, T)
    for i = 1:length(pyobject)
        arr[i, :] = get(pyobject, PyVector{T}, i-1) # i-1 because of Python 0 indexing
    end
end

function delaunay_tesselation{T<:Real}(points::Array{T})
    py = spatial[:Delaunay](points)
    #pts = convert(Array{T}, py["points"])

    # Get simplices forming the triangulation
    s = convert(Vector{Vector{Int}}, py["simplices"])

    # Convert simplices to array
    n_simplices = length(s)
    dim = size(points, 2)
    simplices_arr = zeros(Int, n_simplices, dim + 1)

    for i in 1:n_simplices
      simplices_arr[i, :] = s[i]
    end
    #return pts, ceil.(Int64, simplices_arr) + 1 # Add 1 to account for base difference in indices
    return simplices_arr .+ 1 # Add 1 to account for base difference in indices

end

end

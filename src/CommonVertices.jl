#module commonvertices

"""
    CommonVertices()

Finds the common vertices of simplices X and Y.

First check how many vertices of simplex X is contained in the circumsphere of simplex Y,
 then vice versa. Common simplices exist only if both circumspheres contain simplices of
 each other.

Input arguments
---------------
Simplex1::Array{Float64, 2}         Vertices of simplex 1. Columns represent points.
Simplex2::Array{Float64, 2}         Vertices of simplex 2. Columns represent points.
index1in2::Vector{Int64}     Either a row vector of ints, or 0. Contains all the vertices of simplex X contained in the circumsphere of simplex Y.
index2in1::Vector{Int64}     Either a row vector of ints, or 0. Contains all the vertices of simplex Y contained in the circumsphere of simplex X.
numof1in2::Int               Number of elements in 'index1in2'
numof2in1::Int               Number of elements in 'index2in1'

Returns
-------
Ncomm                              Number of common vertices.
IndexComVert1::Vector{Int64}     Row vector containing indices simplex 1 vertices that are shared by both simplices. This is used in function 'ContainedVertices' later.
IndexComVert2::Vector{Int64}     Row vector containing indices simplex 2 vertices that are shared by both simplices. This is used in function 'ContainedVertices' later.
InternalComIndex1::Vector{Int64}   Indices of the common vertices in the vector 'index1in2'
InternalComIndex2::Vector{Int64}    Indices of the common vertices in the vector 'index2in1'
"""
function CommonVertices(X::Array{Float64, 2},
                        Y::Array{Float64, 2},
                        index1in2,
                        index2in1,
                        numof1in2::Int,
                        numof2in1::Int,
                        n::Int)

    Ncomm = 0
    InternalComIndex1 = [0]
    InternalComIndex2 = [0]
    IndexComVert1 = [0]
    IndexComVert2 = [0]

    if numof1in2 * numof2in1 > 0
        #COMMON VERTICES  (kron(A,B))iaja ibjb = (A)iaja (B)ibjb
        M =  numof1in2 * numof2in1

        comvert = kron(X[:, index1in2], ones(1, numof2in1))
        ComVert = kron(ones(1, numof1in2), Y[:, index2in1])

        tmp = heaviside0(-ones(1, n) * ((comvert - ComVert).^2))
        tmp = tmp .* transpose(collect(1:M))
        tmp = find(transpose(tmp))

        Ncomm = size(tmp, 1)

        if Ncomm > 0
            InternalComIndex1 = ceil.(Int64, tmp / numof2in1)
            InternalComIndex2 = tmp - numof2in1 * (InternalComIndex1 - 1)

            IndexComVert1 = index1in2[InternalComIndex1]
            IndexComVert2 = index2in1[InternalComIndex2]
        end
    end
    InternalComIndex1 = ifelse(length(InternalComIndex1) == 0, [0], InternalComIndex1)
    InternalComIndex2 = ifelse(length(InternalComIndex2) == 0, [0], InternalComIndex2)
    IndexComVert1 = ifelse(length(IndexComVert1) == 0, [0], IndexComVert1)
    IndexComVert2 = ifelse(length(IndexComVert2) == 0, [0], IndexComVert2)

    return Ncomm::Int64,
            InternalComIndex1::Vector{Int64},
            InternalComIndex2::Vector{Int64},
            IndexComVert1::Vector{Int64},
            IndexComVert2::Vector{Int64}
end

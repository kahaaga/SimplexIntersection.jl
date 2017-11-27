"""
IntersectionPolytopeVertices()

Updates IntVert and ConvexExpIntVert according to the following:

In the intersection of boundaries, we have excluded all the vertices of the original that actually lie inside the other simplex.
This function computes the intersecting points between boundaries that are not original vertices any of the simplices.

Input arguments
---------------
X
Y
IntVert
ConvexExpIntVert
Vert1Inside2
Vert2Inside1
numof1in2NotCom
numof2in1NotCom
Ncomm
n

Returns
-------
IntVert
ConvexExpIntVert

"""

function IntersectionPolytopeVertices(X,Y,IntVert,ConvexExpIntVert,Vert1Inside2,Vert2Inside1,numof1in2NotCom,numof2in1NotCom,Ncomm,n)

    X = copy(X)
    Y = copy(Y)

    if Ncomm > 0
        IndexComVert1 = round.(Int64, Vert1Inside2[1, 1:Ncomm])
        IntVert = [IntVert; transpose(X[:, IndexComVert1])]
        Aux1 = transpose(Vert2Inside1[2:(n + 2), 1:Ncomm])
        Aux2 = transpose(Vert1Inside2[2:(n + 2), 1:Ncomm])
        ConvexExpIntVert = [ConvexExpIntVert; [Aux1 Aux2]]
    end

    if numof1in2NotCom * numof2in1NotCom > 0
        IntInd1in2 = Ncomm+1:size(Vert1Inside2, 2)
        IntInd2in1 = Ncomm+1:size(Vert2Inside1, 2)
        IndexVert1in2 = round.(Int64, Vert1Inside2[1, IntInd1in2])
        IndexVert2in1 = round.(Int64, Vert2Inside1[1, IntInd2in1])
        IntVert = [IntVert;transpose([X[:, IndexVert1in2] Y[:, IndexVert2in1]])]
        Aux1in2_1 = zeros(numof1in2NotCom * (n + 1), 1)
        Aux1in2_1[collect(0:numof1in2NotCom - 1) * (n + 1) + IndexVert1in2] = ones(numof1in2NotCom, 1)
        Aux1in2_1 = transpose(reshape(Aux1in2_1, n + 1, numof1in2NotCom))
        Aux1in2_2 = transpose(Vert1Inside2[2:n+2, IntInd1in2])
        Aux2in1_2 = zeros(numof2in1NotCom*(n+1),1)
        Aux2in1_2[(0:numof2in1NotCom - 1) * (n + 1) + IndexVert2in1] = ones(numof2in1NotCom, 1)
        Aux2in1_2 = transpose(reshape(Aux2in1_2, n + 1, numof2in1NotCom))
        Aux2in1_1 = transpose(Vert2Inside1[2:(n + 2), IntInd2in1])
        ConvexExpIntVert = [ConvexExpIntVert; [[Aux1in2_1 Aux1in2_2]; [Aux2in1_1 Aux2in1_2]]]

    elseif numof1in2NotCom > 0
        IntInd1in2 = (Ncomm + 1):size(Vert1Inside2, 2)
        IndexVert1in2 = round.(Int64, Vert1Inside2[1, IntInd1in2])
        IntVert = [IntVert; transpose(X[:, IndexVert1in2])]
        Aux1in2_1 = zeros(numof1in2NotCom * (n + 1), 1)
        Aux1in2_1[(0:numof1in2NotCom-1) * (n + 1) + IndexVert1in2] = ones(numof1in2NotCom, 1)
        Aux1in2_1 = transpose(reshape(Aux1in2_1, n + 1, numof1in2NotCom))
        Aux1in2_2 = transpose(Vert1Inside2[2:(n + 2), IntInd1in2])
        ConvexExpIntVert = [ConvexExpIntVert; [Aux1in2_1 Aux1in2_2]]
    elseif numof2in1NotCom > 0
        IntInd2in1 = (Ncomm + 1):size(Vert2Inside1, 2)
        IndexVert2in1 = round.(Int64, Vert2Inside1[1, IntInd2in1])
        IntVert = [IntVert; transpose(Y[:, IndexVert2in1])]
        Aux2in1_2 = zeros(numof2in1NotCom * (n + 1), 1)
        Aux2in1_2[(0:numof2in1NotCom - 1) * (n + 1) + IndexVert2in1] = ones(numof2in1NotCom, 1)
        Aux2in1_2 = transpose(reshape(Aux2in1_2, n + 1, numof2in1NotCom))
        Aux2in1_1 = transpose(Vert2Inside1[2:(n + 2), IntInd2in1])

        ConvexExpIntVert = [ConvexExpIntVert; [Aux2in1_1 Aux2in1_2]]
    end

    return IntVert, ConvexExpIntVert
end

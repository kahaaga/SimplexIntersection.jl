include("ConvexExpAndIntVert.jl")

function BoundaryxBoundary(X,Y,BoundaryBinaryLabels1,BoundaryBinaryLabels2,IndexVert1,IndexVert2,Beta1in2,Beta2in1,tolerance)
    n = size(X, 1)
    N = size(BoundaryBinaryLabels1, 1)

    ConvexExpIntVert = 0
    BoundaryVertexIndices1 = BoundaryBinaryLabels1 .* repmat(transpose(collect(1:n+1)), N)
    BoundaryVertexIndices2 = BoundaryBinaryLabels2 .* repmat(transpose(collect(1:n+1)), N)
    SearchingIndex = collect(1:N)
    NonZeroSearchingIndex = copy(SearchingIndex)
    IntVert = 0

    for a = 1:N
        if SearchingIndex[a] > 0
            IntIndexVertB1 = transpose(find(BoundaryVertexIndices1[a, :]))
            IntIndexVertB2 = transpose(find(BoundaryVertexIndices2[a, :]))
            #Indices of the vertices of the boundaries in the order given by IndexVert1
            #and IndexVert2, respectively
            IndexVertB1 = IndexVert1[IntIndexVertB1]
            IndexVertB2 = IndexVert2[IntIndexVertB2]


            #Indices of the vertices of the boundaries in the order given by X and Y, respectively
            r = size(IndexVertB1, 1)
            s = size(IndexVertB2, 1)
            TargetVertices = copy(Y)
            ReferenceBoundary = copy(IndexVertB1)
            TargetBoundary = copy(IndexVertB2)
            beta = Beta2in1
            Switch = 0

            if r < s
                TargetVertices = copy(X)
                ReferenceBoundary = copy(IndexVertB2)
                TargetBoundary = copy(IndexVertB1)
                beta = Beta1in2
                aux = copy(r)
                r = s
                s = aux
                Switch = 1
            end

            ExtraIndices = complementary(ReferenceBoundary, n + 1)

            # A matrix
            Gamma = beta[vec(ExtraIndices), vec(TargetBoundary)]
            Rank = rank(Gamma)
            Rank0 = rank([Gamma; ones(1, s)])

            aux = minimum(ones(1, size(Gamma, 1)) * (Gamma .^2))

            if (Rank0 - Rank == 1 && Rank == s - 1 && aux > 0)
                lambda = QR(Gamma, tolerance)
                alpha = zeros(r, s)
                alpha[2:r, :] = beta[ReferenceBoundary[2:r], TargetBoundary]
                alpha[1, :] = 1 - ones(1, r - 1) * alpha[2:r, :]
                alpha = alpha * lambda

                # A column vector
                alpha = alpha .* heaviside(abs.(alpha) - tolerance)
                aux = minimum([minimum(alpha) minimum(lambda)])
                Aux = size(find(alpha), 1)

                #if some of this coefficients are negative the boundaries simply do
                #not intersect.
                #if some are zero, then the boundaries are not minimal and the intersecting point has already been
                #computed or will be computed as the intersection of minimal
                #boundaries.
                #This also rules out duplication of points
                if aux >= 0
                    NonZeroSearchingIndex = UpdateNonZeroSearchingIndex(NonZeroSearchingIndex, a)
                    NonZeroSearchingIndex, SearchingIndex = MinimalBoundaries(Switch, aux, alpha, lambda, IntIndexVertB1, IntIndexVertB2, BoundaryBinaryLabels1, BoundaryBinaryLabels2, NonZeroSearchingIndex, SearchingIndex)
                    if (Aux > 1)

                        # Calculate the new point. Also, check if the new point is so close to
                        # some other point that they are practically indistinguishable. If
                        # so, disregard it.

                        NewPoint = TargetVertices[:, vec(TargetBoundary)] * lambda

                        point_already_exists = false
                        if !(IntVert == 0)
                            for i = 1:size(IntVert, 1)
                                if NewPoint ≈ IntVert[i, :]
                                    point_already_exists = true
                                end
                            end
                        end

                        if !point_already_exists
                            IntVert,ConvexExpIntVert = ConvexExpAndIntVert(IntVert, transpose(NewPoint), transpose(alpha), transpose(lambda), ReferenceBoundary, TargetBoundary, Switch, ConvexExpIntVert, n)
                        end

                    end
                end
            end
        end
    end
    return IntVert, ConvexExpIntVert
end

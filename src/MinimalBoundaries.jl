function MinimalBoundaries(Switch,aux,alpha,lambda,IntIndexVertB1,IntIndexVertB2,BoundaryBinaryLabels1,BoundaryBinaryLabels2,NonZeroSearchingIndex,SearchingIndex)
    #OUTCOME
    # Colum vector of size (NumberOfPairsOfBoundaries) x 1  containing 1 if the
    # corresponding pair of boundaris contain the respective minimal
    # boundaries, 0 otherwise
    M = size(NonZeroSearchingIndex, 1)
    ReferenceBoundary = IntIndexVertB1
    TargetBoundary = IntIndexVertB2
    RefBoundaryBinaryLabels = BoundaryBinaryLabels1
    TargBoundaryBinaryLabels = BoundaryBinaryLabels2

    if Switch == 1
        ReferenceBoundary = IntIndexVertB2
        TargetBoundary = IntIndexVertB1
        RefBoundaryBinaryLabels = BoundaryBinaryLabels2
        TargBoundaryBinaryLabels = BoundaryBinaryLabels1
    end

    IndexMinimalRefBoundary = ReferenceBoundary
    #Indices of the vertices of the minimal reference boundary in
    #the order given by the corresponding IndexVert
    IndexMinimalTargBoundary = TargetBoundary
    if aux == 0
        #There are zero coefficients
        Minimal_alpha = heaviside(transpose(alpha))
        Minimal_lambda = heaviside(transpose(lambda))

        IndexMinimalRefBoundary = ReferenceBoundary .* Minimal_alpha
        #Indices of the vertices of the minimal reference boundary in
        #the order given by the corresponding IndexVert
        IndexMinimalTargBoundary = TargetBoundary .* Minimal_lambda
        #Indices of the vertices of the minimal reference boundary in
        #the order given by the corresponding IndexVert
        IndexMinimalRefBoundary = transpose(find(IndexMinimalRefBoundary))
        IndexMinimalTargBoundary = transpose(find(IndexMinimalTargBoundary))
    end
    r = size(IndexMinimalRefBoundary, 2)
    s = size(IndexMinimalTargBoundary, 2)

    Aux = RefBoundaryBinaryLabels[NonZeroSearchingIndex, transpose(IndexMinimalRefBoundary)] * ones(r, 1)
    #Column vector with as many entries as NonZeroSearchingIndex.
    #Contains the number of 1s of the binary expansion of all the reference
    #boundaries to be searched in the components of the minimal reference boundary.
    #Aux<=r. Aux=r when the corresponding boundary contains the
    #minimal one
    Aux = heaviside(Aux - r)
    #1 if the corresponding reference boundary contains the minimal boundary, 0
    #otherwise
    NonMinimalBoundaries1 = Aux;

    Aux = TargBoundaryBinaryLabels[NonZeroSearchingIndex, transpose(IndexMinimalTargBoundary)] * ones(s, 1)
    #Exactly analogous to the reference case above
    Aux = heaviside(Aux - s)
    #1 if the corresponding target boundary contains the minimal boundary, 0
    #otherwise
    NonMinimalBoundaries2 = Aux
    NonMinimalBoundaries = NonMinimalBoundaries1 .* NonMinimalBoundaries2
    #contains 1 if the pair of boundaries corresponding to the current NonZeroSearchingIndex contains the respective minimal
    #boundaries, 0 otherwise. In case of 1, the corresponding
    #NonZeroSearchingIndex entry must be set to 0
    Aux = find(NonMinimalBoundaries .* collect(1:M))
    #contains the entries of NonZeroSearchingIndex that need to be set to zero
    SearchingIndex[NonZeroSearchingIndex[Aux]] = zeros(size(Aux))
    NonZeroSearchingIndex[Aux] = zeros(size(Aux))
    NonZeroSearchingIndex = find(NonZeroSearchingIndex)

    return  NonZeroSearchingIndex, SearchingIndex
end

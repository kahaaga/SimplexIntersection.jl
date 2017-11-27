function UpdateNonZeroSearchingIndex(NonZeroSearchingIndex, a)
    N = size(NonZeroSearchingIndex, 1)


    #println()
    #println("N ", N)
    #println("a ", a)

    #println("NonZeroSearchingIndex ", NonZeroSearchingIndex)
    j = find(NonZeroSearchingIndex .== a)
    #println("j = ", j)
    NonZeroSearchingIndex = NonZeroSearchingIndex[NonZeroSearchingIndex .>= a]
    #println("NonZeroSearchingIndex ", NonZeroSearchingIndex)
    #println("----------------------------")
    return NonZeroSearchingIndex
end

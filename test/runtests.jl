using SimplexIntersection
using Base.Test

# write your own tests here
#@time include("test_simplexintersection.jl")
@time include("test_NDintersection.jl")

#include("intersectiontests.jl")
#nreps = 1000
# intersectiontest(2, 10, "sharingvertices")
#intersectiontest(3, 10, "sharingvertices")
#intersectiontest(4, 10, "sharingvertices")
# intersectiontest(5, 10, "sharingvertices")
# intersectiontest(2, 10, "nontrivial")
#intersectiontest(3, 10, "nontrivial")
# intersectiontest(4, 10, "nontrivial")
# intersectiontest(5, 10, "nontrivial")
#
# shareverts_discrep₁dim₂ = intersectiontest(2, nreps, "sharingvertices")
#shareverts_discrep₁dim₃ = intersectiontest(3, nreps, "sharingvertices")
#shareverts_discrep₁dim₄ = intersectiontest(4, nreps, "sharingvertices")
# shareverts_discrep₁dim₅ = intersectiontest(5, nreps, "sharingvertices")
#
# nontrivial_discrep₂dim₂ = intersectiontest(2, nreps, "nontrivial")
#nontrivial_discrep₂dim₃ = intersectiontest(3, nreps, "nontrivial")
# nontrivial_discrep₂dim₄ = intersectiontest(4, nreps, "nontrivial")
# nontrivial_discrep₂dim₅ = intersectiontest(5, nreps, "nontrivial")

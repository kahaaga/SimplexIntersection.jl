# SimplexIntersection

Julia package for computing the intersection between simplices in N
dimensions. Computes both the intersecting volume and the points
in which the simplices intersect.

## Requirements
SimplexIntersection.jl requires Julia 0.6 to run. The package has currently only been tested on macOS.

## Installation

To install the package, run the following command in the Julia console:

Pkg.clone("git@github.com:kahaaga/SimplexIntersection.jl.git")

This will create a local copy of the repository in your ~/.julia/v0.6/ directory. Make sure to install the dependencies that get listed during installation.

## R bindings
An R package, rSimplexIntersection, that wraps the functions in SimplexIntersection.jl is available [here](https://github.com/kahaaga/rSimplexIntersection).

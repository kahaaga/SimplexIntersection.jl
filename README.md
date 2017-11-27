# SimplexIntersection

Julia package for computing the intersection between simplices in N
dimensions. Computes both the intersecting volume and the points
in which the simplices intersect.

## Requirements
SimplexIntersection.jl requires Julia 0.6 to run. The package has currently only been tested on macOS. You can download the latest version of Julia from the official Julia [website](https://julialang.org/downloads/). 

## Installation

To install the package, run the following command in the Julia console:

```
Pkg.clone("git@github.com:kahaaga/SimplexIntersection.jl.git")
```

This will create a local copy of the repository in your `~/.julia/v0.6/` directory. Make sure to install the dependencies that get listed during installation.

## R bindings
An R package, rSimplexIntersection, that wraps the functions in SimplexIntersection.jl is available [here](https://github.com/kahaaga/r-simplexintersection).


# Examples
```julia
# Load package
using SimplexIntersection

# Define two simplices that have som overlapping volume
simplex1 = [1.12545 0.0978862 0.401808 1.91521;
        0.979468 -0.171237 0.0978862 1.12545;
        0.491734 -0.452719 -0.171237 0.979468];
simplex2 = [0.979468 -0.171237 0.0978862 1.12545;
        0.491734 -0.452719 -0.171237 0.979468;
        0.401808 -0.537488 -0.452719 0.491734];

# Compute the volume. By default, a tuple of volume and vertices is returned.
# Here, we keep only the volume.
vol = simplexintersection(simplex1, simplex2)[1]
```

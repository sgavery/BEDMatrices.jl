__precompile__()

"""
    module BEDMatrices

A collection of tools for efficiently reading or indexing into .bed
files for calculation or other manipulations. Functionality primarily
revolves around the `immutable` `BEDMatrix`; see `BEDMatrix`
documentation for more details.

"""
module BEDMatrices
include("bedmatrices.jl")

export BEDintomatrix, BEDintomatrix!,
    BEDMatrix, path, rownames, colnames, NArep

end

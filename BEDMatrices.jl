__precompile__()

"""
    module BEDMatrices

A collection of tools for efficiently reading or indexing into .bed
files for calculation or other manipulations. Functionality primarily
revolves around the `immutable` `BEDMatrix`; see `BEDMatrix`
documentation for more details.

"""
module BEDMatrices
include("constants.jl")
include("bedmatrices.jl")
include("coltools.jl")

export BEDintomatrix, BEDintomatrix!, BEDMatrix, path, rownames,
    colnames, NArep, hasNAs, countNAs, column_sum, column_mapreduce,
    column_moments, column_dot

end

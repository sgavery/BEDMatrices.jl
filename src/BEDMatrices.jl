__precompile__()

"""
    module BEDMatrices

A collection of tools for efficiently reading or indexing into .bed
files for calculation or other manipulations. Functionality primarily
revolves around the `immutable` `BEDMatrix`; see `BEDMatrix`
documentation for more details.

"""
module BEDMatrices
using Base.LinAlg.axpy!

module Consts
include("constants.jl")
end

const NA_byte = Consts.NA_byte

include("bedmatrix.jl")
include("coltools.jl")

export BEDintomatrix,
    BEDintomatrix!,
    BEDMatrix,
    path,
    rownames,
    colnames,
    NArep,
    hasNAs,
    countNAs,
    BEDdot,
    column_dist,
    column_sum,
    column_dot,
    column_NAsup_dot,
    column_sumabs2,
    column_norm,
    column_mean_and_std,
    column_mean_and_var

end

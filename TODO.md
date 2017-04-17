#### To Implement
* change `NA_byte` to something more reasonable; for Float types change to `NaN`?
* other constructors: `n` and `p` instead of `.fam` and `.bim`, take an existing matrix `X`
* Base.show to display missing value as `NA`?
* Tests

#### Bugs/Known Issues
* indexing by array of rownames and/or colnames fails:

```julia
julia> bed[["per22_per22", "per45_per45"], 2]
ERROR: ArgumentError: unable to check bounds for indices of type String
 in checkindex at ./abstractarray.jl:378 [inlined]
 in checkbounds_indices at ./abstractarray.jl:309 [inlined]
 in checkbounds at ./abstractarray.jl:270 [inlined]
 in checkbounds at ./abstractarray.jl:284 [inlined]
 in _getindex at ./multidimensional.jl:270 [inlined]
 in getindex(::BEDMatrices.BEDMatrix{UInt8,Array{UInt8,2}}, ::Array{String,1}, ::Int64) at ./abstractarray.jl:752
```

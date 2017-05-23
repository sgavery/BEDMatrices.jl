
#### To Implement

* add rows argument for `BEDdot`
* make column tool interface consistent with StatsBase.jl?
* `column_any` and `column_all`?
* `column_moment`?
* use tuple of tuples for BEDMatrix._bytemap?
* write benchmarks on simulated `BEDMatrix`
* Benchmark usage of `tocontiguous`
* use mean imputation instead of setting NA to 0?
* other constructors: take an existing matrix `X`?
* Test exceptions/invalid input
* note: one could have a LinkedMatrix as the `X` of a `BEDMatrix`
* output (serially) to high-performance (or other) format
* write bed files?
* read all information in .fam and .bim files?


#### Major Design Decision: NA values

While returning a special value, `NA_byte`, for missing values is a
faithful representation of the underlying `BED` format, this is
problematic. Specifically, it requires extreme diligence on users of
`BEDMatrices` to not (silently) make incorrect calculations.

Possible Solutions:
1. Require `BEDMatrices` to be formed only for types with `NaN`:
   `Float16`, `Float32`, `Float64`, and `BigFloat` (`nan`).
2. Use `DataArrays` approach (current backend of `DataFrames`)
3. Use `NullableArrays` approach and `Nullable`s (current backend of `DataTables`)
4. Custom implementation

The last option is ridiculously out of scope; this leaves the first
three. Note that there seems to be a bit of a divide in the Stats
julia community between the second and third approaches. (Hence the
schism between `DataFrames` and `DataTables`. See
eg. [this discussion](https://discourse.julialang.org/t/datatables-or-dataframes/3160/15).)

Advantages of the `Nullable` approach:
1. `Nullable` is in julia base, so `NullableArrays` seems more native julian; unlike `DataArrays.NA`
2. `NullableArrays` are (currently) more performant, which will likely matter with large .bed files
3. As compared with `NaN` approach, allows one to work with integers.
4. `NullableArrays` seem to have more support from julia's core developers.

On the other hand, `NullableArrays` and `Nullable`s are more awkward
to use; `NullableArrays` are currently less well-developed. So this
likely requires implementing any standard statistics tools like
`mean`, `std`, `quantile`, etc. that might be required; although that
was likely to be true anyway.

The `NaN` approach is the simplest and computationally fastest
approach; however, it means one cannot work with integers. This may
not be that big a deal, since basically any calculation requires
converting to floats at some point anyway.

## Implementation
* make `NA` behavior customizable, with `navalue` field
* relegate the question to "outer constructors"


#### Bugs/Known Issues

* Cannot seem to overload `Base.dot` for `BEDColumns`, hence using `BEDdot`.

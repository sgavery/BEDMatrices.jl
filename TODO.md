
#### To Implement

* when flipped SNPs should have the column name changed

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
* Test Union treatment of `NA` values


## Implementation
* make `NA` behavior customizable, with `navalue` field
* relegate the question to "outer constructors"


#### Bugs/Known Issues

* Cannot seem to overload `Base.dot` for `BEDColumns`, hence using `BEDdot`.

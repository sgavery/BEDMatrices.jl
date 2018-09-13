BEDMatrices.jl
==============

Tools for efficiently reading, memory-mapping, and
manipulating [PLINK](http://zzz.bwh.harvard.edu/plink/) (see
also [PLINK1.9](https://www.cog-genomics.org/plink2)) BED-formatted
genotype data (.bed, .fam, .bim files)
in [julia](http://julialang.org/). This package is, in part, based on the R
package [BEDMatrix](https://github.com/QuantGen/BEDMatrix).


Installation
-----------

Note that there are SSL issues with julia v0.5 and earlier.

### From github

To install run

```julia
julia> Pkg.clone("https://github.com/sgavery/BEDMatrices.jl.git")
```

### From MSU gitlab (deprecated)

The package is available through
MSU [gitlab](https://gitlab.msu.edu/QuantGen/BEDMatrices.jl). To
install, run the following

```julia
julia> Pkg.clone("https://[username]@gitlab.msu.edu/QuantGen/BEDMatrices.jl.git")
```

You will be prompted for your gitlab credentials.


Features
--------
* Read plink bed files into standard `Matrix{T}` (use `BEDintomatrix`)
* Memory-mapped `BEDMatrix` type with efficient indexing
* `BEDMatrix` can be indexed with `rownames` and `colnames` read from
  accompanying .fam and .bim files
* RAW-formatted output
* Customizable NA behavior
* Tools for efficient column/SNP-wise calculations. (Parallelized
  multi-column operations will be provided elsewhere, along with
  regression and other statistical analysis tools.)


Examples/Basic Usage
-------------------

If we have a .bed file, "example.bed" with accompanying "example.fam"
and "example.bim" files in the current working directory, then we may
create a `BEDMatrix` object over a memory-mapping of "example.bed" as
follows:

```julia
julia> using BEDMatrices

julia> bed = BEDMatrix("example.bed");
```

For most purposes, we may now treat `bed` as an ordinary matrix. For
example, we may index into it as

```julia
julia> bed[2:12, 1:10]
11×10 Array{Int8,2}:
 1  1  1  1  3  2  2  2  1  1
 1  0  0  2  0  0  1  2  0  1
 2  0  0  0  1  0  2  1  1  2
 0  1  0  0  0  1  1  0  1  0
 1  1  1  0  0  0  0  2  1  1
 1  0  2  0  3  0  1  2  3  0
 1  2  2  0  1  2  1  0  2  0
 1  1  0  1  0  1  1  1  0  1
 1  2  1  1  2  0  1  1  0  1
 2  1  1  0  1  0  1  0  2  0
 2  0  0  1  1  2  0  1  0  1
```

If you prefer to use a different numeric type like `Float64`, you can
define `bed` as follows

```julia
julia> bed = BEDMatrix("example", datatype=Float64);

julia> bed[2:12, 1:10]
11×10 Array{Float64,2}:
 1.0  1.0  1.0  1.0  3.0  2.0  2.0  2.0  1.0  1.0
 1.0  0.0  0.0  2.0  0.0  0.0  1.0  2.0  0.0  1.0
 2.0  0.0  0.0  0.0  1.0  0.0  2.0  1.0  1.0  2.0
 0.0  1.0  0.0  0.0  0.0  1.0  1.0  0.0  1.0  0.0
 1.0  1.0  1.0  0.0  0.0  0.0  0.0  2.0  1.0  1.0
 1.0  0.0  2.0  0.0  3.0  0.0  1.0  2.0  3.0  0.0
 1.0  2.0  2.0  0.0  1.0  2.0  1.0  0.0  2.0  0.0
 1.0  1.0  0.0  1.0  0.0  1.0  1.0  1.0  0.0  1.0
 1.0  2.0  1.0  1.0  2.0  0.0  1.0  1.0  0.0  1.0
 2.0  1.0  1.0  0.0  1.0  0.0  1.0  0.0  2.0  0.0
 2.0  0.0  0.0  1.0  1.0  2.0  0.0  1.0  0.0  1.0
```

Note that `3.0` (or `3` in the previous example) currently is the
default indicator for missing values; this may change. (While this is
in many ways an unfortunate choice, it is a literal translation of the
BED format.) This is set by `BEDMatrices.NA_byte`:

```julia
julia> BEDMatrices.NA_byte
0x03
```

The `BEDMatrix` may be created with different NA behavior. For
example, when working with `Float`s it is probably more desirable to
use `NaN`:

```julia
julia> bed = BEDMatrix("example", datatype=Float64, navalue=NaN);

julia> bed[2:12, 1:10]
11×10 Array{Float64,2}:
 1.0  1.0  1.0  1.0  NaN    2.0  2.0  2.0    1.0  1.0
 1.0  0.0  0.0  2.0    0.0  0.0  1.0  2.0    0.0  1.0
 2.0  0.0  0.0  0.0    1.0  0.0  2.0  1.0    1.0  2.0
 0.0  1.0  0.0  0.0    0.0  1.0  1.0  0.0    1.0  0.0
 1.0  1.0  1.0  0.0    0.0  0.0  0.0  2.0    1.0  1.0
 1.0  0.0  2.0  0.0  NaN    0.0  1.0  2.0  NaN    0.0
 1.0  2.0  2.0  0.0    1.0  2.0  1.0  0.0    2.0  0.0
 1.0  1.0  0.0  1.0    0.0  1.0  1.0  1.0    0.0  1.0
 1.0  2.0  1.0  1.0    2.0  0.0  1.0  1.0    0.0  1.0
 2.0  1.0  1.0  0.0    1.0  0.0  1.0  0.0    2.0  0.0
 2.0  0.0  0.0  1.0    1.0  2.0  0.0  1.0    0.0  1.0

julia> NArep(bed)
NaN
```

`NArep` returns the representation of missing values for the
`BEDMatrix`. One can also work with `Nullable`s
(see
[julia docs](https://docs.julialang.org/en/stable/stdlib/base/#nullables))
in this way:

```julia
julia> bed = BEDMatrix("example.bed", datatype=Nullable{Int}, navalue=Nullable{Int}());

julia> NArep(bed)
Nullable{Int64}()

julia> bed[2:12, 1:10]
11×10 Array{Nullable{Int64},2}:
 1  1  1  1  #NULL  2  2  2  1      1
 1  0  0  2  0      0  1  2  0      1
 2  0  0  0  1      0  2  1  1      2
 0  1  0  0  0      1  1  0  1      0
 1  1  1  0  0      0  0  2  1      1
 1  0  2  0  #NULL  0  1  2  #NULL  0
 1  2  2  0  1      2  1  0  2      0
 1  1  0  1  0      1  1  1  0      1
 1  2  1  1  2      0  1  1  0      1
 2  1  1  0  1      0  1  0  2      0
 2  0  0  1  1      2  0  1  0      1
```

Note that it may be preferable to work
with [NullableArrays](https://github.com/JuliaStats/NullableArrays.jl)
instead of the above slicings for computationally intensive
calculations. We leave that to the user or another module to
implement.

The `eltype` is exposed as the first parameter:

```julia
julia> typeof(bed)
BEDMatrices.BEDMatrix{Nullable{Int64},Array{UInt8,2}}

julia> eltype(bed)
Nullable{Int64}
```

The second parameter, `Array{UInt8,2}`, refers to the internal
byte-level BED representation.

One can also get basic metadata about the `BEDMatrix` with various
functions:

```julia
julia> bed = BEDMatrix("example.bed", datatype=Float64, navalue=NaN);

julia> path(bed)
"[...]/example.bed"

julia> size(bed)
(50,1000)

julia> sizeof(bed)
400000

julia> Base.summarysize(bed)
48615

julia> rownames(bed)[1:10]
10-element Array{String,1}:
 "per0_per0"
 "per1_per1"
 "per2_per2"
 "per3_per3"
 "per4_per4"
 "per5_per5"
 "per6_per6"
 "per7_per7"
 "per8_per8"
 "per9_per9"

julia> colnames(bed)[1:10]
10-element Array{String,1}:
 "snp0_A"
 "snp1_C"
 "snp2_G"
 "snp3_G"
 "snp4_G"
 "snp5_T"
 "snp6_G"
 "snp7_A"
 "snp8_A"
 "snp9_C"
```

The row and column names may be used for indexing:

```julia
julia> bed["per0_per0", 1]
0.0

julia> bed["per0_per0", "snp0_A"]
0.0
```

See help on `BEDMatrix` for information on handling .fam and .bim
files that are missing or have different names/paths. There is also
information on how to flip the major-minor allele reprentation of
selected SNPs. For the latter, see also `setflips!`. One can also use
a custom representation of SNP data.

```julia
help?> BEDMatrix
[...]
```

Column-wise Operations
----------------------

In general, working with slices, while convenient and intuitive, is
bad for performance; each slice has to allocate memory. Additionally,
because of the compact BED format, there are (roughly four-fold)
performance gains to be had by working at the byte level, _if_ one
works with large contiguous subsets of columns.

There are several tools for efficiently performing basic functions on
columns of a `BEDMatrix`. Here are some examples:

```julia
julia> bed = BEDMatrix("example.bed");

julia> hasNAs(view(bed, :, 2))
false

julia> hasNAs(view(bed, :, 3))
true

julia> hasNAs(bed)
true

julia> countNAs(view(bed, :, 2))
0

julia> countNAs(view(bed, :, 3))
2

julia> countNAs(bed)
975

julia> column_dist(bed, 2)  # number of (0s, 1s, 2s, NAs)
(13,26,11,0)

julia> column_dist(bed, 3)
(16,25,7,2)

julia> column_sum(bed, 2)
48

julia> column_sum(bed, 3)
39

julia> column_sum(x -> (x - 0.5)^5, bed, 3)
53.4375

julia> column_norm(bed, 3)
7.280109889280518

julia> column_norm(bed, 3, :, 1)
39.0

julia> column_dot(bed, 1, 2)
41

julia> column_sumabs2(bed, 3)
53

julia> column_mean_and_std(bed, 3)
(0.8125,0.6733924909059431)

```

As demonstrated above, many functions support using `SubArray`s
(formed with `view`) for convenience. While these methods are all
faster and more memory efficient than working with slices, there is
some overhead. For optimal performance use the `column_<methodname>`
form, that these functions call. See the documentation for usage and
method details, including treatment of missing values.

Of these, the most important is `column_dist(bed, col, rows)`, which
returns the column distribution: the number of 0s, 1s, 2s, and NAs in
`B[rows, col]`. This is generally the complete set of information that
one needs about an individual SNP on (a subset of) the sample
population; indeed this is used by many of the other methods.


BED format details
------------------

See [PLINK1 description](http://zzz.bwh.harvard.edu/plink/binary.shtml) and [PLINK2 description](https://www.cog-genomics.org/plink2/formats#bed) for details about the format.


Other julia packages supporting BED format in some capacity
-----------------------------------------------------------

* [PLINK.jl](https://github.com/klkeys/PLINK.jl)
* [OpenMendel](https://openmendel.github.io/), see [SnpArrays.jl](https://openmendel.github.io/SnpArrays.jl/latest/)
* [JWAS.jl](https://github.com/reworkhow/JWAS.jl)
* [StatGenData.jl](https://github.com/dmbates/StatGenData.jl)
* [VarianceComponentTest.jl](https://github.com/Tao-Hu/VarianceComponentTest.jl)
* [Bio.jl](http://biojulia.net/Bio.jl/latest/)

#### See [Julia.jl](https://github.com/svaksha/Julia.jl/blob/master/Biology.md#genomics) for a listing of genomics tools in julia.

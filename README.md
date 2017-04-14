BEDMatrices.jl
==============

Tools for efficiently reading, memory-mapping, and
manipulating [PLINK](http://zzz.bwh.harvard.edu/plink/) (see
also [PLINK1.9](https://www.cog-genomics.org/plink2)) BED formatted
genotype data (.bed, .fam, .bim files)
in [julia](http://julialang.org/). This package is, in part, based on the R
package [BEDMatrix](https://github.com/QuantGen/BEDMatrix).

Features
--------

Example/Basic Usage
-------------------

If we have a .bed file, "example.bed" with accompanying "example.fam"
and "example.bim" files, in the current working directory, then we may
load and create a `BEDMatrix` object as follows:

```julia
julia> include("BEDMatrices.jl");
julia> bed = BEDMatrices.BEDMatrix("example.bed");
```

For most purposes, we may now treat `bed` as an ordinary matrix. For
example, we may index into it as

```julia
julia> bed[2:12, 1:10]
11×10 Array{UInt8,2}:
 0x01  0x01  0x01  0x01  0x03  0x02  0x02  0x02  0x01  0x01
 0x01  0x00  0x00  0x02  0x00  0x00  0x01  0x02  0x00  0x01
 0x02  0x00  0x00  0x00  0x01  0x00  0x02  0x01  0x01  0x02
 0x00  0x01  0x00  0x00  0x00  0x01  0x01  0x00  0x01  0x00
 0x01  0x01  0x01  0x00  0x00  0x00  0x00  0x02  0x01  0x01
 0x01  0x00  0x02  0x00  0x03  0x00  0x01  0x02  0x03  0x00
 0x01  0x02  0x02  0x00  0x01  0x02  0x01  0x00  0x02  0x00
 0x01  0x01  0x00  0x01  0x00  0x01  0x01  0x01  0x00  0x01
 0x01  0x02  0x01  0x01  0x02  0x00  0x01  0x01  0x00  0x01
 0x02  0x01  0x01  0x00  0x01  0x00  0x01  0x00  0x02  0x00
 0x02  0x00  0x00  0x01  0x01  0x02  0x00  0x01  0x00  0x01
```

If you prefer to use a different numeric type like `Float64`, you can
define `bed` as follows

```julia
julia> bed = BEDMatrices.BEDMatrix("example.bed", Float64);
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

Note that `3.0` (or `0x03` in the previous example) currently
indicates missing value; this will likely change. This is set by
`BEDMatrices.NA_byte`:
```julia
julia> BEDMatrices.NA_byte
0x03
```

The return type is exposed as the first parameter:

```
julia> typeof(bed)
BEDMatrices.BEDMatrix{Float64,Array{UInt8,2}}
```

The second parameter `Array{UInt8,2}` refers to the internal .bed
representation, and will generally be the same as above.

More Documentation
------------------

BED format details
------------------

See [PLINK1 description](http://zzz.bwh.harvard.edu/plink/binary.shtml) and [PLINK2 description](https://www.cog-genomics.org/plink2/formats#bed) for details about the format.

Other julia packages supporting BED format in some capacity
-----------------------------------------------------------
* [JWAS.jl](https://github.com/reworkhow/JWAS.jl)
* [StatGenData.jl](https://github.com/dmbates/StatGenData.jl)
* [VarianceComponentTest.jl](https://github.com/Tao-Hu/VarianceComponentTest.jl)

#### See [Julia.jl](https://github.com/svaksha/Julia.jl/blob/master/Biology.md#genomics) for a listing of genomics tools in julia.

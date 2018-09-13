

######################### Basic BED file tools #########################

## .bed file format ##################################################
#
# See http://pngu.mgh.harvard.edu/~purcell/plink/binary.shtml
# http://zzz.bwh.harvard.edu/plink/binary.shtml
# and https://www.cog-genomics.org/plink2/formats#bed
#
######################################################################

"""
    rawformat(snp::Integer, quartermap=quarterstohuman)

Returns RAW format from the BED format quarter-byte:

| BED    | RAW (default) | general         | Meaning                 |
|:------ |:------------- |:--------------- |:----------------------- |
| `0b00` | `0b10`        | `quartermap[1]` | homozygous minor allele |
| `0b11` | `0b00`        | `quartermap[2]` | homozygous major allele |
| `0b01` | `NA_byte`     | `quartermap[3]` | missing value           |
| `0b10` | `0b01`        | `quartermap[4]` | heterozygous            |

"""
@inline function rawformat(snp::Integer, quartermap=Consts.quarterstohuman)
    @inbounds return quartermap[snp + 1]
end

"""
    breakbyte(byte::UInt8, bytemap=Consts.defaultbytemap)

Return a `Tuple{Vararg{UInt8, 4}}` of the RAW-formatted SNP quarters, as
determined by `bytemap`, in `byte`.

"""
@inline function breakbyte(byte::UInt8, bytemap=Consts.defaultbytemap, flip::Bool=false)
    @inbounds return bytemap[byte + 1, flip + 1]
end

@inline function breakbyte(byte::UInt8, bytemap, snpind::Int, flip::Bool=false)
    @inbounds return breakbyte(byte, bytemap, flip)[snpind]
end

# Approach from Quantgen/BEDMatrix: shift by 2(n - 1) and then mask
# result.
#
# Performance-wise, the lookup table approach of breakbyte is a bit
# faster.
"""
    quarter(byte::UInt8, n::Integer, quartermap=quarterstohuman)

Returns the RAW-formatted `n`th quarter of `byte`, equivalent to
`breakbyte(byte)[n]`.

"""
function quarter(byte::UInt8, n::Integer, quartermap=Consts.quarterstohuman)
    @inbounds return quartermap[((byte >> 2(n - 1)) & 0b00000011) + 1]
end

@inline function flipBEDbyte(byte::UInt8)
    # 0b00 -> 0b11
    # 0b11 -> 0b00
    # 0b10 -> 0b10
    # 0b01 -> 0b01

    @inbounds return Consts.flipBEDbytes[byte + 1]
end

"""
    BEDsize(filebase::AbstractString)

Returns the number of rows and columns (as a tuple) of the BED matrix
in `filebase*".bed"` as determined by linecounts of `filebase*".fam"`
and `filebase*".bim"`.

"""
function BEDsize(filebase::AbstractString)
    famfile, bimfile = filebase*".fam", filebase*".bim"
    isfile(famfile) || error("missing FAM file \"$famfile\"")
    isfile(bimfile) || error("missing BIM file \"$bimfile\"")

    n = countlines(famfile)
    p = countlines(bimfile)
    return n, p
end

"""
    checkmagic(bytes)

Determine if first two "magic" bytes match plink format. If not throws
an error, else returns `true`.

If bytes are UInt16, then it should be in big endian format.

"""
# Puts first two bytes into big endian format
checkmagic(bytes::Vector{UInt8}) = checkmagic((convert(UInt16, bytes[1]) << 8) + bytes[2])

checkmagic(twobytes::UInt16) = (twobytes == Consts.plinkmagic ||
                                error("Bad magic: not a plink bed file"))

function checkmagic(bedstream::IO)
    seekstart(bedstream)

    # Note that I don't just read(bedstream, UInt16) to avoid system
    # endian dependence, although I believe julia currently only runs
    # on little endian machines.
    return checkmagic((convert(UInt16, read(bedstream, UInt8)) << 8) + read(bedstream, UInt8))
end


"""
    BEDmode(byte)

Returns the BED mode as determined by the third byte of the file:
either the current standard, `:SNPmajor`, or older plink formats,
`:SNPminor`.

"""
BEDmode(byte::UInt8) = Consts.modes[byte]

function BEDmode(bedstream::IO)
    seek(bedstream, 2)
    return BEDmode(read(bedstream, UInt8))
end

BEDmode(bytevector::Vector{UInt8}) = BEDmode(bytevector[3])


getbytemap(navalue::T) where {T} = getbytemap((convert(T, 0b10), navalue,
                                     convert(T, 0b01), convert(T, 0b00)))
function getbytemap(quartermap::Tuple{T, T, T, T}) where {T}
    bytemap = hcat([(rawformat(byte & 0b00_00_00_11, quartermap),
                     rawformat((byte & 0b00_00_11_00) >> 2, quartermap),
                     rawformat((byte & 0b00_11_00_00) >> 4, quartermap),
                     rawformat((byte & 0b11_00_00_00) >> 6, quartermap)) for
                    byte in 0x00:0xff],
                   [(rawformat(byte & 0b00_00_00_11, quartermap),
                     rawformat((byte & 0b00_00_11_00) >> 2, quartermap),
                     rawformat((byte & 0b00_11_00_00) >> 4, quartermap),
                     rawformat((byte & 0b11_00_00_00) >> 6, quartermap)) for
                    byte in Consts.flipBEDbytes])

    bytemap
end


"""
    unsafe_breakbyte!(vector::AbstractArray, byte::UInt8,
                      bytemap=bytetoquarters, vecstart=1,
                      quarterstart=1, quarterstop=4,
                      flip::Bool=false)

Fills `vector` starting at `vecstart` with the `quarterstart`th snp
through the `quarterstop`th snp in `byte`, as determined by `bytemap`.

## (Unchecked) Constraints:
* 1 <= `quarterstart` <= `quarterstop` <= 4
* 1 <= `vecstart` <= `end - (quarterstop - quarterstart)`

"""
function unsafe_breakbyte!(vector::AbstractArray, byte::UInt8,
                           bytemap=Consts.bytetoquarters, flip::Bool=false,
                           vecstart::Integer=1, quarterstart::Integer=1, quarterstop::Integer=4)

    # Use `unsafe_copyto!`?
    @inbounds copyto!(vector, vecstart,
                             breakbyte(byte, bytemap, flip), quarterstart,
                             quarterstop - quarterstart + 1)
    vector
end


"""
    unsafe_copybytestosnps!(snparray::AbstractArray, bytearray::AbstractArray{UInt8},
                            bytestart::Integer, quarterstart::Integer,
                            bytestop::Integer, quarterstop::Integer,
                            deststart::Integer=1, bytemap=bytetoquarters,
                            flip::Bool=false)

Fills `snparray[deststart:(deststart + num_snps)]` with snps, where

    num_snps = 4*(bytestop - bytestart) + (quarterstop - quarterstart) + 1

The SNPs start with the `quarterstart`th snp in the `bytestart`th
byte, and end with `quarterstop` and `bytestop`.  The SNP
representation is determined by the `bytemap`. Currently only supports
unit stride for both `snparray` and `bytearray`, since at present
there is no compelling use-case involving non-unit strides.

Utility function used by more front-facing functions.

If `flip` is `true`, then minor and major allele are flipped when copying.

"""
function unsafe_copybytestosnps!(snparray::AbstractArray, bytearray::AbstractArray{UInt8},
                                 bytestart::Integer, quarterstart::Integer,
                                 bytestop::Integer, quarterstop::Integer,
                                 deststart::Integer=1, bytemap=Consts.bytetoquarters,
                                 flip::Bool=false)

    # First byte
    if quarterstart != 1
        stop = ifelse(bytestart == bytestop, quarterstop, 4)
        @inbounds unsafe_breakbyte!(snparray, bytearray[bytestart], bytemap, flip,
                                    deststart, quarterstart, stop)
        deststart += stop - quarterstart + 1
        bytestart += 1
    end

    if bytestart <= bytestop
        # Last byte
        if quarterstop != 4
            @inbounds unsafe_breakbyte!(snparray, bytearray[bytestop], bytemap, flip,
                                        deststart + 4*(bytestop - bytestart), 1, quarterstop)
            bytestop -= 1
        end

        # Main course
        @simd for bytej in bytestart:bytestop
            @inbounds unsafe_breakbyte!(snparray, bytearray[bytej], bytemap, flip,
                                        deststart + 4*(bytej - bytestart))
        end
    end

    return snparray
end


######################### Reading .bed into a Matrix #########################

"""
    BEDintomatrix{T<:Real}(bedfilename::AbstractString; datatype::Type{T}=UInt8,
                           navalue=NA_byte, n::Integer=0, p::Integer=0, use_mmap=true)

Returns a `Matrix{T}` representation of the BED file `bedfilename`. If
the number of rows or columns, `n` and `p`, are not provided, then
they will attempt to be determined from corresponding .fam and .bim
files. `use_mmap` determines whether to memory map the bedfile and
then read into the matrix for potential speedups. Missing values are
set to `navalue`.

"""
function BEDintomatrix(bedfilename::AbstractString;
                       datatype::Type{T}=UInt8,
                       navalue=NA_byte,
                       n::Integer=0,
                       p::Integer=0,
                       use_mmap::Bool=true) where {T<:Real}
    if !isfile(bedfilename)
        isfile(bedfilename*".bed") || error("Cannot find file \"$bedfilename\"")
        bedfilename = bedfilename*".bed"
    end
    filebase = splitext(bedfilename)[1]

    if n == p == 0
        n, p = BEDsize(filebase)
    end

    A = Matrix{T}(undef, n, p)

    fsz = filesize(bedfilename)
    if use_mmap && 0 < fsz < typemax(Int)
        bedfilevector = open(bedfilename, "r") do bedfile
            bvec = Mmap.mmap(bedfile, Vector{UInt8}, (Int(fsz),))
            checkmagic(bvec)
            BEDmode(bvec) == :SNPmajor || error("SNPminor mode not supported")
            bvec
        end
        BEDintomatrix!(A, bedfilevector, convert(T, navalue))
    else
        open(bedfilename, "r") do bedfin
            checkmagic(bedfin)
            BEDmode(bedfin) == :SNPmajor || error("SNPminor mode not supported")

            BEDintomatrix!(A, bedfin, convert(T, navalue))
        end
    end

    return A::Matrix{T}
end

function BEDintomatrix!(A::AbstractMatrix{T}, bedvector::Vector{UInt8}, navalue=NA_byte) where {T<:Real}
    n, p = size(A)
    byteheight = ceil(Int, n/4)
    quarterstop = n - 4*(byteheight - 1)

    bytemap = getbytemap(convert(T, navalue))

    @inbounds for col in 1:p
        unsafe_copybytestosnps!(A, bedvector, byteheight*(col - 1) + 3 + 1, 1, byteheight*col + 3, quarterstop, n*(col - 1) + 1, bytemap)
    end

    return A
end

"""
    BEDintomatrix!{T<:Real}(A::AbstractMatrix{T}, bedstream::IO, navalue=NA_byte)

Fills `A` with type `T` representation of `bedstream` that corresponds
to .bed file format. `A` must have correct dimensions, as determined
via `BEDMatrices.BEDsize`, for example.

"""
function BEDintomatrix!(A::AbstractMatrix{T}, bedstream::IO, navalue=NA_byte) where {T<:Real}
    n, p = size(A)
    bytestop = ceil(Int, n/4)
    quarterstop = n - 4*(bytestop - 1)

    bytemap = getbytemap(convert(T, navalue))

    bytecol = Vector{UInt8}(bytestop)
    @inbounds for col in 1:p
        read!(bedstream, bytecol)
        unsafe_copybytestosnps!(A, bytecol, 1, 1, bytestop, quarterstop, n*(col - 1) + 1, bytemap)
    end

    return A
end


######################### BEDMatrix type #########################
# See implementation of BitArray for inspiration julia/base/bitarray.jl
# in /usr/share/julia

## Things we get for "free", but may want to overload for better performance, etc.:
# * checkbounds(B, i...)
# * eltype(B)
# * length(B)
# * ndims(B)
# * strides(B)
# * view(B, i...) -- not sure if it is memory and computationally efficient...
# * matrix operations: v * B and B * v, norm(B), etc.


# Allowing for `X` to be a general subtype `S` of `AbstractMatrix`
# requires exposing the type as below. See discussion at
# http://docs.julialang.org/en/stable/manual/performance-tips/#type-declarations

struct BEDMatrix{T, S<:AbstractMatrix} <: DenseArray{T, 2}
    n::Int                               # number of samples (rows)
    p::Int                               # number of SNPs (columns)
    X::S                                 # Matrix{UInt8} bed file
    navalue::T                           # representation of missing

    path::String                         # file location
    colnames::Vector{String}             # SNP labels, can be used for indexing
    rownames::Vector{String}             # sample labels, can be used for indexing
    colnamemap::Dict{String, Int}        # for quick name indexing by column
    rownamemap::Dict{String, Int}        # for quick name indexing by row

    _byteheight::Int                     # number of bytes in each column
    _lastrowSNPheight::Int               # number of SNPs in last byte of each column

    _bytemap::Matrix{Tuple{T, T, T, T}}  # quarters and flipped quarters for 0x00:0xff

    _flip::BitVector                     # whether to flip SNP major--minor
                                         # allele encoding for each SNP

    function BEDMatrix{T, S}(n::Integer, p::Integer, X::S, navalue::T,
                       path::AbstractString, colnames::AbstractVector, rownames::AbstractVector,
                       bytemap, flip::BitVector) where {T, S<:AbstractMatrix{UInt8}}
        byteheight = ceil(Int, n/4)
        lastrowheight = n - 4*(byteheight - 1)

        size(X) == (byteheight, p) || throw(DimensionMismatch("Matrix dimensions $(size(X)) do not agree with supplied BED dimensions (n = $n, p = $p)"))
        length(colnames) == p || throw(DimensionMismatch("colnames has incorrect length"))
        length(rownames) == n || throw(DimensionMismatch("rownames has incorrect length"))
        size(bytemap) == (256, 2) || error("bytemap must be of size (256, 2)")
        length(flip) == p || throw(DimensionMismatch("flip has length $(length(flip)); expected length of p = $p"))

        colnamemap = Dict{String, Int}()
        for (cidx, colname) in enumerate(colnames)
            colnamemap[colname] = cidx
        end

        rownamemap = Dict{String, Int}()
        for (ridx, rowname) in enumerate(rownames)
            rownamemap[rowname] = ridx
        end


        return new{T, S}(n, p, X, navalue, path, colnames, rownames, colnamemap, rownamemap,
                         byteheight, lastrowheight, bytemap, flip)
    end
end

"""
    BEDMatrix(bedfilename::AbstractString;
              datatype::DataType=Int8, nsamples::Integer=0, nSNPs::Integer=0, navalue=NA_byte
              famfile::AbstractString="", bimfile::AbstractString="",
              quartermap::Tuple=Consts.quarterstohuman, flip::AbstractVector=[])

Create a `BEDMatrix` of type `datatype` using memory mapping of BED
file, `bedfilename`. Use `navalue` for missing values. `famfile` and
`bimfile` are only required if the .fam and .bim file are not in same
location as or do not share the same base name as the
`bedfilename`. `nsamples` and `nSNPs` are only used if .bim and .fam
file cannot be found and are not provided; in which case the
rownames/colnames will be generic.

`convert(datatype, x)` must work for `x` in `[0b00, 0b01, 0b10, navalue]`.

## Optional arguments
 - `datatype`: the return type for indexing.
 - `nsamples`: number of samples if .fam file cannot be found
 - `nSNPs`: the number of SNPs if .bim file cannot be found
 - `navalue`: the representation of missing values; this
   should not be used with quartermap.
 - `famfile`: the .fam file, only needed if the .fam file does
   not have the same base name as `bedfilename`.
 - `bimfile`: the .bim file, only needed if the .bim file does
   not have the same base name as `bedfilename`.
 - `quartermap`: a `Tuple` of the respective representation of
   homozygous minor, homozygous major, missing, and heterozygous
   values. Useful if you don't want to use standard RAW format
   encoding.
 - `flip`: a `Vector` of indices that should have the major--minor
   representation flipped. `flip` may either be a list of indices, a
   logical vector with `false` representing no flipping, or a list of
   columnnames. This may be recovered or changed later using
   `getflips` and `setflips!`, respectively.

## Examples
```julia
julia> bed = BEDMatrix("test/data/example.bed");

julia> bed[1:5, 1:5]
5×5 Array{Int8,2}:
 0  1  1  1  0
 1  1  1  1  3
 1  0  0  2  0
 2  0  0  0  1
 0  1  0  0  0

julia> bed = BEDMatrix("test/data/example", datatype=Nullable{Int}, navalue=Nullable{Int}());

julia> bed[1:5, 1:5]
5×5 Array{Nullable{Int64},2}:
 0  1  1  1  0    
 1  1  1  1  #NULL
 1  0  0  2  0    
 2  0  0  0  1    
 0  1  0  0  0    

julia> bed = BEDMatrix("test/data/example.bed", datatype=Float64, navalue=NaN);

julia> bed[1:5, 1:5]
5×5 Array{Float64,2}:
 0.0  1.0  1.0  1.0    0.0
 1.0  1.0  1.0  1.0  NaN  
 1.0  0.0  0.0  2.0    0.0
 2.0  0.0  0.0  0.0    1.0
 0.0  1.0  0.0  0.0    0.0

julia> bed = BEDMatrix("test/data/example", famfile="badfilename", nsamples=50);

julia> rownames(bed)[1:5]
5-element Array{String,1}:
 "sample_1"
 "sample_2"
 "sample_3"
 "sample_4"
 "sample_5"
```

"""
function BEDMatrix(bedfilename::AbstractString;
                   datatype::DataType=Int8,
                   nsamples::Integer=0, nSNPs::Integer=0, navalue=NA_byte,
                   famfile::AbstractString="", bimfile::AbstractString="",
                   quartermap::Tuple=Consts.quarterstohuman, flip::AbstractVector=Int[])

    if !isfile(bedfilename)
        isfile(bedfilename*".bed") || error("Cannot find file \"$bedfilename\"")
        bedfilename = bedfilename*".bed"
    end
    filebase = splitext(bedfilename)[1]

    rownames = create_rownames(filebase, famfile, nsamples)
    colnames = create_colnames(filebase, bimfile, nSNPs)
    n, p = length(rownames), length(colnames)

    X = open(filebase*".bed", "r") do bedfile
        checkmagic(bedfile)
        BEDmode(bedfile) == :SNPmajor || error("Old-style SNP-minor mode bed"*
                                               " files not supported. Use plink"*
                                               " to convert to SNP-major format")

        Mmap.mmap(bedfile, Matrix{UInt8}, (ceil(Int, n/4), p))
    end

    # validate quartermap and navalue
    if quartermap == Consts.quarterstohuman
        bytemap = getbytemap(convert(datatype, navalue))
    else
        bytemap = getbytemap(quartermap)
        navalue = quartermap[2]
    end

    # Process and validate flip
    flip = process_flip(flip, colnames)

    return BEDMatrix{datatype, typeof(X)}(n, p, X, convert(datatype, navalue),
                                          abspath(bedfilename), colnames, rownames, bytemap, flip)
end

function parsefamline(line)
    fields = split(line)

    (fields[1], fields[2], parse(Int, fields[3]), parse(Int, fields[4]), parse(Int, fields[5]), parse(Int, fields[6]))
end

famrowname(line) = join(split(line)[1:2], '_')

readrownames(famfile::AbstractString) = map(famrowname, readlines(famfile))

function create_rownames(filebase::AbstractString, famfile::AbstractString="", nsamples::Integer=0)
    famfile = famfile == "" ? filebase*".fam" : famfile
    if isfile(famfile)
        rownames = readrownames(famfile)
    elseif nsamples > 0
        rownames = [string("sample_", j) for j in 1:nsamples]
    else
        error("Cannot find FAM file \"$famfile\" and nsamples not provided")
    end

    return rownames
end

function bimcolname(line)
    fields = split(line)
    string(fields[2], '_', fields[5])
end

readcolnames(bimfile::AbstractString) = map(bimcolname, readlines(bimfile))

function create_colnames(filebase::AbstractString, bimfile::AbstractString="", nSNPs::Integer=0)
    bimfile = bimfile == "" ? filebase*".bim" : bimfile

    if isfile(bimfile)
        colnames = readcolnames(bimfile)
    elseif nSNPs > 0
        colnames = [string("SNP_", j) for j in 1:nSNPs]
    else
        error("Cannot find BIM file \"$bimfile\" and nSNPs not provided")
    end
end

function process_flip(flip, colnames)
    p = length(colnames)

    if eltype(flip) == Bool
        length(flip) == p || throw(DimensionMismatch("flip has length $(length(flip)) expected length of p = $p"))
        flip = BitVector(flip)
    elseif eltype(flip) <: Integer
        checkbounds(Bool, 1:p, flip) || error("Invalid index in flip vector")
        flip = BitVector(map(x -> x in flip, 1:p))
    elseif eltype(flip) <: AbstractString
        for name in flip
            name in colnames || error("\"$name\" in flip is not a column_name")
        end
        flip = BitVector(map(x -> x in flip, colnames))
    else
        error("Invalid flip vector: expected logical, index, or column name vector")
    end

    return flip
end


######################## Getters, Properties, and Setters ########################

Base.size(B::BEDMatrix) = (B.n, B.p)
function Base.size(B::BEDMatrix, k::Integer)
    k > 0 || Base.throw_boundserror(size(B), k)
    # This is the same behavior as existing Base.size methods
    return k == 1 ? B.n : ifelse(k == 2, B.p, 1)
end


Base.IndexStyle(::Type{T}) where {T<:BEDMatrix} = Base.IndexCartesian()

# not sure if this is the right definition; see
# https://github.com/JuliaLang/julia/issues/16614.
# This is the size of `B[:, :]`.
Base.sizeof(B::BEDMatrix{T, S}) where {T, S} = B.n*B.p*sizeof(T)

"""
    path(B::BEDMatrix)

Returns the location of the .bed file used to create `X`.

"""
path(B::BEDMatrix) = B.path

"""
    rownames(B::BEDMatrix)

Returns `Vector{String}` of row names. These may be used for indexing.

"""
rownames(B::BEDMatrix) = B.rownames

"""
    colnames(B::BEDMatrix)

Returns `Vector{String}` of column names. These may be used for indexing.

"""
colnames(B::BEDMatrix) = B.colnames

"""
    NArep(B::BEDMatrix)

Returns value used for missing entries.

"""
NArep(B::BEDMatrix) = B.navalue

getrow(B::BEDMatrix, rowname::AbstractString) = B.rownamemap[rowname]
getcol(B::BEDMatrix, colname::AbstractString) = B.colnamemap[colname]

function _get_quartermap_from_bytemap(bytemap::Matrix)
    # 0b11_10_01_00 is the BED encoding for (homozygous minor,
    # homozygous major, missing, heterozygous). Recall the reverse
    # (little endian-like) order of the SNPs.
    bytemap[0b11_10_01_00 + 1, false + 1]
end

"""
    getquartermap(bed::BEDMatrix)

Returns the `quartermap` for `bed`. This is a `Tuple` of four values
giving the representation of respectively (homozygous minor,
homozygous major, missing value, heterozygous).

"""
function getquartermap(bed::BEDMatrix)
    _get_quartermap_from_bytemap(bed._bytemap)
end

"""
    getflips(B::BEDMatrix)

Returns a `Vector{Int}` of the columns which have a flipped major--minor allele
representation. Use `setflips!` to change this.

"""
getflips(B::BEDMatrix) = find(B._flip)


"""
    setflips!(B::BEDMatrix, indices)

(Re)sets which columns of the BEDMatrix use a flipped major--minor
representation. For example, `setflips!(bed, [])` would make all columns use the
standard encoding and `setflips!(bed, indices(bed, 2))` would make all columns
use a flipped representation.

"""
function setflips!(B::BEDMatrix, flip_logical::AbstractVector{Bool})
    length(flip_logical) == size(B, 2) || throw(DimensionsMismatch("flip had incorrect length"))
    copyto!(B._flip, flip_logical)
    B._flip
end

function setflips!(B::BEDMatrix, flip_indices::AbstractVector{K}) where {K<:Integer}
    checkbounds(1:B.p, flip_indices)

    fill!(B._flip, false)

    @inbounds for col in flip_indices
        B._flip[col] = true
    end

    B._flip
end

function setflips!(B::BEDMatrix, flip_names::AbstractVector{S}) where {S<:AbstractString}
    # I don't do this in the same loop as below, since I want an
    # atomic operation.
    for name in flip
        name in colnames || error("\"$name\" in flip is not a column_name")
    end
    fill!(B._flip, false)

    for name in flip
        B._flip[getcol(B, name)] = true
    end

    B._flip
end


#################### Indexing ####################

Base.getindex(B::BEDMatrix, rownames::AbstractVector{T}, col) where {T<:AbstractString} = B[[getrow(B, rowname) for rowname in rownames], col]
Base.getindex(B::BEDMatrix, row, colnames::AbstractVector{T}) where {T<:AbstractString}= B[row, [getcol(B, colname) for colname in colnames]]
Base.getindex(B::BEDMatrix, rownames::AbstractVector{S}, colnames::AbstractVector{T}) where {T<:AbstractString, S<:AbstractString} = [B[row, col] for row in rownames, col in colnames]

Base.getindex(B::BEDMatrix, rowname::AbstractString, colname::AbstractString) = B[getrow(B, rowname), getcol(B, colname)]
Base.getindex(B::BEDMatrix, rowname::AbstractString, col) = B[getrow(B, rowname), col]
Base.getindex(B::BEDMatrix, row, colname::AbstractString) = B[row, getcol(B, colname)]

Base.getindex(B::BEDMatrix, rownames::AbstractVector{T}, colname::AbstractString) where {T<:AbstractString} = [B[row, colname] for row in rownames]
Base.getindex(B::BEDMatrix, rowname::AbstractString, colnames::AbstractVector{T}) where {T<:AbstractString} = [B[rowname, col] for col in colnames]

# This is is the only getindex method we _need_, the other methods are
# provided for better performance, or convenience.
function Base.getindex(B::BEDMatrix, row::Integer, col::Integer)
    @boundscheck checkbounds(B, row, col)

    unsafe_getindex(B, row, col)
end

function Base.getindex(B::BEDMatrix, ::Colon, col::Integer)
    @boundscheck checkbounds(B, :, col)

    unsafe_getcol(B, col)
end

function Base.getindex(B::BEDMatrix{T, S}, ::Colon, cols::AbstractVector{K}) where {T, S, K<:Integer}
    @boundscheck checkbounds(B, :, cols)

    matrix = Matrix{T}(undef, B.n, length(cols))

    for (cidx, col) in enumerate(cols)
        unsafe_getcol!(matrix, B, col, (cidx - 1)*B.n + 1)
    end
    return matrix
end

function Base.getindex(B::BEDMatrix, rrange::AbstractUnitRange, col::Integer)
    @boundscheck checkbounds(B, rrange, col)

    unsafe_getrowrange(B, rrange, col)
end

function Base.getindex(B::BEDMatrix{T, S}, rrange::AbstractUnitRange, cols::AbstractVector{K}) where {T, S, K <: Integer}
    @boundscheck checkbounds(B, rrange, cols)

    matrix = Matrix{T}(undef, length(rrange), length(cols))
    for (cidx, col) in enumerate(cols)
        unsafe_getrowrange!(matrix, B, rrange, col, (cidx - 1)*length(rrange) + 1)
    end
    return matrix
end

function unsafe_getindex(B::BEDMatrix{T, S}, row::Integer, col::Integer) where {T, S}
    byterow, snpind = rowtobytequarter(row)

    @inbounds return breakbyte(B.X[byterow, col], B._bytemap, snpind, B._flip[col])
end

function unsafe_getcol(B::BEDMatrix{T, S}, col::Integer) where {T, S}
    column = Vector{T}(undef, B.n)

    unsafe_getcol!(column, B, col)
    return column
end

function unsafe_getcol!(column::AbstractArray{T}, B::BEDMatrix{T, S}, col::Integer, deststart=1) where {T, S}
    @inbounds unsafe_copybytestosnps!(column, B.X, B._byteheight*(col - 1) + 1, 1, B._byteheight*col,
                                      B._lastrowSNPheight, deststart, B._bytemap, B._flip[col])
    return column
end

function unsafe_getrowrange(B::BEDMatrix{T, S}, rrange::AbstractUnitRange, col::Integer) where {T, S}
    vector = Vector{T}(undef, length(rrange))
    unsafe_getrowrange!(vector, B, rrange, col)
    return vector
end

function unsafe_getrowrange!(vector::AbstractArray{T}, B::BEDMatrix{T, S}, rrange::AbstractUnitRange, col::Integer, deststart=1) where {T, S}
    bytestart, quarterstart = rowtobytequarter(first(rrange))
    bytestart += B._byteheight*(col - 1)

    bytestop, quarterstop = rowtobytequarter(last(rrange))
    bytestop += B._byteheight*(col - 1)

    @inbounds unsafe_copybytestosnps!(vector, B.X, bytestart, quarterstart, bytestop, quarterstop, deststart, B._bytemap, B._flip[col])
    return vector
end

function rowtobytequarter(row::Integer)
    # Curse you julia for your foolish 1-indexing!
    return (((row - 1) >> 2) + 1, ((row - 1) & 3) + 1)
end

"""
    getquarterblock(B::BEDMatrix, row::Integer, col::Integer)

Returns `(fourquarters, snpind)`, where `fourquarters` is the block of
quarters that `B[row, col]` belongs to: `B[row, col] =
fourquarters[snpind]`.

"""
function getquarterblock(B::BEDMatrix, row::Integer, col::Integer)
    byterow, snpind = rowtobytequarter(row)

    @inbounds return (breakbyte(B.X[byterow, col], B._bytemap, B._flip[col]), byterow, snpind)
end

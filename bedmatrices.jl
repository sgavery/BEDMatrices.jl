

######################### Plink Constants #########################

## .bed file format ##################################################
#
# See http://pngu.mgh.harvard.edu/~purcell/plink/binary.shtml
# http://zzz.bwh.harvard.edu/plink/binary.shtml
# and https://www.cog-genomics.org/plink2/formats#bed
#
######################################################################

const plinkmagic = [0b0110_1100, 0b0001_1011]
const modes = Dict(0b0000_0001 => :SNPmajor,
                   0b0000_0000 => :SNPminor)

const NA_byte = 0b11  # not recommended, 0xff? 0x9?
const quarterstohuman = (0b10, NA_byte, 0b01, 0b00)

"""
    rawformat(snp::Integer)

Returns RAW format from the BED format quarter-byte:

| BED    | RAW       | Meaning                 |
|:------ |:--------- |:----------------------- |
| `0b00` | `0b10`    | homozygous minor allele |
| `0b11` | `0b00`    | homozygous major allele |
| `0b01` | `NA_byte` | missing value           |
| `0b10` | `0b01`    | heterozygous            |

"""
@inline function rawformat(snp::Integer)
    @inbounds return quarterstohuman[snp + 1]
end

"""
    const bytetoquarters::Vector{Vector{UInt8}}

A constant array storing all 256 bytes split into 4 quarters. Note
that the mapping is given by `bytetoquarters[BEDbyte + 1]` because of
julia's 1-indexing.

# Example
```julia
julia> BEDbyte = 0b11101001
0xe9
julia> BEDMatrices.bytetoquarters[BEDbyte + 1]
4-element Array{UInt8,1}:
 0x03
 0x01
 0x01
 0x00
```

This can be understood as breaking the byte into quarters, in reverse
order: `0b01`, `0b10`, `0b10`, `0b11`, and then changing from BED
format to `rawformat` with `0x03` representing missing value.

"""
const bytetoquarters = [[rawformat(snp1), rawformat(snp2), rawformat(snp3), rawformat(snp4)] for
                        snp4 in 0b00:0b11 for
                        snp3 in 0b00:0b11 for
                        snp2 in 0b00:0b11 for
                        snp1 in 0b00:0b11]

"""
    breakbyte(byte::UInt8)

Return length-4 `Vector{UInt8}` of the RAW-formatted SNP quarters in `byte`.

"""
function breakbyte(byte::UInt8)
    @inbounds quarters = bytetoquarters[byte + 1]
    return quarters
end

# Approach from Quantgen/BEDMatrix: shift by 2(n - 1) and then mask
# result.
#
# Performance-wise, there does not seem to be a clear advantage of one
# approach over the other.
"""
    quarter(byte::UInt8, n::Integer)

Returns the RAW-formatted `n`th quarter of `byte`, equivalent to
`breakbyte(byte)[n]`.

"""
quarter(byte::UInt8, n::Integer) = rawformat((byte >> 2(n - 1)) & 0b00000011)


"""
    BEDsize(filebase::AbstractString)

Returns the number of rows and columns (as a tuple) of the BED matrix
in `filebase*".bed"` as determined by linecounts of `filebase*".fam"`
and `filebase*".bim"`.

"""
function BEDsize(filebase::AbstractString)
    famfile, bimfile = filebase*".fam", filebase*".bim"
    isfile(famfile) || error("missing .fam file")
    isfile(bimfile) || error("missing .bim file")
    n = countlines(famfile)
    p = countlines(bimfile)
    return n, p
end

"""
    checkmagic(bytes)

Determine if first two "magic" bytes match plink format. If not throws
an error, else returns `true`.

"""
checkmagic(bytes::Vector{UInt8}) = (bytes[1:2] == plinkmagic || error("Bad magic"))

function checkmagic(bedstream::IO)
    seekstart(bedstream)
    return checkmagic(read(bedstream, 2))
end


"""
    BEDmode(byte)

Returns the BED mode as determined by the third byte of the file:
either the current standard, `:SNPmajor`, or older plink formats,
`:SNPminor`.

"""
BEDmode(byte::UInt8) = modes[byte]

function BEDmode(bedstream::IO)
    seek(bedstream, 2)
    return BEDmode(read(bedstream, UInt8))
end

BEDmode(bytevector::Vector{UInt8}) = BEDmode(bytevector[3])


"""
    unsafe_breakbyte!(vector::AbstractArray, byte::UInt8, vecstart=1, quarterstart=1, quarterstop=4)

Fills `vector` starting at `vecstart` with the `quarterstart`th snp
through the `quarterstop`th snp in `byte`.

## (Unchecked) Constraints:
* 1 <= `quarterstart` <= `quarterstop` <= 4
* 1 <= `vecstart` <= `end - (quarterstop - quarterstart)`

"""
function unsafe_breakbyte!(vector::AbstractArray, byte::UInt8, vecstart=1, quarterstart=1, quarterstop=4)
    @inbounds copy!(vector, vecstart,
                    breakbyte(byte), quarterstart,
                    quarterstop - quarterstart + 1)
    vector
end


"""
    unsafe_copybytestosnps!(snparray::AbstractArray, bytearray::AbstractArray{UInt8}, bytestart::Integer, quarterstart::Integer, bytestop::Integer, quarterstop::Integer, deststart::Integer=1)

Fills `snparray[deststart:(deststart + num_snps)]` with snps, where

    num_snps = 4*(bytestop - bytestart) + (quarterstop - quarterstart) + 1

The snps start with the `quarterstart`th snp in the `bytestart`th
byte, and end with `quarterstop` and `bytestop`.  Currently only
supports unit stride for both `snparray` and `bytearray`, since at
present there is no compelling use-case involving non-unit strides.

Utility function used by more front facing functions.

"""
function unsafe_copybytestosnps!(snparray::AbstractArray, bytearray::AbstractArray{UInt8},
                                 bytestart::Integer, quarterstart::Integer, bytestop::Integer, quarterstop::Integer, deststart::Integer=1)
    # First byte
    if quarterstart != 1
        stop = bytestart == bytestop ? quarterstop : 4
        @inbounds unsafe_breakbyte!(snparray, bytearray[bytestart], deststart, quarterstart, stop)
        deststart += stop - quarterstart + 1
        bytestart += 1
    end

    if bytestart <= bytestop
        # Last byte
        if quarterstop != 4
            @inbounds unsafe_breakbyte!(snparray, bytearray[bytestop], deststart + 4*(bytestop - bytestart), 1, quarterstop)
            bytestop -= 1
        end

        # Main course
        @simd for bytej in bytestart:bytestop
            @inbounds unsafe_breakbyte!(snparray, bytearray[bytej], deststart + 4*(bytej - bytestart))
        end
    end

    return snparray
end


######################### Reading .bed into a Matrix #########################

"""
    BEDintomatrix{T<:Real}(bedfilename::AbstractString, ::Type{T}=UInt8, n::Integer=0, p::Integer=0; use_mmap=true)

Returns a `Matrix{T}` representation of the BED file `bedfilename`. If
the number or rows and columns, `n` and `p`, are not provided, then
they will attempt to be determined from corresponding .fam and .bim
files. `use_mmap` determines whether to memory map the bedfile and
then read into the matrix for potential speedups.

"""
function BEDintomatrix{T<:Real}(bedfilename::AbstractString, ::Type{T}=UInt8, n::Integer=0, p::Integer=0; use_mmap=true)
    n == p == 0 && !endswith(bedfilename, ".bed") && error("Need .bed file extension or dimensions n and p provided")
    if n == p == 0
        filebase = join((split(bedfilename, '.')[1:(end - 1)]), '.')
        n, p = BEDsize(filebase)
    end

    A = Matrix{T}(n, p)

    fsz = filesize(bedfilename)
    if use_mmap && 0 < fsz < typemax(Int)
        bedfilevector = Mmap.mmap(bedfilename, Vector{UInt8}, (Int(fsz),))
        checkmagic(bedfilevector)
        BEDmode(bedfilevector) == :SNPmajor || error("SNPminor mode not supported")
        BEDintomatrix!(A, bedfilevector)
    else
        open(bedfilename, "r") do bedfin
            checkmagic(bedfin)
            BEDmode(bedfin) == :SNPmajor || error("SNPminor mode not supported")

            BEDintomatrix!(A, bedfin)
        end
    end

    return A::Matrix{T}
end

function BEDintomatrix!{T<:Real}(A::AbstractMatrix{T}, bedvector::Vector{UInt8})
    n, p = size(A)
    byteheight = ceil(Int, n/4)
    quarterstop = n - 4*(byteheight - 1)

    @inbounds for col in 1:p
        unsafe_copybytestosnps!(A, bedvector, byteheight*(col - 1) + 3 + 1, 1, byteheight*col + 3, quarterstop, n*(col - 1) + 1)
    end

    return A
end

"""
    BEDintomatrix!{T<:Real}(A::AbstractMatrix{T}, bedstream::IO)

Fills `A` with type `T` representation of `bedstream` that corresponds
to .bed file format. `A` must have correct dimensions, as determined
via `BEDsize`, for example.

"""
function BEDintomatrix!{T<:Real}(A::AbstractMatrix{T}, bedstream::IO)
    n, p = size(A)
    bytestop = ceil(Int, n/4)
    quarterstop = n - 4*(bytestop - 1)

    bytecol = Vector{UInt8}(bytestop)
    @inbounds for col in 1:p
        read!(bedstream, bytecol)
        unsafe_copybytestosnps!(A, bytecol, 1, 1, bytestop, quarterstop, n*(col - 1) + 1)
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

"""
    immutable BEDMatrix{T<:Real, S<:AbstractMatrix} <: DenseArray{T, 2}

`Matrix{T}`-like representation of data in a literal BED matrix of
type `S`. Intended usage is for the BED matrix `X::S` to be memory
mapping of .bed file. That is `X` should be in compact, 4 SNPs per
byte, format with padded columns.

## Fields:
    n::Int
    p::Int
    X::S
    path::String
    colnames::Vector{String}
    rownames::Vector{String}

    _byteheight::Int  # number of bytes in each column
    _lastrowSNPheight::Int  # number of SNPs in last byte of each column ∈ (1, 2, 3, 4)

"""
immutable BEDMatrix{T<:Real, S<:AbstractMatrix} <: DenseArray{T, 2}
    n::Int
    p::Int
    X::S
    path::String
    colnames::Vector{String}
    rownames::Vector{String}

    _byteheight::Int  # number of bytes in each column
    _lastrowSNPheight::Int  # number of SNPs in last byte of each column ∈ (1, 2, 3, 4)

    function BEDMatrix(n::Integer, p::Integer, X::AbstractMatrix{UInt8}, path::AbstractString, colnames::AbstractVector, rownames::AbstractVector)
        byten = ceil(Int, n/4)
        lastrowheight = n - 4*(byten - 1)

        size(X) == (byten, p) || throw(DimensionMismatch("Matrix dimensions $(size(X)) do not agree with supplied BED dimensions (n = $n, p = $p)"))
        length(colnames) == p || throw(DimensionMismatch("colnames has incorrect length"))
        length(rownames) == n || throw(DimensionMismatch("rownames has incorrect length"))

        return new(n, p, X, path, colnames, rownames, byten, lastrowheight)
    end
end

"""
    BEDMatrix(bedfilename::AbstractString, U::DataType=UInt8)

Create a `BEDMatrix` of type `U` using memory mapping of BED file,
`bedfilename`.

"""
function BEDMatrix(bedfilename::AbstractString, U::DataType=UInt8)
    isfile(bedfilename) || error("Cannot find file \"$bedfilename\"")
    endswith(bedfilename, ".bed") || error("File does not have .bed extension")

    filebase = splitext(bedfilename)[1]
    rownames = readrownames(filebase*".fam")
    colnames = readcolnames(filebase*".bim")
    n, p = length(rownames), length(colnames)

    X = open(filebase*".bed", "r") do bedfile
        checkmagic(bedfile)
        if BEDmode(bedfile) != :SNPmajor
            error("SNPminor mode not supported")
        end

        Mmap.mmap(bedfile, Matrix{UInt8}, (ceil(Int, n/4), p))
    end

    return BEDMatrix{U, typeof(X)}(n, p, X, abspath(bedfilename), colnames, rownames)
end

function parsefamline(line)
    fields = split(line)

    (fields[1], fields[2], parse(Int, fields[3]), parse(Int, fields[4]), parse(Int, fields[5]), parse(Int, fields[6]))
end

famrowname(line) = join(split(line)[1:2], '_')

readrownames(famfile::AbstractString) = map(famrowname, readlines(famfile))

function bimcolname(line)
    fields = split(line)
    string(fields[2], '_', fields[5])
end

readcolnames(bimfile::AbstractString) = map(bimcolname, readlines(bimfile))

Base.size(B::BEDMatrix) = (B.n, B.p)
function Base.size(B::BEDMatrix, k::Integer)
    k > 0 || Base.throw_boundserror(size(B), k)
    # This is the same behavior as existing Base.size methods
    return k == 1 ? B.n : ifelse(k == 2, B.p, 1)
end

Base.linearindexing{T<:BEDMatrix}(::Type{T}) = Base.LinearSlow()

Base.sizeof(B::BEDMatrix) = sizeof(B.X)

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

getrow(B::BEDMatrix, rowname::AbstractString) = findfirst(B.rownames, rowname)
getcol(B::BEDMatrix, colname::AbstractString) = findfirst(B.colnames, colname)


#################### Indexing ####################

Base.getindex{T<:AbstractString}(B::BEDMatrix, rownames::AbstractVector{T}, col) = B[map(name -> getrow(B, name), rownames), col]
Base.getindex{T<:AbstractString}(B::BEDMatrix, row, colnames::AbstractVector{T}) = B[row, map(name -> getcol(B, name), colnames)]
Base.getindex{T<:AbstractString, S<:AbstractString}(B::BEDMatrix, rownames::AbstractVector{S}, colnames::AbstractVector{T}) = B[map(name -> getrow(B, name), rownames), map(name -> getcol(B, name), colnames)]


Base.getindex(B::BEDMatrix, rowname::AbstractString, col::Integer) = B[getrow(B, rowname), col]

Base.getindex(B::BEDMatrix, row::Integer, colname::AbstractString) = B[row, getcol(B, colname)]

Base.getindex(B::BEDMatrix, rowname::AbstractString, colname::AbstractString) = B[getrow(B, rowname), getcol(B, colname)]

# This is is the only getindex method we _need_, the other methods are
# provided for better performance, or convenience.
function Base.getindex(B::BEDMatrix, row::Integer, col::Integer)
    @boundscheck checkbounds(B, row, col)

    unsafe_getindex(B, row, col)
end

function unsafe_getindex{T, S}(B::BEDMatrix{T, S}, row::Integer, col::Integer)
    # Curse you julia for your foolish 1-indexing!
    byterow = (row - 1) >> 2 + 1
    snpind = (row - 1) & 3 + 1

    # @inbounds snp = convert(T, quarter(B.X[byterow, col], snpind))
    @inbounds snp = convert(T, breakbyte(B.X[byterow, col])[snpind])
    return snp
end

function Base.getindex(B::BEDMatrix, ::Colon, col::Integer)
    @boundscheck checkbounds(B, :, col)

    unsafe_getcol(B, col)
end

function Base.getindex{T, S, K<:Integer}(B::BEDMatrix{T, S}, ::Colon, cols::AbstractVector{K})
    @boundscheck checkbounds(B, :, cols)

    matrix = Matrix{T}(B.n, length(cols))

    for (cidx, col) in enumerate(cols)
        unsafe_getcol!(matrix, B, col, (cidx - 1)*B.n + 1)
    end
    return matrix
end

function unsafe_getcol{T, S}(B::BEDMatrix{T, S}, col::Integer)
    column = Vector{T}(B.n)

    unsafe_getcol!(column, B, col)
    return column
end

function unsafe_getcol!{T, S}(column::AbstractArray{T}, B::BEDMatrix{T, S}, col::Integer, deststart=1)
    unsafe_copybytestosnps!(column, B.X, B._byteheight*(col - 1) + 1, 1, B._byteheight*col, B._lastrowSNPheight, deststart)
    return column
end

function Base.getindex(B::BEDMatrix, rrange::UnitRange, col::Integer)
    @boundscheck checkbounds(B, rrange, col)

    unsafe_getrowrange(B, rrange, col)
end

function Base.getindex{T, S, K <: Integer}(B::BEDMatrix{T, S}, rrange::UnitRange, cols::AbstractVector{K})
    @boundscheck checkbounds(B, rrange, cols)

    matrix = Matrix{T}(length(rrange), length(cols))
    for (cidx, col) in enumerate(cols)
        unsafe_getrowrange!(matrix, B, rrange, col, (cidx - 1)*length(rrange) + 1)
    end
    return matrix
end

function unsafe_getrowrange{T, S}(B::BEDMatrix{T, S}, rrange::UnitRange, col::Integer)
    vector = Vector{T}(length(rrange))
    unsafe_getrowrange!(vector, B, rrange, col)
    return vector
end

function unsafe_getrowrange!{T, S}(vector::AbstractArray{T}, B::BEDMatrix{T, S}, rrange::UnitRange, col::Integer, deststart=1)
    bytestart = B._byteheight*(col - 1) + (rrange.start - 1) >> 2 + 1
    quarterstart = (rrange.start - 1) & 3 + 1
    bytestop = B._byteheight*(col - 1) + (rrange.stop - 1) >> 2 + 1
    quarterstop = (rrange.stop - 1) & 3 + 1

    unsafe_copybytestosnps!(vector, B.X, bytestart, quarterstart, bytestop, quarterstop, deststart)
    return vector
end


######################### Other stuff #########################

# Some ideas from bitarray.jl

## iterator ?
# Base.start
# Base.next
# Base.done

## map ?

## filter

## hcat, vcat -- LinkedMatrix?

## write to other formats?

# Base.show() to display NA_byte as NA

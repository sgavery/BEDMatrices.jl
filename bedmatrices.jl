

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
@inline function rawformat(snp::Integer, quartermap=quarterstohuman)
    @inbounds return quartermap[snp + 1]
end

"""
    breakbyte(byte::UInt8, bytemap=bytetoquarters)

Return length-4 `Vector{UInt8}` of the RAW-formatted SNP quarters, as
determined by `bytemap`, in `byte`.

"""
function breakbyte(byte::UInt8, bytemap=bytetoquarters)
    @inbounds return bytemap[byte + 1]
end

# Approach from Quantgen/BEDMatrix: shift by 2(n - 1) and then mask
# result.
#
# Performance-wise, there does not seem to be a clear advantage of one
# approach over the other.
"""
    quarter(byte::UInt8, n::Integer, quartermap=quarterstohuman)

Returns the RAW-formatted `n`th quarter of `byte`, equivalent to
`breakbyte(byte)[n]`.

"""
function quarter(byte::UInt8, n::Integer, quartermap=quarterstohuman)
    @inbounds return quartermap[((byte >> 2(n - 1)) & 0b00000011) + 1]
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

"""
checkmagic(bytes::Vector{UInt8}) = (bytes[1:2] == plinkmagic || error("Bad magic: not a plink bed file"))

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
    unsafe_breakbyte!(vector::AbstractArray, byte::UInt8, bytemap=bytetoquarters, vecstart=1, quarterstart=1, quarterstop=4)

Fills `vector` starting at `vecstart` with the `quarterstart`th snp
through the `quarterstop`th snp in `byte`, as determined by `bytemap`.

## (Unchecked) Constraints:
* 1 <= `quarterstart` <= `quarterstop` <= 4
* 1 <= `vecstart` <= `end - (quarterstop - quarterstart)`

"""
function unsafe_breakbyte!(vector::AbstractArray, byte::UInt8, bytemap=bytetoquarters, vecstart=1, quarterstart=1, quarterstop=4)
    @inbounds copy!(vector, vecstart,
                    breakbyte(byte, bytemap), quarterstart,
                    quarterstop - quarterstart + 1)
    vector
end


"""
    unsafe_copybytestosnps!(snparray::AbstractArray, bytearray::AbstractArray{UInt8},
                            bytestart::Integer, quarterstart::Integer,
                            bytestop::Integer, quarterstop::Integer,
                            deststart::Integer=1, bytemap=bytetoquarters)

Fills `snparray[deststart:(deststart + num_snps)]` with snps, where

    num_snps = 4*(bytestop - bytestart) + (quarterstop - quarterstart) + 1

The SNPs start with the `quarterstart`th snp in the `bytestart`th
byte, and end with `quarterstop` and `bytestop`.  The SNP
representation is determined by the `bytemap`. Currently only supports
unit stride for both `snparray` and `bytearray`, since at present
there is no compelling use-case involving non-unit strides.

Utility function used by more front facing functions.

"""
function unsafe_copybytestosnps!(snparray::AbstractArray, bytearray::AbstractArray{UInt8},
                                 bytestart::Integer, quarterstart::Integer,
                                 bytestop::Integer, quarterstop::Integer,
                                 deststart::Integer=1, bytemap=bytetoquarters)
    # First byte
    if quarterstart != 1
        stop = ifelse(bytestart == bytestop, quarterstop, 4)
        @inbounds unsafe_breakbyte!(snparray, bytearray[bytestart], bytemap, deststart, quarterstart, stop)
        deststart += stop - quarterstart + 1
        bytestart += 1
    end

    if bytestart <= bytestop
        # Last byte
        if quarterstop != 4
            @inbounds unsafe_breakbyte!(snparray, bytearray[bytestop], bytemap, deststart + 4*(bytestop - bytestart), 1, quarterstop)
            bytestop -= 1
        end

        # Main course
        @simd for bytej in bytestart:bytestop
            @inbounds unsafe_breakbyte!(snparray, bytearray[bytej], bytemap, deststart + 4*(bytej - bytestart))
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
via `BEDMatrices.BEDsize`, for example.

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

immutable BEDMatrix{T, S<:AbstractMatrix} <: DenseArray{T, 2}
    n::Int
    p::Int
    X::S
    navalue::T

    path::String
    colnames::Vector{String}
    rownames::Vector{String}

    _byteheight::Int  # number of bytes in each column
    _lastrowSNPheight::Int  # number of SNPs in last byte of each column ∈ (1, 2, 3, 4)

    _bytemap::Vector{Vector{T}}  # quarters for 0x00:0xff

    function BEDMatrix(n::Integer, p::Integer, X::AbstractMatrix{UInt8}, navalue, path::AbstractString, colnames::AbstractVector, rownames::AbstractVector)
        byteheight = ceil(Int, n/4)
        lastrowheight = n - 4*(byteheight - 1)

        size(X) == (byteheight, p) || throw(DimensionMismatch("Matrix dimensions $(size(X)) do not agree with supplied BED dimensions (n = $n, p = $p)"))
        length(colnames) == p || throw(DimensionMismatch("colnames has incorrect length"))
        length(rownames) == n || throw(DimensionMismatch("rownames has incorrect length"))

        quartermap = (convert(T, 0b10), navalue, convert(T, 0b01), convert(T, 0b00))
        bytemap = [[rawformat(snp1, quartermap), rawformat(snp2, quartermap),
                    rawformat(snp3, quartermap), rawformat(snp4, quartermap)] for
                        snp4 in 0b00:0b11 for
                        snp3 in 0b00:0b11 for
                        snp2 in 0b00:0b11 for
                        snp1 in 0b00:0b11]

        return new(n, p, X, navalue, path, colnames, rownames,
                   byteheight, lastrowheight, bytemap)
    end
end

"""
    BEDMatrix(bedfilename::AbstractString;
              datatype::DataType=UInt8, nsamples::Integer=0, nSNPs::Integer=0, navalue=NA_byte
              famfile::AbstractString="", bimfile::AbstractString="")

Create a `BEDMatrix` of type `datatype` using memory mapping of BED
file, `bedfilename`. Use `navalue` for missing values. `famfile` and
`bimfile` are only required if the .fam and .bim file are not in same
location as or do not share the same base name as the
`bedfilename`. `nsamples` and `nSNPs` are only used if .bim and .fam
file cannot be found and are not provided; in which case the
rownames/colnames will be generic.

`convert(datatype, x)` must work for `x` in `[0b00, 0b01, 0b10, navalue]`.

## Examples
```julia
julia> bed = BEDMatrix("test/data/example.bed");

julia> bed[1:5, 1:5]
5×5 Array{UInt8,2}:
 0x00  0x01  0x01  0x01  0x00
 0x01  0x01  0x01  0x01  0x03
 0x01  0x00  0x00  0x02  0x00
 0x02  0x00  0x00  0x00  0x01
 0x00  0x01  0x00  0x00  0x00

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
                   datatype::DataType=UInt8, nsamples::Integer=0, nSNPs::Integer=0, navalue=NA_byte,
                   famfile::AbstractString="", bimfile::AbstractString="")
    if !isfile(bedfilename)
        isfile(bedfilename*".bed") || error("Cannot find file \"$bedfilename\"")
        bedfilename = bedfilename*".bed"
    end
    filebase = splitext(bedfilename)[1]

    famfile = famfile == "" ? filebase*".fam" : famfile
    if isfile(famfile)
        rownames = readrownames(famfile)
    elseif nsamples > 0
        rownames = [string("sample_", j) for j in 1:nsamples]
    else
        error("Cannot find FAM file \"$famfile\" and nsamples not provided")
    end

    bimfile = bimfile == "" ? filebase*".bim" : bimfile
    if isfile(bimfile)
        colnames = readcolnames(bimfile)
    elseif nSNPs > 0
        colnames = [string("SNP_", j) for j in 1:nSNPs]
    else
        error("Cannot find BIM file \"$bimfile\" and nSNPs not provided")
    end

    n, p = length(rownames), length(colnames)

    X = open(filebase*".bed", "r") do bedfile
        checkmagic(bedfile)
        BEDmode(bedfile) == :SNPmajor || error("Old-style SNP-minor mode bed files not supported. Use plink to convert to SNP-major format")

        Mmap.mmap(bedfile, Matrix{UInt8}, (ceil(Int, n/4), p))
    end

    return BEDMatrix{datatype, typeof(X)}(n, p, X, convert(datatype, navalue),
                                          abspath(bedfilename), colnames, rownames)
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

# not sure if this is the right definition; see
# https://github.com/JuliaLang/julia/issues/16614.
# This is the size of `B[:, :]`.
Base.sizeof{T, S}(B::BEDMatrix{T, S}) = B.n*B.p*sizeof(T)

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

function Base.getindex(B::BEDMatrix, rrange::AbstractUnitRange, col::Integer)
    @boundscheck checkbounds(B, rrange, col)

    unsafe_getrowrange(B, rrange, col)
end

function Base.getindex{T, S, K <: Integer}(B::BEDMatrix{T, S}, rrange::AbstractUnitRange, cols::AbstractVector{K})
    @boundscheck checkbounds(B, rrange, cols)

    matrix = Matrix{T}(length(rrange), length(cols))
    for (cidx, col) in enumerate(cols)
        unsafe_getrowrange!(matrix, B, rrange, col, (cidx - 1)*length(rrange) + 1)
    end
    return matrix
end

function unsafe_getindex{T, S}(B::BEDMatrix{T, S}, row::Integer, col::Integer)
    byterow, snpind = rowtobytequarter(row)

    # @inbounds snp = convert(T, quarter(B.X[byterow, col], snpind))
    @inbounds return breakbyte(B.X[byterow, col], B._bytemap)[snpind]
end

function unsafe_getcol{T, S}(B::BEDMatrix{T, S}, col::Integer)
    column = Vector{T}(B.n)

    unsafe_getcol!(column, B, col)
    return column
end

function unsafe_getcol!{T, S}(column::AbstractArray{T}, B::BEDMatrix{T, S}, col::Integer, deststart=1)
    unsafe_copybytestosnps!(column, B.X, B._byteheight*(col - 1) + 1, 1, B._byteheight*col, B._lastrowSNPheight, deststart, B._bytemap)
    return column
end

function unsafe_getrowrange{T, S}(B::BEDMatrix{T, S}, rrange::AbstractUnitRange, col::Integer)
    vector = Vector{T}(length(rrange))
    unsafe_getrowrange!(vector, B, rrange, col)
    return vector
end

function unsafe_getrowrange!{T, S}(vector::AbstractArray{T}, B::BEDMatrix{T, S}, rrange::AbstractUnitRange, col::Integer, deststart=1)
    bytestart, quarterstart = rowtobytequarter(first(rrange))
    bytestart += B._byteheight*(col - 1)

    bytestop, quarterstop = rowtobytequarter(last(rrange))
    bytestop += B._byteheight*(col - 1)

    unsafe_copybytestosnps!(vector, B.X, bytestart, quarterstart, bytestop, quarterstop, deststart, B._bytemap)
    return vector
end

function rowtobytequarter(row::Integer)
    # Curse you julia for your foolish 1-indexing!
    return ((row - 1) >> 2 + 1, (row - 1) & 3 + 1)
end

"""
    getquarterblock(B::BEDMatrix, row::Integer, col::Integer)

Returns `(fourquarters, snpind)`, where `fourquarters` is the block of
quarters that `B[row, col]` belongs to: `B[row, col] =
fourquarters[snpind]`.

"""
function getquarterblock(B::BEDMatrix, row::Integer, col::Integer)
    byterow, snpind = rowtobytequarter(row)
    @inbounds return (breakbyte(B.X[byterow, col], B._bytemap), byterow, snpind)
end

"""
    tocontiguous(index::AbstractVector{Bool})

Takes a logical index (vector of `Bool`s for each index) and converts
it into a `Vector` of `UnitRange`s. This is especially useful for
multicolumn or repeated uses of the same `index` to exploit the
performance benefits of byte-wise calculations.

"""
function tocontiguous(index::AbstractVector{Bool})
    n = length(index)

    ranges = Vector{UnitRange{Int}}(n)
    start = findfirst(index)

    # If there are no indices return []
    if start == 0
        return resize!(ranges, 0)
    end

    stop = findlast(index)

    oldidx = idx = start - 1
    rangeidx = 0
    while idx < stop
        oldidx = idx
        idx = findnext(!, index, idx + 1)
        idx == 0 && break

        if oldidx + 1 < idx
            rangeidx += 1
            ranges[rangeidx] = (oldidx + 1):(idx - 1)
        end
    end

    if idx == 0 && oldidx < stop
        rangeidx += 1
        ranges[rangeidx] = (oldidx + 1):stop
    end

    resize!(ranges, rangeidx)
end

"""
    tocontiguous{T<:Integer}(index::AbstractVector{T})

Takes a `Vector` of indices and converts it into a `Vector` of
`UnitRange`s.

"""
function tocontiguous{T<:Integer}(index::AbstractVector{T})
    n = length(index)
    ranges = Vector{UnitRange{Int}}(n)

    oldidx = idx = 1
    rangeidx = 0
    while idx <= n
        oldidx = idx
        idx += 1
        while idx <= n && (index[idx] == index[oldidx] + (idx - oldidx))
            idx += 1
        end
        rangeidx += 1
        ranges[rangeidx] = index[oldidx]:index[idx - 1]
    end

    resize!(ranges, rangeidx)
end

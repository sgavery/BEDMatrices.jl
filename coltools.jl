

#################### BED byte to RAW math ####################

"""
    maskbyte(bedbyte::UInt8, num_quarters::Integer, qoffset::Integer=0)

Reads `byte` in BED format and returns `byte`, after masking quarters
such that `(qoffset + 1):(num_quarters + qoffset)` are preserved and
other quarters are "zeroed" (with respect to RAW format). Assumes RAW
format!!

"""
function maskbyte(bedbyte::UInt8, num_quarters::Integer, qoffset::Integer=0)
    @inbounds return (bedbyte | partialbytemasks[num_quarters, qoffset + 1])
end

"""
    hasNAs(bedbyte::UInt8)

Returns `true` if `bedbyte` has any missing values, and `false`
otherwise.

"""
function hasNAs(bedbyte::UInt8)
    @inbounds return hasNAmap[bedbyte + 1]
end

"""
    zeroNAs(bedbyte::UInt8)

Reads `bedbyte` in BED format and returns `bedbyte` with any missing
values set to zero. Assumes RAW format!!

"""
function zeroNAs(bedbyte::UInt8)
    @inbounds return natozeromap[bedbyte + 1]
end

"""
    countNAs(bedbyte::UInt8)

Returns the number of missing values in `bedbyte`.

"""
function countNAs(bedbyte::UInt8)
    @inbounds return nacountmap[bedbyte + 1]
end


function bytedot(b1::UInt8, b2::UInt8)
    @inbounds return bytebytemulttable[b1 + 1, b2 + 1]
end

"""
    bytedot(byteL::UInt8, byteR::UInt8, num_quarters::Integer=4, qoffset::Integer=0)

Return RAW format dot product of two bytes in BED format over quarters
in `(qoffset + 1):(qoffset + num_quarters)`. Missing values are
treated as zero.

"""
function bytedot(byteL::UInt8, byteR::UInt8, num_quarters::Integer, qoffset::Integer=0)
    @boundscheck qoffset + num_quarters <= 4

    bytedot(maskbyte(byteL, num_quarters, qoffset), byteR)
end

function bytedot{T}(byte::UInt8, v::AbstractArray{T}, voffset::Integer=0)
    @boundscheck checkbounds(voffset + 4)

    dotsum = zero(T)
    quarters = breakbyte(zeroNAs(byte))

    for j in 1:4
        @inbounds q = quarters[j]
        @inbounds dotsum += q*v[j + voffset]
    end

    dotsum
end

"""
    bytedot(byte::UInt8, v::AbstractArray, voffset::Integer=0, num_quarters::Integer=4, qoffset::Integer=0)

Return RAW format dot product of `byte[(qoffset + 1):(num_quarters +
qoffset)]` in BED format with `v[(voffset + 1):(voffset +
num_quarters)]`. Missing values (in `byte`) are treated as zero.

"""
function bytedot(byte::UInt8, v::AbstractArray, voffset::Integer, num_quarters::Integer, qoffset::Integer=0)
    @boundscheck num_quarters + qoffset <= 4 && checkbounds(v, voffset + 4)

    bytedot(maskbyte(byte, num_quarters, qoffset), v, voffset)
end

bytesum(byte::UInt8) = bytedot(byte, onebyte)

function bytesum(byte::UInt8, num_quarters::Integer, qoffset::Integer=0)
    @boundscheck num_quarters + qoffset <= 4

    @inbounds return bytedot(byte, onebyte, num_quarters, qoffset)
end


#################### Column tools ####################


funcvector(B::BEDMatrix, func::Function) = map(quarters -> sum(map(func, quarters)), B._bytemap)
funcvector(B::BEDMatrix, func::Function, red_op::Function) = map(quarters -> reduce(red_op, map(func, quarters)), B._bytemap)

function column_mapreduce(B::BEDMatrix, col::Integer, rows::AbstractUnitRange, func::Function, red_op, red_id)
end

function column_mapreduce(B::BEDMatrix, col::Integer, func::Function, red_op, red_id)
end

# fallback method
function column_mapreduce(B::BEDMatrix, col::Integer, rows, func, red_op, red_id)
    @boundscheck checkbounds(B, rows, col)
    total = red_id

    for r in rows
       @inbounds total = red_op(total, func(B[r, col]))
    end

    total
end

function createormask(inds::AbstractArray)
    mask = 0b11_11_11_11

    for qidx in inds
        mask -= 0b11 << 2*(qidx - 1)
    end

    mask
end

function column_mapsum{T}(B::BEDMatrix, col::Integer, funcvec::Vector{T}, lastmask::UInt8)
    total = zero(promote_type(T, Int))

    X = B.X
    nbytes = B._byteheight

    @simd for x in 1:(nbytes - 1)
        @inbounds total += funcvec[X[x, col] + 1]
    end
    @inbounds total += funcvec[(X[nbytes, col] | lastmask) + 1]

    total
end

function column_mapreduce{T}(B::BEDMatrix, col::Integer,
                             funcvec::Vector{T}, lastfuncvec::Vector{T},
                             red_op::Function, id_el::T)
    total = id_el

    X = B.X
    nbytes = B._byteheight

    for x in 1:(nbytes - 1)
        @inbounds total = red_op(total, funcvec[X[x, col] + 1])
    end
    @inbounds total = red_op(total, lastfuncvec[X[x, col] + 1])

    total
end


"""
    columnNAcount(B::BEDMatrix{T, S}, col::Integer)

Returns number of NAs in `B[:, col]`. Use `columnmoments` if you need
 other moments.

"""
function columnNAcount(B::BEDMatrix, col::Integer)
    columnmoments(B, Val{0})
end

# function columnmoments{T, S, N}(B::BEDMatrix{T, S}, col::Integer, ::Type{Val{N}})
# end

function columnmoments(B::BEDMatrix, col::Integer, ::Type{Val{0}})
    nacount = 0

    X = B.X
    nbytes = B._byteheight
    @inbounds begin
        for x in 1:(nbytes - 1)
            nacount += countNAs(X[x, col])
        end
        nacount += countNAs(maskbyte(X[nbytes, col], B._lastrowSNPheight))
    end

    nacount
end

function columnmoments{T, S}(B::BEDMatrix{T, S}, col::Integer, ::Type{Val{1}})
    nacount = 0
    total = zero(promote_type(T, Int))  # Avoid overflow when T === UInt8

    X = B.X
    nbytes = B._byteheight

    @inbounds begin
        for x in 1:(nbytes - 1)
            b = X[x, col]
            nacount += countNAs(b)
            total += bytesum(b)
        end

        b = maskbyte(X[nbytes, col], B._lastrowSNPheight)
        nacount += countNAs(b)
        total += bytesum(b)
    end

    (nacount, total)
end

function columnmoments{T, S}(B::BEDMatrix{T, S}, col::Integer, ::Type{Val{2}})
    nacount = 0
    total = zero(promote_type(T, Int))  # Avoid overflow when T === UInt8
    squares = total

    X = B.X
    nbytes = B._byteheight

    @inbounds begin
        for x in 1:(nbytes - 1)
            b = X[x, col]
            nacount += countNAs(b)
            total += bytesum(b)
            squares += bytedot(b, b)
        end

        b = maskbyte(X[nbytes, col], B._lastrowSNPheight)
        nacount += countNAs(b)
        total += bytesum(b)
        squares += bytedot(b, b)
    end

    (nacount, total, squares)
end

"""
    columnmoments{T, S}(B::BEDMatrix{T, S}, col::Integer)

Returns number of NAs, sum of non-NAs, sum of non-NAs squared for
`B[:, col]`, in a computational efficient manner. These can be thought
of as the 0th, 1st, and 2nd moments of the column.

"""
function columnmoments{T, S}(B::BEDMatrix{T, S}, col::Integer)
    columnmoments(B, col, Val{2})
end

function columndot{T, S}(B::BEDMatrix{T, S}, col::Integer, v::AbstractVector)
    @boundscheck B.n == length(v)

    dotsum = zero(promote_type(T, eltype(v), Int))  # Avoid overflow when T == UInt8

    X = B.X
    nbytes = B._byteheight

    for x in 1:(nbytes - 1)
        @inbounds dotsum += bytedot(X[x, col], v, 4*(x-1))
    end
    @inbounds dotsum += bytedot(X[nbytes, col], v, 4*(nbytes-1), B._lastrowSNPheight)

    dotsum
end

# https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance may be of interest

function colmean(B::BEDMatrix, col::Integer)
    nacount, total = columnmoments(B, col, Val{1})

    total/(B.n - nacount)
end

function colmeanstd(B::BEDMatrix, col::Integer)
    total, totsq, nacount = columnmoments(B, col, Val{2})

    μ = total/(B.n - nacount)
    var = (totsq - 2μ*total + (B.n - nacount)*abs2(μ))/(B.n - nacount - 1)

    return μ, var
end

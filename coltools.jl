

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

function byteabs2(b::UInt8)
    @inbounds return bytebytemulttable[b + 1, b + 1]
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

"""
    bytesum(byte::UInt8, num_quarters::Integer=4, qoffset::Integer=0)

Returns sum of non-missing values in RAW format in `byte[(qoffset +
1):(num_quarters + qoffset)]`.

"""
function bytesum(byte::UInt8, num_quarters::Integer, qoffset::Integer=0)
    @boundscheck num_quarters + qoffset <= 4

    @inbounds return bytedot(byte, onebyte, num_quarters, qoffset)
end

"""
    funcvector(func::Function, B::BEDMatrix)

Convert `func` with return type `T` into a 256-length `Vector{T}`
corresponding to `sum(func, quarters)` for all bytes. This is used for
efficiently `column_sum`-ing on large or multiple columns.

"""
funcvector(func::Function, B::BEDMatrix) = map(quarters -> sum(func, quarters), B._bytemap)


"""
    funcvector(func::Function, red_op::Function, B::BEDMatrix)

Convert `func` with return type `T` into a 256-length `Vector{T}`
corresponding to `mapreduce(func, red_op, quarters)` for all
bytes. This is used for efficiently `mapreduce`-ing on large or
multiple columns.

"""
funcvector(func::Function, red_op::Function, B::BEDMatrix) = map(quarters -> mapreduce(func, red_op, quarters), B._bytemap)

"""
    funcquartervector{T}(func::Function, navalue::T)

Convert `func` with return type `S` into a length-4 `Vector{S}`
corresponding to how `func` acts on the RAW/`navalue` representation
of the BED quarters `[0b00, 0b01, 0b10, 0b11]`. This fully specifies
the behavior of `func` on a `BEDMatrix{T, Q}`.

"""
funcquartervector{T}(func::Function, navalue::T) = map(func, [convert(T, 0b00), convert(T, 0b01), convert(T, 0b10), navalue])


#################### SubArray Interface ####################

typealias Column{T, K, B} SubArray{T, 1, K, Tuple{Colon, Int}, B}
typealias BEDColumn{T, K <: BEDMatrix} Column{T, K, false}
typealias ColumnUnitRange{T, K, R <: Union{AbstractUnitRange, Colon}, B} SubArray{T, 1, K, Tuple{R, Int}, B}
typealias BEDColumnUnitRange{T, K <: BEDMatrix, R <: Union{AbstractUnitRange, Colon}} ColumnUnitRange{T, K, R, false}

"""
    hasNAs(B::BEDMatrix)

Returns `true` if there are missing values in `B`.

"""
function hasNAs(B::BEDMatrix)
    for col in indices(B, 2)
        if column_hasNAs(B, col)
            return true
        end
    end

    return false
end

hasNAs(BC::BEDColumn) = column_hasNAs(parent(BC), parentindexes(BC)[2])

"""
    hasNAs(BC::BEDColumn)

Returns `true` if there are missing values in `BC`.

"""
function hasNAs(BC::BEDColumnUnitRange)
    rrange, col = parentindexes(BC)
    column_hasNAs(parent(BC), col, rrange)
end

function countNAs(B::BEDMatrix)
    total = 0
    for col in indices(B, 2)
        total += column_countNAs(B, col)
    end
    total
end

countNAs(BC::BEDColumn) = column_countNAs(parent(BC), parentindexes(BC)[2])

Base.sum(BC::BEDColumn) = column_sum(parent(BC), parentindexes(BC)[2])
function Base.sum(BC::BEDColumnUnitRange)
    rrange, col = parentindexes(BC)
    column_sum(parent(BC), col, rrange)
end
# Base.sum(func, ...)

Base.sumabs(BC::BEDColumnUnitRange) = sum(BC)

Base.sumabs2(BC::BEDColumn) = column_sumabs2(parent(BC), parentindexes(BC)[2])
function Base.sumabs2(BC::BEDColumnUnitRange)
    rrange, col = parentindexes(BC)
    column_sumabs2(parent(BC), col, rrange)
end

# Base.dot
# Base.mapreduce
# Base.norm


#################### Column tools ####################

"""
    column_sum{T, S}(B::BEDMatrix{T, S}, col::Integer)

Returns the sum of the non-missing entries in `B[:, col]`.

"""
function column_sum{T, S}(B::BEDMatrix{T, S}, col::Integer)
    total = zero(promote_type(T, Int))  # Avoid overflow when T === UInt8

    X = B.X
    nbytes = B._byteheight

    @inbounds begin
        for x in 1:(nbytes - 1)
            total += bytesum(X[x, col])
        end

        total += bytesum(X[nbytes, col], B._lastrowSNPheight)
    end

    total
end

"""
    column_sum{T, S}(B::BEDMatrix{T, S}, col::Integer, rows)

Returns the sum of non-missing entries in `B[rows, col]`, optimized
for large `UnitRange`s of rows.

"""
function column_sum{T, S}(B::BEDMatrix{T, S}, col::Integer, rrange::AbstractUnitRange)
    total = zero(promote_type(T, Int))  # Avoid overflow when T === UInt8

    X = B.X
    bytestart, quarterstart = rowtobytequarter(first(rrange))
    bytestop, quarterstop = rowtobytequarter(last(rrange))

    @inbounds begin
        if quarterstart > 1
            stop = bytestart == bytestop ? quarterstop : 4
            total += bytesum(maskbyte(X[bytestart, col], stop - quarterstart + 1, quarterstart - 1))
            bytestart += 1
        end

        if bytestop >= bytestart
            total += bytesum(maskbyte(X[bytestop, col], quarterstop))
            bytestop -= 1

            for x in bytestart:bytestop
                total += bytesum(X[x, col])
            end
        end
    end

    total
end

function column_sumabs2{T, S}(B::BEDMatrix{T, S}, col::Integer)
    squares = zero(promote_type(T, Int))  # Avoid overflow when T === UInt8

    X = B.X
    nbytes = B._byteheight

    @inbounds begin
        for x in 1:(nbytes - 1)
            squares += byteabs2(X[x, col])
        end

        squares += byteabs2(maskbyte(X[nbytes, col], B._lastrowSNPheight))
    end

    squares
end

function column_sumabs2{T, S}(B::BEDMatrix{T, S}, col::Integer, rrange::AbstractUnitRange)
    squares = zero(promote_type(T, Int))  # Avoid overflow when T === UInt8

    X = B.X
    bytestart, quarterstart = rowtobytequarter(first(rrange))
    bytestop, quarterstop = rowtobytequarter(last(rrange))

    @inbounds begin
        if quarterstart > 1
            stop = bytestart == bytestop ? quarterstop : 4
            squares += byteabs2(maskbyte(X[bytestart, col], stop - quarterstart + 1, quarterstart - 1))
            bytestart += 1
        end

        if bytestop >= bytestart
            squares += byteabs2(maskbyte(X[bytestop, col], quarterstop))
            bytestop -= 1

            for x in bytestart:bytestop
                squares += byteabs2(X[x, col])
            end
        end
    end

    squares
end

function qvec_column_sum(func::Function, B::BEDMatrix, col::Integer)
    func_qvec = funcquartervector(func, NArep(B))

    qvec_column_sum(func_qvec, B.X, col, 1, 1, B._byteheight, B._lastrowSNPheight)
end

function qvec_column_sum{T}(func_qvec::Vector{T}, X::Matrix{UInt8}, col::Integer, bytestart, quarterstart, bytestop, quarterstop)
    accumulator = zero(promote_type(T, Int))

    @inbounds begin
        if quarterstart > 1
            quarters = breakbyte(X[bytestart, col])
            stop = bytestart == bytestop ? quarterstop : 4
            for q in quarterstart:stop
                accumulator += func_qvec[quarters[q] + 1]
            end
            bytestart += 1
        end

        if bytestop >= bytestart
            quarters = breakbyte(X[bytestop, col])
            for q in 1:quarterstop
                accumulator += func_qvec[quarters[q] + 1]
            end
            bytestop -= 1

            for row in bytestart:bytestop
                quarters = breakbyte(X[row, col])
                for q in 1:4
                    accumulator += func_qvec[quarters[q] + 1]
                end
            end
        end
    end

    return accumulator
end

# This approach (vec_[...]) is optimized for multiple columns or very
# long columns, where the fixed cost of computing funcvec starts being
# offset by ≈ 4x speed up on main computation. It _never_ makes sense
# if total number of rows is less than ≈ 1024. Although, it starts
# making sense soon after. (This is because it precomputes func on all
# 256 bytes, each broken into 4 quarters.)
function vec_column_sum(func::Function, B::BEDMatrix, col::Integer)
    funcvec = funcvector(func, B)

    vec_bedcolumn_sum(funcvec, B.X, col, 1, 1, B._byteheight, B._lastrowSNPheight)
end

vec_column_sum(func, B::BEDMatrix, col::Integer, ::Colon) = vec_column_sum(func, B, col)

function vec_column_sum(func::Function, B::BEDMatrix, col::Integer, rrange::AbstractUnitRange)
    funcvec = funcvector(func, B)
    bytestart, quarterstart = rowtobytequarter(first(rrange))
    bytestop, quarterstop = rowtobytequarter(last(rrange))

    vec_bedcolumn_sum(funcvec, B.X, col, bytestart, quarterstart, bytestop, quarterstop)
end

function vec_bedcolumn_sum{T}(funcvec::Vector{T}, X::AbstractMatrix{UInt8}, col::Integer, bytestart, quarterstart, bytestop, quarterstop)
    accumulator = zero(T)

    if quarterstart > 1
        stop = bytestart == bytestop ? quarterstop : 4
        @inbounds accumulator += funcvec[maskbyte(X[bytestart, col], stop - quarterstart + 1, quarterstart - 1) + 1]
        bytestart += 1
    end

    if bytestart <= bytestop
        @inbounds accumulator += funcvec[maskbyte(X[bytestop, col], quarterstop) + 1]
        bytestop -= 1

        # @simd ?
        for row in bytestart:bytestop
            @inbounds accumulator += funcvec[X[row, col] + 1]
        end
    end

    accumulator
end

"""
    column_sum(func::Function, B::BEDMatrix, col::Integer)

Apply `func` to `B[:, col]`, and return the sum. It is the
responsibility of `func` to appropriately handle `NArep(B)` values.

"""
function column_sum(func::Function, B::BEDMatrix, col::Integer)
    if B.n > 2048
        return vec_column_sum(func, B, col)
    end
    bedcolumn_sum(func, B.X, B._bytemap, col, 1, 1, B._byteheight, B._lastrowSNPheight)
end

"""
    column_sum(func::Function, B::BEDMatrix, col::Integer, rows)

Apply `func` to `B[rows, col]` and return the sum. It is the
responsibility of `func` to appropriately handle `NArep(B)` values.

"""
function column_sum(func::Function, B::BEDMatrix, col::Integer, rrange::AbstractUnitRange)
    if length(rrange) > 2048
        return vec_column_sum(func, B, col, rrange)
    end
    bytestart, quarterstart = rowtobytequarter(first(rrange))
    bytestop, quarterstop = rowtobytequarter(last(rrange))

    bedcolumn_sum(func, B.X, B._bytemap, col, bytestart, quarterstart, bytestop, quarterstop)
end

# list indexing
# note: if one is just excluding a small fraction of row indices, it
# is better to convert to multiple `UnitRange`s instead and use vec_ version.
function column_sum{T <: Integer}(func::Function, B::BEDMatrix, col::Integer, rows::AbstractVector{T})
    @boundscheck checkbounds(B, rows, col)

    @inbounds accumulator = func(B[rows[1], col])
    for ridx in 2:length(rows)
        accumulator += func(B[rows[ridx], col])
    end

    accumulator
end

# logical indexing
# note: if one is just excluding a small fraction of row indices, it
# is better to convert to multiple `UnitRange`s instead and use vec_ version.
function column_sum{T <: Bool}(func::Function, B::BEDMatrix, col::Integer, rows::AbstractVector{T})
    @boundscheck checkbounds(B, rows, col)
    accumulator = zero(func(zero(eltype(B))))

    ridx = findnext(rows, 1)
    while ridx > 0
        @inbounds accumulator += func(B[rows[ridx], col])
        ridx = findnext(rows, ridx + 1)
    end

    accumulator
end

function bedcolumn_sum{T}(func::Function, X::AbstractMatrix{UInt8}, bytemap::Vector{Vector{T}}, col::Integer, bytestart, quarterstart, bytestop, quarterstop)
    accumulator = zero(func(breakbyte(X[bytestart, col], bytemap)[1]))

    if quarterstart > 1
        stop = bytestart == bytestop ? quarterstop : 4
        @inbounds quarters = breakbyte(X[bytestart, col], bytemap)
        for j in quarterstart:stop
            @inbounds accumulator += func(quarters[j])
        end
        bytestart += 1
    end

    if bytestart <= bytestop
        @inbounds quarters = breakbyte(X[bytestop, col], bytemap)
        for j in 1:quarterstop
            @inbounds accumulator += func(quarters[j])
        end
        bytestop -= 1

        for row in bytestart:bytestop
            @inbounds accumulator += sum(func, breakbyte(X[row, col], bytemap))
        end
    end

    accumulator
end

function column_mapreduce(func::Function, red_op, B::BEDMatrix, col::Integer, rows::AbstractUnitRange)
end

function column_mapreduce(func::Function, red_op, B::BEDMatrix, col::Integer)
    funcv = funcvector(B, func, red_op)

    B._byteheight
    X = B.X
    if B._lastrowSNPheight < 4
        bytestop -= 1
        accumulator = bedcolumn_mapreduce(X, col, 1, bytestop, funcv, red_op)

        @inbounds lastquarters = breakbyte(X[bytestop + 1, col], B._bytemap)
        for j in 1:B._lastrowSNPheight
            accumulator = red_op(accumulator, func(lastquarters[j]))
        end
    else
        accumulator = bedcolumn_mapreduce(X, col, 1, bytestop, funcv, red_op)
    end

    accumulator
end

function bedcolumn_mapreduce(funcvec::Vector, red_op::Function, X::Matrix{UInt8}, col::Integer, bytestart::Integer, bytestop::Integer)

    @inbounds accumulator = funcvec[X[bytestart, col] + 1]

    for row in (rowstart+1):bytestop
        @inbounds accumulator = red_op(accumulator, funcvec[X[row, col] + 1])
    end

    accumulator
end

# fallback method
function column_mapreduce(func, red_op, B::BEDMatrix, col::Integer, rows::AbstractArray)
    @boundscheck checkbounds(B, rows, col)

    @inbounds accumulator = func(B[rows[1], col])

    for ridx in 2:length(rows)
       @inbounds accumulator = red_op(accumulator, func(B[rows[ridx], col]))
    end

    accumulator
end

function createormask(inds::AbstractArray)
    mask = 0b11_11_11_11

    for qidx in inds
        mask -= 0b11 << 2*(qidx - 1)
    end

    mask
end

function column_mapreduce{T}(funcvec::Vector{T}, lastfuncvec::Vector{T}, red_op::Function, id_el::T, B::BEDMatrix, col::Integer)
    total = id_el

    X = B.X
    nbytes = B._byteheight

    for x in 1:(nbytes - 1)
        @inbounds total = red_op(total, funcvec[X[x, col] + 1])
    end
    @inbounds total = red_op(total, lastfuncvec[X[x, col] + 1])

    total
end

function column_hasNAs(B::BEDMatrix, col::Integer, rrange::AbstractUnitRange)
    X = B.X
    nbytes = B._byteheight
    na = NArep(B)
    rowstart, rowend = first(rrange), last(rrange)

    @inbounds begin
        firstquarters, bytestart, qstart = getquarterblock(B, rowstart, col)
        lastquarters, bytestop, qstop = getquarterblock(B, rowend, col)

        if qstart > 1
            stop = bytestart == bytestop ? qstop : 4
            for j in qstart:stop
                firstquarters[j] === na && return true
            end
            bytestart += 1
        end

        if bytestart <= bytestop
            if qstop < 4
                for j in 1:qstop
                    lastquarters[j] === na && return true
                end
                bytestop -= 1
            end

            for r in bytestart:bytestop
                hasNAs(X[r, col]) && return true
            end
        end
    end

    return false
end

"""
    column_hasNAs(B::BEDMatrix, col::Integer)

Returns `true` if there are missing values in `B[:, col]`.

"""
function column_hasNAs(B::BEDMatrix, col::Integer)
    X = B.X
    nbytes = B._byteheight

    for x in 1:nbytes
        @inbounds hasNAs(X[x, col]) && return true
    end

    return false
end

function countNAs(BC::BEDColumnUnitRange)
    B = parent(BC)
    X = B.X
    nacount = 0

    rrange, col = parentindexes(BC)
    rowstart, rowend = first(rrange), last(rrange)

    @inbounds begin
        bytestart, qstart = rowtobytequarter(rowstart)
        bytestop, qstop = rowtobytequarter(rowend)

        if qstart > 1
            stop = ifelse(bytestart == bytestop, qstop, 4)
            nacount += countNAs(maskbyte(X[bytestart, col], stop - qstart + 1, qstart - 1))
            bytestart += 1
        end

        if bytestart <= bytestop
            if qstop < 4
                nacount += countNAs(maskbyte(X[bytestop, col], qstop))
                bytestop -= 1
            end

            for x in bytestart:bytestop
                nacount += countNAs(X[x, col])
            end
        end
    end

    nacount
end

"""
    column_countNAs(B::BEDMatrix{T, S}, col::Integer)

Returns number of NAs in `B[:, col]`. Use `column_moments` if you need
 other moments.

"""
function column_countNAs(B::BEDMatrix, col::Integer)
    column_moments(B, col, Val{0})
end

# function column_moments{T, S, N}(B::BEDMatrix{T, S}, col::Integer, ::Type{Val{N}})
# end

function column_moments(B::BEDMatrix, col::Integer, ::Type{Val{0}})
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

function column_moments{T, S}(B::BEDMatrix{T, S}, col::Integer, ::Type{Val{1}})
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

function column_moments{T, S}(B::BEDMatrix{T, S}, col::Integer, ::Type{Val{2}})
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
    column_moments{T, S}(B::BEDMatrix{T, S}, col::Integer)

Returns number of NAs, sum of non-NAs, sum of non-NAs squared for
`B[:, col]`, in a computational efficient manner. These can be thought
of as the 0th, 1st, and 2nd moments of the column.

"""
function column_moments{T, S}(B::BEDMatrix{T, S}, col::Integer)
    column_moments(B, col, Val{2})
end

function column_dot{T, S}(B::BEDMatrix{T, S}, col::Integer, v::AbstractVector)
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

function column_dot{T, S}(B::BEDMatrix{T, S}, colL::Integer, colR::Integer)
    @boundscheck checkbounds(B, :, colL) && checkbounds(B, :, colR)

    dotsum = zero(promote_type(T, Int))
    X = B.X
    nbytes = B._byteheight

    for x in 1:(nbytes - 1)
        @inbounds dotsum += bytedot(X[x, colL], X[x, colR])
    end
    @inbounds dotsum += bytedot(X[nbytes, colL], X[nbytes, colR], B._lastrowSNPheight)

    dotsum
end

# https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance may be of interest

function colmean(B::BEDMatrix, col::Integer)
    nacount, total = column_moments(B, col, Val{1})

    total/(B.n - nacount)
end

function colmeanstd(B::BEDMatrix, col::Integer)
    nacount, total, totsq = column_moments(B, col)

    μ = total/(B.n - nacount)
    var = (totsq - 2μ*total + (B.n - nacount)*abs2(μ))/(B.n - nacount - 1)

    return μ, var
end



#################### BED byte to RAW math ####################

"""
    maskbyte(bedbyte::UInt8, num_quarters::Integer, qoffset::Integer=0)

Reads `byte` in BED format and returns `byte`, after masking quarters
such that `(qoffset + 1):(num_quarters + qoffset)` are preserved and
other quarters are "zeroed" (with respect to RAW format). Assumes RAW
format!!

"""
function maskbyte(bedbyte::UInt8, num_quarters::Integer, qoffset::Integer=0)
    @inbounds return (bedbyte | Consts.partialbytemasks[num_quarters, qoffset + 1])
end

"""
    hasNAs(bedbyte::UInt8)

Returns `true` if `bedbyte` has any missing values, and `false`
otherwise.

"""
function hasNAs(bedbyte::UInt8)
    @inbounds return Consts.hasNAmap[bedbyte + 1]
end

"""
    zeroNAs(bedbyte::UInt8)

Reads `bedbyte` in BED format and returns `bedbyte` with any missing
values set to zero. Assumes RAW format!!

"""
function zeroNAs(bedbyte::UInt8)
    @inbounds return Consts.natozeromap[bedbyte + 1]
end

"""
    countNAs(bedbyte::UInt8)

Returns the number of missing values in `bedbyte`.

"""
function countNAs(bedbyte::UInt8)
    @inbounds return Consts.nacountmap[bedbyte + 1]
end


function bytedot(b1::UInt8, b2::UInt8)
    @inbounds return Consts.bytebytemulttable[b1 + 1, b2 + 1]
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
    @boundscheck checkbounds(v, voffset + 4)

    dotsum = zero(T)
    quarters = breakbyte(zeroNAs(byte))

    for j in 1:4
        @inbounds dotsum += quarters[j]*v[j + voffset]
    end

    dotsum
end

"""
    bytedot(byte::UInt8, v::AbstractArray, voffset::Integer=0, num_quarters::Integer=4, qoffset::Integer=0)

Return RAW format dot product of `byte[(qoffset + 1):(num_quarters +
qoffset)]` in BED format with `v[(voffset + 1):(voffset +
num_quarters)]`. Missing values (in `byte`) are treated as zero.

"""
function bytedot{T}(byte::UInt8, v::AbstractArray{T}, voffset::Integer, num_quarters::Integer, qoffset::Integer=0)
    @boundscheck begin
        num_quarters + qoffset <= 4 || error("too many quarters")
        checkbounds(v, voffset + num_quarters)
    end

    quarters = breakbyte(zeroNAs(byte))
    dotsum = zero(T)

    for j in 1:num_quarters
        @inbounds dotsum += v[j + voffset]*quarters[j + qoffset]
    end

    dotsum
end


function byteabs2(b::UInt8)
    @inbounds return Consts.bytebytemulttable[b + 1, b + 1]
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

function bytedist(bedbyte::UInt8)
    @inbounds return Consts.distributionmap[bedbyte + 1]
end

"""
    bytedist(bedbyte::UInt8, num_quarters::Integer=4, qoffset::Integer=0)

Returns number of (RAW format) `0x00`s, `0x01`s, `0x02`s, and missing
values in `bedbyte`'s `(1 + qoffset):(num_quarters + qoffset)`
quarters.

!!! warning
    The excluded quarters are included in the `0x00` count. Thus if
    `num_quarters < 4`, one likely wants to use
    `bytedist(bedbyte)[1] - 4 + num_quarters` instead of `bytedist(bedbyte)[1]`.

"""
function bytedist(bedbyte::UInt8, num_quarters::Integer, qoffset::Integer=0)
    @boundscheck num_quarters + qoffset <= 4
    @inbounds dist = Consts.distributionmap[maskbyte(bedbyte, num_quarters, qoffset) + 1]
    return dist
end

"""
    quarterfuncvector{T}(func::Function, ::Type{T}=Int, func_na=0x0)

Convert `func` with return type `T` into a length-4 `Tuple`
corresponding to how `func` acts on the RAW representation of the BED
quarters `(0b00, 0b01, 0b10, 0b11)` converted to type `T`. `func_na`
should be the desired behavior on missing entries, defaults to
`zero(T)`. This fully specifies the behavior of `func` on a
`BEDMatrix{T, Q}`.

"""
quarterfuncvector{T}(func::Function, ::Type{T}=Int, func_na=0x00) = ntuple(q -> q < 4 ? func(convert(T, q - 1)) : convert(T, func_na), 4)


#################### SubArray Interface ####################

typealias Column{T, K, B} SubArray{T, 1, K, Tuple{Colon, Int}, B}
typealias BEDColumn{T, K <: BEDMatrix} Column{T, K, false}

typealias ColumnUnitRange{T, K, R <: Union{AbstractUnitRange, Colon}, B} SubArray{T, 1, K, Tuple{R, Int}, B}
typealias BEDColumnUnitRange{T, K <: BEDMatrix, R <: Union{AbstractUnitRange, Colon}} ColumnUnitRange{T, K, R, false}

typealias SubColumn{T, K, J, B} SubArray{T, 1, J, Tuple{J, Int}, B}
typealias BEDSubColumn{T, K <: BEDMatrix, J, B} SubArray{T, 1, J, Tuple{J, Int}, B}

# Strictly two dimensional view of BEDMatrix
typealias BEDSubMatrix{T, K <: BEDMatrix, R} SubArray{T, 2, K, R, false}

"""
    hasNAs(B)

Returns `true` if there are missing values in `B`.

"""
function hasNAs(B::BEDMatrix)
    for col in indices(B, 2)
        column_hasNAs(B, col) && return true
    end

    return false
end

function hasNAs(B::BEDSubMatrix)
    for col in indices(B, 2)
        hasNAs(view(B, :, col)) && return true
    end

    return false
end

hasNAs(BC::BEDColumn) = column_hasNAs(parent(BC), parentindexes(BC)[2])

function hasNAs(BC::BEDSubColumn)
    rows, col = parentindexes(BC)
    column_hasNAs(parent(BC), col, rows)
end

"""
    countNAs(B)

Returns the number of missing values in `B`. `B` may be a `BEDMatrix`
or a single column `SubArray` of a `BEDMatrix`.

"""
function countNAs(B::BEDMatrix)
    total = 0

    # parallizable
    for col in indices(B, 2)
        total += column_countNAs(B, col)
    end
    total
end

countNAs(BC::BEDColumn) = column_countNAs(parent(BC), parentindexes(BC)[2])

function countNAs(BC::BEDSubColumn)
    B = parent(BC)
    rows, col = parentindexes(BC)

    @inbounds return column_dist(B, col, rows)[4]
end

Base.sum(BC::BEDColumn) = column_sum(parent(BC), parentindexes(BC)[2])
function Base.sum(BC::BEDSubColumn)
    rows, col = parentindexes(BC)
    column_sum(parent(BC), col, rows)
end

"""
    sum(func::Function, BC::BEDSubColumn; skipna=true)

Applies `func` element-wise to `BC` and returns the sum. If one sets
`skipna=false`, then `func(NArep(parent(BC)))` is used for missing
values; otherwise missing values are not included.

"""
function Base.sum(func::Function, BC::BEDSubColumn; skipna=true)
    B = parent(BC)
    rows, col = parentindexes(BC)

    if skipna
        return column_sum(func, B, col, rows)
    end

    func_na = func(NArep(B))
    column_sum(func, func_na, B, col, rows)
end

# RAW format is positive:
Base.sumabs(BC::BEDSubColumn) = sum(BC)

Base.sumabs2(BC::BEDColumn) = column_sumabs2(parent(BC), parentindexes(BC)[2])
function Base.sumabs2(BC::BEDSubColumn)
    rows, col = parentindexes(BC)
    column_sumabs2(parent(BC), col, rows)
end

"""
    BEDdot(BCL::BEDColumn, BCR::BEDColumn)

Computes the dot product of two BEDMatrix columns, skipping missing
values.

"""
function BEDdot(BCL::BEDColumn, BCR::BEDColumn)
    BL = parent(BCL)
    BR = parent(BCR)
    if BL === BR
        return column_dot(BL, parentindexes(BCL)[2], parentindexes(BCR)[2])
    end

    column_dot(BL, parentindexes(BCL)[2], BR, parentindexes(BCR)[2])
end

"""
    BEDdot(BC::BEDColumn, v::Array)
    BEDdot(v::Array, BC::BEDColumn)

Computes the dot product of a BEDMatrix column `BC` with `v`, skipping
missing values.

"""
BEDdot(BC::BEDColumn, v::AbstractArray) = column_dot(parent(BC), parentindexes(BC)[2], v)
BEDdot(v::AbstractArray, BC::BEDColumn) = BEDdot(BC, v)

# ? Base.mapreduce
# Base.norm


#################### Column tools ####################

# ? column_any
# ? column_all

# Roughly equivalent to StatsBase.countmap
function column_dist(B::BEDMatrix, col::Integer, row::Integer)
    @boundscheck checkbounds(B, row, col)
    byterow, snpind = rowtobytequarter(row)
    @inbounds q = breakbyte(B.X[byterow, col])[snpind]

    if q === 0x00
        return (1, 0, 0, 0)
    elseif q === 0x01
        return (0, 1, 0, 0)
    elseif q === 0x02
        return (0, 0, 1, 0)
    end

    return (0, 0, 0, 1)
end

function column_dist(B::BEDMatrix, col::Integer, indices::AbstractVector{Bool})
    # Heuristic for converting to contiguous ranges
    # move the constants 0.9 and 100 to constants.jl
    if countnz(indices) > max(0.9*B.n, 100)
        return column_dist(B, col, tocontiguous(indices))
    end

    zero_count = 0
    one_count = 0
    two_count = 0
    na_count = 0
    X = B.X
    @inbounds for (row, b) in enumerate(indices)
        if b
            byterow, snpind = rowtobytequarter(row)
            q = breakbyte(X[byterow, col])[snpind]
            if q === 0x00
                zero_count += 1
            elseif q === 0x01
                one_count += 1
            elseif q === 0x02
                two_count += 1
            else
                na_count += 1
            end
        end
    end

    return (zero_count, one_count, two_count, na_count)
end

# TODO: Compare performance of two approaches, and settle on heuristic
function column_dist{T <: Integer}(B::BEDMatrix, col::Integer, indices::AbstractVector{T})
    # Heuristic for converting to contiguous ranges
    # move the constants 0.9 and 100 to constants.jl
    if length(indices) > max(0.9*B.n, 100)
        return column_dist(B, col, tocontiguous(indices))
    end

    zero_count = 0
    one_count = 0
    two_count = 0
    na_count = 0
    X = B.X
    @inbounds for row in indices
        byterow, snpind = rowtobytequarter(row)
        q = breakbyte(X[byterow, col])[snpind]
        if q === 0x00
            zero_count += 1
        elseif q === 0x01
            one_count += 1
        elseif q === 0x02
            two_count += 1
        else
            na_count += 1
        end
    end

    return (zero_count, one_count, two_count, na_count)
end

function column_dist{R <: AbstractUnitRange}(B::BEDMatrix, col::Integer, ranges::Vector{R})
    zero_count = 0
    one_count = 0
    two_count = 0
    na_count = 0

    for rrange in ranges
        counttup = column_dist(B, col, rrange)

        zero_count += counttup[1]
        one_count += counttup[2]
        two_count += counttup[3]
        na_count += counttup[4]
    end

    return (zero_count, one_count, two_count, na_count)
end

"""
    column_dist(B::BEDMatrix, col::Integer, rows=(:))

Returns a length-4 `Tuple` of the number of `0`s, `1`s, `2`s, and
`NArep(B)`s in `B[rows, col]`. `rows` may be a `UnitRange`, `Vector` of
`UnitRange`s, `Vector` of indices, or a logical index
(`Vector{Bool}`).

In the latter two cases, `rows` may be converted into
`Vector{UnitRange}` for better performance.

"""
function column_dist(B::BEDMatrix, col::Integer, rrange::AbstractUnitRange)
    bytestart, quarterstart = rowtobytequarter(first(rrange))
    bytestop, quarterstop = rowtobytequarter(last(rrange))

    unsafe_column_dist(B.X, col, bytestart, quarterstart, bytestop, quarterstop)
end

column_dist(B::BEDMatrix, col::Integer, ::Colon) = column_dist(B, col)
column_dist(B::BEDMatrix, col::Integer) = unsafe_column_dist(B.X, col, 1, 1, B._byteheight, B._lastrowSNPheight)

function unsafe_column_dist(X::Matrix{UInt8}, col::Integer, bytestart::Integer, quarterstart::Integer, bytestop::Integer, quarterstop::Integer)
    zero_count = 0
    one_count = 0
    two_count = 0
    na_count = 0

    @inbounds begin
        if quarterstart > 1
            num_quarters = ifelse(bytestop == bytestart, quarterstop - quarterstart + 1, 5 - quarterstart)

            bdist = bytedist(X[bytestart, col], num_quarters, quarterstart - 1)
            zero_count += (bdist[1] - 4 + num_quarters)  # correction for zeroed quarters
            one_count += bdist[2]
            two_count += bdist[3]
            na_count += bdist[4]

            bytestart += 1
        end

        if bytestart <= bytestop
            if quarterstop < 4
                bdist = bytedist(X[bytestop, col], quarterstop)
                zero_count += (bdist[1] - 4 + quarterstop)  # correction for zeroed quarters
                one_count += bdist[2]
                two_count += bdist[3]
                na_count += bdist[4]

                bytestop -= 1
            end

            @simd for x in bytestart:bytestop
                bdist = bytedist(X[x, col])
                zero_count += bdist[1]
                one_count += bdist[2]
                two_count += bdist[3]
                na_count += bdist[4]
            end
        end
    end

    return zero_count, one_count, two_count, na_count
end

"""
    column_sum{T, S}(B::BEDMatrix{T, S}, col::Integer, rows=(:))

Returns the sum of non-missing entries in `B[rows, col]`, optimized
for large `UnitRange`s of rows.

"""
function column_sum{T, S}(B::BEDMatrix{T, S}, col::Integer, rows=(:))
    counts = column_dist(B, col, rows)
    @inbounds return convert(promote_type(T, Int), counts[2] + 2*counts[3])
end

"""
    column_sumabs2{T, S}(B::BEDMatrix{T, S}, col::Integer, rows=(:))

Returns the sum of the squares of non-missing entries in `B[rows, col]`.

"""
function column_sumabs2{T, S}(B::BEDMatrix{T, S}, col::Integer, rows=(:))
    counts = column_dist(B, col, rows)
    @inbounds return convert(promote_type(T, Int), counts[2] + 4*counts[3])
end

"""
    column_norm(B::BEDMatrix, col::Integer, rows=(:), p=2)

Computes the `p`-norm of `B[rows, col]`, where the `p`-norm of a
vector `v` is given by ``\LaTeX``

```math
\\|v\\|_p = \\big(\sum_j |v_j|^p\\big)^\\frac{1}{p}
```

Missing values are ignored; use with care.

"""
function column_norm(B::BEDMatrix, col::Integer, rows=(:), p=2)
    if p == 2
        return sqrt(column_sumabs2(B, col, rows))
    end

    if p == 3
        return cbrt(column_sum((0, 1, 8, 0), B, col, rows))
    end

    @inbounds return (column_sum((0, 1, 2^p, 0), B, col, rows))^(1.0/p)
end

"""
    column_sum(func::Function, B::BEDMatrix, col::Integer, rows=(:))

Apply `func` to `B[:, col]`, and return the sum. `NArep(B)` values are
skipped.

"""
function column_sum{T, S}(f::Function, B::BEDMatrix{T, S}, col::Integer, rows=(:))
    column_sum(quarterfuncvector(f, T), B, col, rows)
end

"""
    column_sum(func::Function, func_na, B::BEDMatrix, col::Integer, rows)

Apply `func` to `B[rows, col]` and return the sum. `func_na` is
treated as the value of `func(NArep(B))`. (If omitted it is `0`.)

"""
function column_sum{T, S}(f::Function, func_na, B::BEDMatrix{T, S}, col::Integer, rows=(:))
    column_sum(quarterfuncvector(f, T, func_na), B, col, rows)
end

"""
    column_sum(func_qvec::Tuple, B::BEDMatrix, col::Integer, rows=(:))

`func_qvec` is a length-4 `Tuple` representation of a function `f`:
`func_qvec = (f(0), f(1), f(2), f(NArep(B)))`.  The result of applying
`f` to `B[rows, col]` and computing the sum is returned.

"""
function column_sum(func_qvec::Tuple, B::BEDMatrix, col::Integer, rows=(:))
    counts = column_dist(B, col, rows)
    @inbounds return tuple_dot(func_qvec, counts)
end

"""
    tuple_dot(t1::Tuple, t2::Tuple)

Compute scalar dot product of two `Tuple`s, treating them like `Vector`s.

"""
function tuple_dot(t1::Tuple, t2::Tuple)
    @boundscheck length(t1) == length(t2) || error("Tuples must have identical length")

    acc = 0
    for j in 1:length(t1)
        @inbounds acc += t1[j]*t2[j]
    end
    acc
end

"""
    column_hasNAs(B::BEDMatrix, col::Integer, rows=(:))

Returns `true` if `B[rows, col]` has missing values.

"""
function column_hasNAs(B::BEDMatrix, col::Integer, row::Integer)
    @inbounds return B[row, col] === NArep(B)
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

column_hasNAs(B::BEDMatrix, col::Integer, ::Colon) = column_hasNAs(B, col)

function column_hasNAs(B::BEDMatrix, col::Integer)
    X = B.X
    nbytes = B._byteheight

    for x in 1:nbytes
        @inbounds hasNAs(X[x, col]) && return true
    end

    return false
end

# Fall-back method
function column_hasNAs(B::BEDMatrix, col::Integer, rows)
    na = NArep(B)

    @inbounds for row in rows
        B[row, col] === na && return true
    end

    return false
end


"""
    column_countNAs(B::BEDMatrix{T, S}, col::Integer, rows=(:))

Returns number of NAs in `B[:, col]`.

"""
function column_countNAs(B::BEDMatrix, col::Integer)
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

column_countNAs(B::BEDMatrix, col::Integer, ::Colon) = column_countNAs(B, col)

# Fall-back method
function column_countNAs(B::BEDMatrix, col::Integer, rows)
    @inbounds return column_dist(B, col, rows)[4]
end


"""
    column_dot{T, S}(B::BEDMatrix{T, S}, col::Integer, v::AbstractVector)

Computes the scalar dot product between `B[:, col]` and `v`, with
missing values in `B` treated as zero.

"""
function column_dot{T, S}(B::BEDMatrix{T, S}, col::Integer, v::AbstractArray)
    @boundscheck B.n == length(v) || throw(DimensionMismatch("v has incorrect length"))

    dotsum = zero(promote_type(T, eltype(v), Int))  # Avoid overflow when T == UInt8

    X = B.X
    nbytes = B._byteheight

    for x in 1:(nbytes - 1)
        @inbounds dotsum += bytedot(X[x, col], v, 4*(x-1))
    end
    @inbounds dotsum += bytedot(X[nbytes, col], v, 4*(nbytes-1), B._lastrowSNPheight)

    dotsum
end

"""
    column_dot{T, S}(B::BEDMatrix{T, S}, colL::Integer, colR::Integer)

Computes the scalar dot product between `B[:, colL]` and `B[:, colR]`,
treating missing values as zero.

"""
function column_dot{T, S}(B::BEDMatrix{T, S}, colL::Integer, colR::Integer)
    @boundscheck begin
        checkbounds(B, :, colL)
        checkbounds(B, :, colR)
    end

    dotsum = zero(promote_type(T, Int))
    X = B.X
    nbytes = B._byteheight

    for x in 1:(nbytes - 1)
        @inbounds dotsum += bytedot(X[x, colL], X[x, colR])
    end
    @inbounds dotsum += bytedot(X[nbytes, colL], X[nbytes, colR], B._lastrowSNPheight)

    dotsum
end

"""
    column_dot{T, S, U, R}(BL::BEDMatrix{T, S}, colL::Integer, BR::BEDMatric{U, R}, colR::Integer)

Computes the scalar dot product between `BL[:, colL]` and
`BR[:, colR]`, treating missing values as zero.

"""
function column_dot{T, S, U, R}(BL::BEDMatrix{T, S}, colL::Integer, BR::BEDMatrix{U, R}, colR::Integer)
    @boundscheck begin
        checkbounds(BL, :, colL)
        checkbounds(BR, :, colR)
        BL.n == BR.n || throw(DimensionMismatch("BEDMatrices must have same height"))
    end

    dotsum = zero(promote_type(T, U, Int))
    XL = BL.X
    XR = BR.X
    nbytes = BL._byteheight

    for x in 1:(nbytes - 1)
        @inbounds dotsum += bytedot(XL[x, colL], XR[x, colR])
    end
    @inbounds dotsum += bytedot(XL[nbytes, colL], XR[nbytes, colR], BL._lastrowSNPheight)

    dotsum
end


"""
    column_mean(B::BEDMatrix, col::Integer, rows=(:))

Returns the mean of `B[rows, col]`, skipping missing values.

"""
function column_mean(B::BEDMatrix, col::Integer, rows=(:))
    counts = column_dist(B, col, rows)

    @inbounds return (counts[2] + 2.0*counts[3])/(B.n - counts[4])
end

# name matches `StatsBase.mean_and_var`, but interfaces does not
"""
    column_mean_and_var(B::BEDMatrix, col::Integer, rows=(:))

Computes the mean and (Bessel-corrected) sample variance of
`B[rows, col]`, skipping missing values.

"""
function column_mean_and_var(B::BEDMatrix, col::Integer, rows=(:))
    counts = column_dist(B, col, rows)

    @inbounds μ = (counts[2] + 2.0*counts[3])/(B.n - counts[4])
    @inbounds σ2 = (counts[1]*abs2(μ) + counts[2]*abs2(1.0 - μ) + counts[3]*abs2(2.0 - μ))/(B.n - counts[4] - 1)

    return μ, σ2
end

# name matches with `StatsBase.mean_and_std`, but interface does not
"""
    column_mean_and_std(B::BEDMatrix, col::Integer, rows=(:))

Computes the mean and (Bessel-corrected) sample standard deviation of
`B[rows, col]`, skipping missing values.

"""
function column_mean_and_std(B::BEDMatrix, col::Integer, rows=(:))
    μ, σ2 = column_mean_and_var(B, col, rows)
    return μ, sqrt(σ2)
end

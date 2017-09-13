# gwas.jl
# this file should ultimately be moved to a different package

## Reference implementation for no missing values, compare output with Base.linreg
# function univariate_ols{T}(x::Vector{T}, y::Vector{T})
#     xbar = mean(x)
#     ybar = mean(y)
#     n = length(x)
#
#     numerator = dot(x, y) - n*xbar*ybar
#     denominator = dot(x, x) - n*abs2(xbar)
#     b = numerator/denominator
#
#     R2 = dot(y, y) - n*abs2(ybar) - b*numerator
#     t = b*sqrt((n-2)*denominator/(n*R2))
#     pval = tTestpvalue(t, n - 2)
#
#     b, ybar - b*xbar, t, pval
# end

import Distributions

# two-sided test
function tTestpvalue(t::Real, dof::Real)
    dof <= zero(dof) && return 1.0

    # Compute using F-distribution test for t^2
    return Distributions.ccdf(Distributions.FDist(1, dof), abs2(t))
end

immutable UnivariateOLSFit
    column::Int
    n::Int
    y0::Float64
    β::Float64

    se::Float64
    t::Float64
    pval::Float64
end

colno(fit::UnivariateOLSFit) = fit.column
npoints(fit::UnivariateOLSFit) = fit.n
ndof(fit::UnivariateOLSFit) = fit.n - 2
intercept(fit::UnivariateOLSFit) = fit.y0
coef(fit::UnivariateOLSFit) = fit.β
stderr(fit::UnivariateOLSFit) = fit.se
tvalue(fit::UnivariateOLSFit) = fit.t
pvalue(fit::UnivariateOLSFit) = fit.pval

function Base.show(io::IO, fit::BEDMatrices.UnivariateOLSFit)
    # print(io, "(", fit.varname, ", ", fit.n, ", ", fit.y0, ", ", fit.β, ", ", fit.pval, ")")
    print(io, "col. ", colno(fit), ": ", "y = ", intercept(fit), ifelse(coef(fit) > 0, " + ", " - "), abs(coef(fit)), "*x", "; p = ", pvalue(fit))
end

function Base.show(io::IO, ::MIME"text/plain", fit::BEDMatrices.UnivariateOLSFit)
    print(io, "BEDMatrices.UnivariateOLSFit:\n")
    label = string("column ", colno(fit), " with ", npoints(fit), " entries\n")
    print(io, "  ", label)
    print(io, "  intercept: ", fit.y0, "\tslope: ", fit.β, "\n")
    print(io, "  std err: ", fit.se, "\tp-value: ", fit.pval)
end

function gwas_writecsv(fout::IO, fits::Vector{UnivariateOLSFit}; header::Bool=true, labels::Vector{String}=Vector{String}(0))
    if header
        if length(labels) > 0
            write(fout, "\"snp_id\",")
        end
        write(fout, "\"coef\",", "\"std err\",", "\"t-value\",", "\"Pr(>|t|)\"\n")
    end
    if length(labels) > 0
        for fit in fits
            @inbounds write(fout, "\"", labels[colno(fit)], "\",", string(coef(fit)), ",", string(stderr(fit)), ",", string(tvalue(fit)), ",", string(pvalue(fit)), "\n")
        end
    else
        for fit in fits
            write(fout, string(coef(fit)), ",", string(stderr(fit)), ",", string(tvalue(fit)), ",", string(pvalue(fit)), "\n")
        end
    end
end

"""
    gwas_writecsv(file, fits::Vector{UnivariateOLSFit}; header::Bool=true, labels::Vector{String}=Vector{String}(0))

Writes `fits` to `file` in CSV format. `file` may either be a filename
or a writeable `IO`. `header` determine whether a header is written to
`file`. `labels` should be a list of labels for the fits; otherwise
omitted.

"""
function gwas_writecsv(filename::AbstractString, fits::Vector{UnivariateOLSFit}; header::Bool=true, labels::Vector{String}=Vector{String}(0))
    open(filename, "w") do fout
        gwas_writecsv(fout, fits; header=header, labels=labels)
    end
end

sub_mean(y::AbstractArray, ::Colon) = mean(y)

function sub_mean{R <: AbstractUnitRange}(y::AbstractArray, rranges::Vector{R})
    n = 0
    acc = 0.0
    for rrange in rranges
        @simd for r in rrange
            @inbounds acc += y[r]
        end
        n += length(rrange)
    end

    acc/n
end
function sub_mean(y::AbstractArray, rows)
    acc = 0.0
    @simd for r in rows
        @inbounds acc += y[r]
    end

    acc/length(rows)
end

sub_dot(x, y, ::Colon) = dot(x, y)

function sub_dot{R <: AbstractUnitRange}(x::AbstractArray, y::AbstractArray, rranges::Vector{R})
    # @boundscheck
    acc = 0.0
    for rrange in rranges
        @simd for r in rrange
            @inbounds acc += x[r]*y[r]
        end
    end

    acc
end

function sub_dot(x, y, rows)
    # @boundscheck
    acc = 0.0
    @simd for r in rows
        @inbounds acc += x[r]*y[r]
    end

    acc
end


"""
    _column_olsfit{T}(B::BEDMatrix, col::Integer, rows, y::Vector{T}, ybar::T, y2::T)

Performs OLS fit of `y` over `B[:, col]` and returns `(β, y0, n, ssr,
varx)`: the univariate OLS slope, `β`; the intercept, `y0`; the number
of non-missing entries, `n`; the sum of squared residuals, `ssr`; and
the variance of `B[:, col]`, `varx`.

"""
function _column_olsfit{T}(B::BEDMatrix, col::Integer, rows, y::Vector{T}, ybar::T, y2::T)
    N_zeros, N_ones, N_twos, N_nas, x_y, nasup_y, nasup_y2 = BEDMatrices.column_dist_dot(B, col, y, rows)
    ntot = N_zeros + N_ones + N_twos + N_nas

    @inbounds begin
        n = ntot
        if N_nas > 0  # there are missing values
            # adjust n
            n -= N_nas

            # adjust ybar and y2 to the sample
            ybar = (ntot*ybar - nasup_y)/n
            y2 = y2 - nasup_y2

            # insufficient data for a fit
            n < 2 && return UnivariateOLSFit(col, n, ybar, 0.0, 0.0, 0.0, 1.0)
        end

        xsum = N_ones + 2*N_twos
        x2 = N_ones + 4*N_twos
        numerator = x_y - xsum*ybar

        # I keep the denominator::Integer by writing the formula this
        # way.
        denominator = n*x2 - abs2(xsum)

        # Handle some singular cases in a defensible way
        if denominator === 0
            if numerator !== zero(numerator)
                # all the points have same x value but different y values
                return UnivariateOLSFit(col, n, -sign(xsum)*sign(numerator)*Inf,
                                        sign(numerator)*Inf, Inf, 0.0, 1.0)
            else
                # all the points have same x AND y values!
                return UnivariateOLSFit(col, n, NaN, NaN, NaN, NaN, NaN)
            end
        end

        β = n*numerator/denominator
        y0 = ybar - β*xsum/n

        R2 = y2 - n*abs2(ybar) - β*numerator

        se = sqrt(n*R2/((n - 2)*denominator))
        t = β/se
        pval = tTestpvalue(t, n - 2)
    end

    return UnivariateOLSFit(col, n, y0, β, se, t, pval)
end


function column_olsfit(B, col::Integer, y::AbstractArray, rows=(:))
    @boundscheck begin
        length(y) == size(B, 1) || throw(DimensionMismatch("vector must have same length as size(B, 1)"))
        checkbounds(B, :, col)
    end

    ybar = sub_mean(y, rows)
    y2 = sub_dot(y, y, rows)

    _column_olsfit(B, col, rows, y, ybar, y2)
end

st_column_olsfit(B, y::AbstractArray, rows=(:)) = st_column_olsfit(B, y, rows, indices(B, 2))
function st_column_olsfit(B, y::AbstractArray, rows, cols)
    ybar = sub_mean(y, rows)
    y2 = sub_dot(y, y, rows)

    gwas_results = Vector{UnivariateOLSFit}(length(cols))

    for col in cols
        @inbounds gwas_results[col] = _column_olsfit(B, col, rows, y, ybar, y2)
    end

    return gwas_results
end

mt_column_olsfit(B, y::AbstractArray, rows=(:)) = mt_column_olsfit(B, y, rows, indices(B, 2))
function mt_column_olsfit(B, y::AbstractArray, rows, cols)
    ybar::Float64 = sub_mean(y, rows)
    y2 = convert(Float64, sub_dot(y, y, rows))

    gwas_results = Vector{UnivariateOLSFit}(length(cols))

    _mt_fill_gwas!(gwas_results, B, y, ybar, y2, rows, cols)

    return gwas_results
end

# workaround for https://github.com/JuliaLang/julia/issues/15276, and https://github.com/JuliaLang/julia/issues/17395, https://github.com/yuyichao/explore/blob/8d52fb6caa745a658f2c9bbffd3b0f0fe4a2cc48/julia/issue-17395/scale.jl#L21
@noinline function _mt_fill_gwas!(results, B, y::AbstractArray, ybar, y2, rows, cols)
    @inbounds begin
        Threads.@threads for col in cols
            results[col] = _column_olsfit(B, col, rows, y, ybar, y2)
        end
    end
    results
end

mp_column_olsfit(B::BEDMatrix, y::AbstractArray, rows=(:)) = mp_column_olsfit(B::BEDMatrix, y::AbstractArray, rows, indices(B, 2))
function mp_column_olsfit(B::BEDMatrix, y::AbstractArray, rows, cols)
    @boundscheck begin
        length(y) == size(B, 1) || throw(DimensionMismatch("vector must have same length as size(B, 1)"))
    end
    ybar = sub_mean(y, rows)
    y2 = sub_dot(y, y, rows)

    gwas_results = pmap(col -> _column_olsfit(B, col, rows, y, ybar, y2), cols)
end

# function mt_column_olsfit(B::BEDSubMatrix, y::AbstractArray)
# end
#
# function mp_column_olsfit(B::BEDSubMatrix, y::AbstractArray)
# end

"""
    GWAS(B::BEDMatrix, y::AbstractArray; mode::Symbol=:multithreaded, outfile::IO=DevNull, verbose::Bool=false, method::Symbol=:ols)

Returns vector of `UnivariateOLSFit` containers with fields `(n, y0, β, se, t, pval)`
 for fits of `y` over _each_ column of `B`.

## Optional keyword arguments
* `mode` determines parallelism; can be one of `:single`, `:multithreaded`, or [Not implemented] `:multiprocessed`
* `method` determines type of regression; currently only `:ols`, ordinary least squares, supported
* `outfile::IO` gives IO stream for writing the results
* `verbose::Bool` determines whether informative messages are printed

"""
function GWAS(B::BEDMatrix, y::AbstractArray; mode::Symbol=:multithreaded, outfile::IO=DevNull, verbose::Bool=true, method::Symbol=:ols)
    mode in (:single, :multithreaded, :multiprocessed) || error("Invalid mode, $mode")
    mode === :multiprocessed && error("Multiprocessing not currently implemented.")

    method in (:ols,) || error("Unsupported method, $method")

    if mode === :multithreaded
        if Threads.nthreads() === 1
            verbose && info("Multithreaded requested, but only a single thread available.",
                            " Number of threads fixed by environment variable `JULIA_NUM_THREADS`,",
                            " see julia documentation for details.")
            mode = :single
        elseif verbose
            info("Using ",Threads.nthreads(), " threads.")
        end
    end

    results = mode === :single ? st_column_olsfit(B, y) : (mode === :multithreaded ? mt_column_olsfit(B, y) : mp_column_olsfit(B, y))

    if outfile !== DevNull
        gwas_writecsv(outfile, results, labels=colnames(B))
    end

    results
end


##################################################
# For general ::Matrix (not ::BEDMatrix)
function mean_impute!{T}(M::AbstractMatrix{T}, na::T)
    n, p = size(M)

    z = zero(T)
    μ = z/one(T)

    na_inds = fill(0, n)

    @inbounds begin
        for col in indices(M, 2)
            # set up
            total = z
            na_count = 0

            # find NAs and total on non-NAs
            for row in indices(M, 1)
                x = M[row, col]
                if x === na
                    na_count += 1
                    na_inds[na_count] = row
                else
                    total += x
                end
            end

            # replace NAs with mean
            μ = total/(n - na_count)
            for na_idx in 1:na_count
                M[na_inds[na_idx], col] = μ

                # clean up as we go
                na_inds[na_idx] = 0
            end
        end
    end

    M
end

function _matrix_column_covvar{T}(M::AbstractMatrix{T}, col::Integer, y::AbstractArray, na::T)
    @inbounds begin
        na_count = 0
        xsum = zero(M[1, col])
        x2 = xsum
        ysum = zero(y[1])
        y2 = ysum
        x_y = zero(M[1, col]*y[1])
        for row in indices(M, 1)
            x = M[row, col]

            if x === na
                na_count += 1
            else
                xsum += x
                x2 += abs2(x)
                ysum += y[row]
                y2 += abs2(y)
                x_y += x*y[row]
            end
        end
    end
    return na_count, xsum, ysum, x2, y2, x_y
end

function _matrix_column_olsfit{T}(M::AbstractMatrix{T}, col::Integer, y::AbstractArray, n::Integer, na::T)
    @inbounds begin
        na_count, xsum, ysum, x2, y2, x_y = _matrix_column_covvar(M, col, y, na)
        n = size(M, 1) - na_count
        n < 2 && return UnivariateOLSFit(col, n, ysum/n, 0.0, 0.0, 0.0, 1.0)

        numerator = n*x_y - xsum*ysum
        denominator = n*x2 - abs2(xsum)

        β = numerator/denominator
        y0 = (ysum - β*xsum)/n

        R2 = y2 - n*abs2(ybar) - β*numerator/n
        se = sqrt(R2/((n-2)*denominator))
        t = β/se
        pval = tTestpvalue(t, n - 2)
    end

    return UnivariateOLSFit(col, n, y0, β, se, t, pval)
end

function st_matrix_column_olsfit{T}(M::AbstractMatrix{T}, y::AbstractArray{T}, na::T)
    gwas_results = Vector{UnivariateOLSFit}(size(M, 2))

    for col in 1:B.p
        @inbounds gwas_results[col] = _matrix_column_olsfit(M, col, y, length(y), na)
    end

    return gwas_results
end

function mt_matrix_column_olsfit{T}(M::AbstractMatrix{T}, y::AbstractArray{T}, na::T)
    n, p = size(M)

    gwas_results = Vector{UnivariateOLSFit}(p)

    _mt_matrix_fill_gwas!(gwas_results, M, y, n, p, na)
end

@noinline function _mt_matrix_fill_gwas!(results, M, y, n, p, na)
    @inbounds begin
        Threads.@threads for col in 1:p
            results[col] = _matrix_column_olsfit(M, col, y, n, na)
        end
    end

    results
end


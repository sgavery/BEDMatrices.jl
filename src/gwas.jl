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
    varname::String
    n::Int
    y0::Float64
    β::Float64

    se::Float64
    t::Float64
    pval::Float64
end

varname(fit::UnivariateOLSFit) = fit.varname
npoints(fit::UnivariateOLSFit) = fit.n
ndof(fit::UnivariateOLSFit) = fit.n - 2
intercept(fit::UnivariateOLSFit) = fit.y0
coef(fit::UnivariateOLSFit) = fit.β
stderr(fit::UnivariateOLSFit) = fit.se
tvalue(fit::UnivariateOLSFit) = fit.t
pvalue(fit::UnivariateOLSFit) = fit.pval

function gwas_writecsv(fout::IO, fits::Vector{UnivariateOLSFit}; header::Bool=true, labels::Bool=true)
    if header
        if labels
            write(fout, "\"snp_id\",")
        end
        write(fout, "\"coef\",", "\"std err\",", "\"t-value\",", "\"Pr(>|t|)\"\n")
    end
    if labels
        for fit in fits
            write(fout, "\"", varname(fit), "\",", string(coef(fit)), ",", string(stderr(fit)), ",", string(tvalue(fit)), ",", string(pvalue(fit)), "\n")
        end
    else
        for fit in fits
            write(fout, string(coef(fit)), ",", string(stderr(fit)), ",", string(tvalue(fit)), ",", string(pvalue(fit)), "\n")
        end
    end
end

"""
    gwas_writecsv(file, fits::Vector{UnivariateOLSFit}; header::Bool=true, labels::Bool=true)

Writes `fits` to `file` in CSV format. `file` may either be a filename
or a writeable `IO`. `header` and `labels` determine whether a header
and SNP labels are written to `file`, respectively.

"""
function gwas_writecsv(filename::AbstractString, fits::Vector{UnivariateOLSFit}; header::Bool=true, labels::Bool=true)
    open(filename, "w") do fout
        gwas_writecsv(fout, fits; header=header, labels=labels)
    end
end

"""
    _column_olsfit{T}(B::BEDMatrix, col::Integer, y::Vector{T}, ybar::T, y2::T)

Performs OLS fit of `y` over `B[:, col]` and returns `(β, y0, n, ssr,
varx)`: the univariate OLS slope, `β`; the intercept, `y0`; the number
of non-missing entries, `n`; the sum of squared residuals, `ssr`; and
the variance of `B[:, col]`, `varx`.

"""
function _column_olsfit{T}(B::BEDMatrix, col::Integer, y::Vector{T}, ybar::T, y2::T)
    @inbounds begin
        N_zeros, N_ones, N_twos, N_nas, x_y, nasup_y, nasup_y2 = BEDMatrices.column_dist_dot(B, col, y)
        ntot = n = B.n
        snp_name = colnames(B)[col]

        if N_nas > 0  # there are missing values
            # adjust n
            n -= N_nas

            # insufficient data for a fit
            n < 2 && return UnivariateOLSFit(snp_name, n, ybar, 0.0, 0.0, 0.0, 1.0)

            # adjust ybar and y2 to the sample
            ybar = (ntot*ybar - nasup_y)/n
            y2 = y2 - nasup_y2
        end

        xsum = N_ones + 2*N_twos
        x2 = N_ones + 4*N_twos
        numerator = x_y - xsum*ybar

        # I keep the denominator::Integer by writing the formula this
        # way.
        denominator = (n*x2 - abs2(xsum))

        β = n*numerator/denominator
        y0 = ybar - β*xsum/n

        R2 = y2 - n*abs2(ybar) - β*numerator

        se = sqrt(n*R2/((n - 2)*denominator))
        t = β/se
        pval = tTestpvalue(t, n - 2)
    end

    return UnivariateOLSFit(snp_name, n, y0, β, se, t, pval)
end

function column_olsfit(B::BEDMatrix, col::Integer, y::AbstractArray)
    @boundscheck begin
        length(y) == size(B, 1) || throw(DimensionMismatch("vector must have same length as size(B, 1)"))
        checkbounds(B, :, col)
    end

    ybar = mean(y)
    y2 = dot(y, y)

    _column_olsfit(B, col, y, ybar, y2)
end

function st_column_olsfit(B::BEDMatrix, y::AbstractArray)
    @boundscheck begin
        length(y) == size(B, 1) || throw(DimensionMismatch("vector must have same length as size(B, 1)"))
    end

    ybar = mean(y)
    y2 = dot(y, y)

    gwas_results = Vector{UnivariateOLSFit}(B.p)

    for col in 1:B.p
        @inbounds gwas_results[col] = _column_olsfit(B, col, y, ybar, y2)
    end

    return gwas_results
end

function mt_column_olsfit(B::BEDMatrix, y::AbstractArray)
    @boundscheck begin
        length(y) == size(B, 1) || throw(DimensionMismatch("vector must have same length as size(B, 1)"))
    end

    ybar = mean(y)
    y2 = dot(y, y)

    gwas_results = Vector{UnivariateOLSFit}(B.p)

    Threads.@threads for col in 1:B.p
        @inbounds gwas_results[col] = _column_olsfit(B, col, y, ybar, y2)
    end

    return gwas_results
end

function mp_column_olsfit(B::BEDMatrix, y::AbstractArray)
    @boundscheck begin
        length(y) == size(B, 1) || throw(DimensionMismatch("vector must have same length as size(B, 1)"))
    end
    ybar = mean(y)
    y2 = dot(y, y)

    gwas_results = pmap(col -> _column_olsfit(B, col, y, ybar, y2), 1:B.p)
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
        gwas_writecsv(outfile, results)
    end

    results
end

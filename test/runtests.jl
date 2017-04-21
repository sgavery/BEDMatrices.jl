## To do:
# * exceptions: no file, missing fam/bim, out-of-bounds indexing, not a bed file, old bed file


include("../BEDMatrices.jl")
using BEDMatrices
using Base.Test

const examplepath = "./data/"
const bedfile = "example.bed"
const rawfile = "example.raw"

const singleindices = [(13, 1), (13, 2), (13, 3), (13, 4),
                       (50, 998), (50, 999), (50, 1000),
                       (47, 3), (48, 3), (49, 3), (50, 3)]
const rowarrays = [[1, 11, 12], [1, 2], [48, 7, 49, 50]]
const colarrays = [[491, 13, 73], [1], [33, 31]]
const rowranges = [1:0, 1:2, 1:5, 1:50, 2:50, 3:50, 4:50]
const colranges = [1:0, 1:2, 1:101, 990:1000]
const linearranges = [1:75, 2:75, 50:51, 49:101, 40_000:50_000]


function readRAW(filename::String)
    n = countlines(filename) - 1

    lines = readlines(filename)

    header = split(lines[1])
    p = length(header) - 6
    colnames = header[7:end]
    data = Matrix{UInt8}(n, p)
    rownames = Vector{String}(n)

    for (row, line) in enumerate(lines[2:end])

        fields = split(line)
        rownames[row] = fields[1]*"_"*fields[2]
        data[row, :] = map(fld -> fld == "NA" ? BEDMatrices.NA_byte : parse(UInt8, fld), fields[7:end])
    end

    return rownames, colnames, data
end

examplerows, examplecols, exampledata = readRAW(examplepath*rawfile)

@testset "Basic Functional Testing" begin
    @testset "Reading BED file directly into matrix" begin
        @test BEDintomatrix(examplepath*bedfile) == exampledata
    end

    @testset "BEDMatrix" begin
        bed = BEDMatrix(examplepath*bedfile)
        @test bed[:, :] == exampledata

        @test path(bed) == abspath(examplepath*bedfile)

        @test rownames(bed) == examplerows
        @test colnames(bed) == examplecols
    end
end

@testset "Numeric Indexing" begin
    bed = BEDMatrix(examplepath*bedfile)
    @testset "single indices" begin
        for (r, c) in singleindices
            @test bed[r, c] == exampledata[r, c]
            @test bed[50*(c-1) + r] == exampledata[50*(c-1) + r]
        end
    end
    @testset "row ranges" begin
        for c in 1:1000
            for rrange in rowranges
                @test bed[rrange, c] == exampledata[rrange, c]
            end
            @test bed[4:end, c] == exampledata[4:end, c]
            @test bed[:, c] == exampledata[:, c]
        end
    end
    @testset "column ranges" begin
        for crange in colranges
            for r in 1:50
                @test bed[r, crange] == exampledata[r, crange]
            end
        end
        @test bed[17, 2:end] == exampledata[17, 2:end]
        @test bed[19, :] == exampledata[19, :]
    end
    @testset "linear ranges" begin
        for lrange in linearranges
            @test bed[lrange] == exampledata[lrange]
        end
    end
    @testset "column arrays" begin
        for cvec in colarrays
            @test bed[11, cvec] == exampledata[11, cvec]
        end
    end
    @testset "row arrays" begin
        for rvec in rowarrays
            @test bed[rvec, 11] == exampledata[rvec, 11]
        end
    end
    @testset "row arrays with column ranges" begin
        for crange in colranges
            for rvec in rowarrays
                @test bed[rvec, crange] == exampledata[rvec, crange]
            end
        end
    end
end

@testset "Name Indexing" begin
    bed = BEDMatrix(examplepath*bedfile)
    @testset "single indices" begin
        for (r, c) in singleindices
            rname = examplerows[r]
            cname = examplecols[c]
            @test bed[rname, c] == exampledata[r, c]
            @test bed[r, cname] == exampledata[r, c]
            @test bed[rname, cname] == exampledata[r, c]
        end
    end

    @testset "arrays" begin
        for cvec in colarrays
            cnames = [examplecols[c] for c in cvec]
            @test_broken bed[11, cnames] == exampledata[11, cvec]

            for rvec in rowarrays
                rnames = [examplerows[r] for r in rvec]
                @test_broken bed[rnames, cnames] == exampledata[rvec, cvec]
            end
        end
        for rvec in rowarrays
            rnames = [examplerows[r] for r in rvec]
            @test_broken bed[rnames, 11] == exampledata[rvec, 11]
        end
    end
end

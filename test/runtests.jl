include("../BEDMatrices.jl")
using Base.Test

const examplepath = "./data/"
const bedfile = "example.bed"
const rawfile = "example.raw"

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
        @test BEDMatrices.BEDintomatrix(examplepath*bedfile) == exampledata
    end

    @testset "BEDMatrix" begin
        bed = BEDMatrices.BEDMatrix(examplepath*bedfile)
        @test bed[:, :] == exampledata

        @test bed.path == abspath(examplepath*bedfile)
    end
end

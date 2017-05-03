## To do:
# * exceptions: no file, missing fam/bim, out-of-bounds indexing, not a bed file, old bed file
# * working directory

include("../BEDMatrices.jl")
using BEDMatrices
using Base.Test

const examplepath = "./data/"
const bedfile = "example.bed"
const rawfile = "example.raw"

const singleindices = [(13, 1), (13, 2), (13, 3), (13, 4),
                       (50, 998), (50, 999), (50, 1000),
                       (47, 3), (48, 3), (49, 3), (50, 3)]
const rowarrays = [[], collect(1:50), [1, 11, 12], [1, 2], [48, 7, 49, 50]]
const colarrays = [[], collect(1:1000), [491, 13, 73], [1], [33, 31]]
const rowranges = [1:0, 1:2, 1:5, 1:50, 2:50, 3:50, 4:50]
const colranges = [1:0, 1:2, 1:101, 990:1000]
const linearranges = [1:75, 2:75, 50:51, 49:101, 40_000:50_000]
const rowlogicals = [fill(true, 50), fill(false, 50), map(x -> x > 25, 1:50)]
const collogicals = [fill(true, 1000), fill(false, 1000), map(x -> x > 777, 1:1000)]

function simulatedBEDMatrix(n::Integer, p::Integer, datatype::DataType=UInt8, navalue=BEDMatrices.NA_byte, missingvalues=true)
    srand(0)  # for reproducibility
    if missingvalues
        M = rand(0b00:0b11, n, p)
    else
        M = rand([0b00, 0b10, 0b11], n, p)
    end

    byten = cld(n, 4)
    numlastsnps = n - 4*(byten - 1)

    X = Matrix{UInt8}(byten, p)
    for col in 1:p
        for r in 1:(byten - 1)
            X[r, col] = quarterstobyte(M[(4r - 3):(4r), col])
        end
        X[byten, col] = quarterstobyte([q <= numlastsnps ? M[end-numlastsnps + q, col] : 0b0 for q in 1:4])
    end

    path = "memory"
    rownames = [string("row_", j) for j in 1:n]
    colnames = [string("col_", j) for j in 1:p]

    bed = BEDMatrix{datatype, typeof(X)}(n, p, X, convert(datatype, navalue),
                                         path, colnames, rownames)
    navalue_typed = convert(datatype, navalue)
    data = map(x -> x === 0b01 ? navalue_typed : convert(datatype, BEDMatrices.rawformat(x)), M)
    return (data, bed)
end

function quarterstobyte(v::Vector{UInt8})
    byte = 0x0
    for (qidx, q) in enumerate(v)
        byte += q << 2*(qidx - 1)
    end
    byte
end

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

const examplerows, examplecols, exampledata = readRAW(examplepath*rawfile)
const bed = BEDMatrix(examplepath*bedfile)
const simbeds = [simulatedBEDMatrix(4, 5), simulatedBEDMatrix(5, 5, Int, 255), simulatedBEDMatrix(6, 5, Float64),
                 simulatedBEDMatrix(7, 5, Int8, -128), simulatedBEDMatrix(1, 5, Int8, -128, false), simulatedBEDMatrix(2, 1, Int8, -128, false),
                 simulatedBEDMatrix(3, 1, UInt, 0b11, false), simulatedBEDMatrix(4, 1, Int8, -128, false)]

@testset "Basic Functional Testing" begin
    @testset "Reading BED file directly into matrix" begin
        @test BEDintomatrix(examplepath*bedfile) == exampledata
    end

    @testset "BEDMatrix" begin
        @test bed[:, :] == exampledata

        @test path(bed) == abspath(examplepath*bedfile)

        @test rownames(bed) == examplerows
        @test colnames(bed) == examplecols
    end

    @testset "simulated BEDMatrix tests" begin
        # Mostly to check for possible mod 4 or type issues
        for (simdata, simbed) in simbeds
            @test simbed[:, :] == simdata
        end
    end
end

@testset "Numeric Indexing" begin
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
            @test bed[11, cnames] == exampledata[11, cvec]
            @test bed[:, cnames] == exampledata[:, cvec]
            @test bed[2:end, cnames] == exampledata[2:end, cvec]

            for rvec in rowarrays
                rnames = [examplerows[r] for r in rvec]
                @test bed[rnames, cnames] == exampledata[rvec, cvec]
            end
        end
        for rvec in rowarrays
            rnames = [examplerows[r] for r in rvec]
            @test bed[rnames, 11] == exampledata[rvec, 11]
            @test bed[rnames, :] == exampledata[rvec, :]
            @test bed[rnames, 2:end] == exampledata[rvec, 2:end]
        end
    end
end

@testset "Logical Indexing" begin
    for (r, c) in singleindices
        for clogic in collogicals
            @test bed[r, clogic] == exampledata[r, clogic]
        end
        for rlogic in rowlogicals
            @test bed[rlogic, c] == exampledata[rlogic, c]
        end
    end
    for clogic in collogicals
        for rlogic in rowlogicals
            @test bed[rlogic, clogic] == exampledata[rlogic, clogic]
        end
    end
end

@testset "Column Tools" begin
    @testset "typealiases" begin
        m = zeros(10, 20)

        v = view(bed, :, 1)
        w = view(bed, 1:25, 1)
        u = view(m, :, 1)
        t = view(m, 1:4, 1)

        @test isa(w, BEDMatrices.Column) == false
        @test isa(w, BEDMatrices.ColumnUnitRange) == true
        @test isa(w, BEDMatrices.BEDColumn) == false
        @test isa(w, BEDMatrices.BEDColumnUnitRange) == true

        @test isa(v, BEDMatrices.Column) == true
        @test isa(v, BEDMatrices.ColumnUnitRange) == true
        @test isa(v, BEDMatrices.BEDColumn) == true
        @test isa(v, BEDMatrices.BEDColumnUnitRange) == true

        @test isa(u, BEDMatrices.Column) == true
        @test isa(u, BEDMatrices.ColumnUnitRange) == true
        @test isa(u, BEDMatrices.BEDColumn) == false
        @test isa(u, BEDMatrices.BEDColumnUnitRange) == false

        @test isa(t, BEDMatrices.Column) == false
        @test isa(t, BEDMatrices.ColumnUnitRange) == true
        @test isa(t, BEDMatrices.BEDColumn) == false
        @test isa(t, BEDMatrices.BEDColumnUnitRange) == false
    end

    @testset "hasNAs" begin
        # extra test to confirm that simulatedBEDMatrix is working correctly
        nomissingdata, nomissingbed = simulatedBEDMatrix(50, 1000, Int, -128, false)
        @test any(e -> e == -128, nomissingdata) == false
        @test hasNAs(nomissingbed) == false

        @test hasNAs(bed) == true
        for (simdata, simbed) in simbeds
            @test hasNAs(simbed) == any(e -> e === NArep(simbed), simdata)
        end

        for col in 1:20
            @test hasNAs(view(bed, :, col)) == any(e -> e === BEDMatrices.NA_byte, exampledata[:, col])
            @test hasNAs(view(bed, 1:25, col)) == any(e -> e === BEDMatrices.NA_byte, exampledata[1:25, col])
        end

        @test hasNAs(view(bed, 1:6, 9)) == any(e -> e === BEDMatrices.NA_byte, exampledata[1:6, 9])
        @test hasNAs(view(bed, 1:7, 9)) == any(e -> e === BEDMatrices.NA_byte, exampledata[1:7, 9])
        @test hasNAs(view(bed, 8:50, 9)) == any(e -> e === BEDMatrices.NA_byte, exampledata[8:end, 9])
        @test hasNAs(view(bed, 7:50, 9)) == any(e -> e === BEDMatrices.NA_byte, exampledata[7:end, 9])
    end

    @testset "countNAs" begin
        @test countNAs(bed) == count(e -> e === BEDMatrices.NA_byte, exampledata)
        for col in 1:25
            @test countNAs(view(bed, :, col)) == count(e -> e === BEDMatrices.NA_byte, exampledata[:, col])
            @test countNAs(view(bed, 1:25, col)) == count(e -> e === BEDMatrices.NA_byte, exampledata[1:25, col])
        end
        for (simdata, simbed) in simbeds
            @test countNAs(simbed) == count(e -> e === NArep(simbed), simdata)
            @test countNAs(view(simbed, :, 1)) == count(e -> e === NArep(simbed), simdata[:, 1])
        end
    end

    @testset "sum" begin
        for col in 1:20
            @test sum(view(bed, :, col)) == sum(e -> e === BEDMatrices.NA_byte ? zero(e) : e, exampledata[:, col])
            @test sum(view(bed, 1:25, col)) == sum(e -> e === BEDMatrices.NA_byte ? zero(e) : e, exampledata[1:25, col])
        end
    end

    @testset "dot" begin
    end

    @testset "sumabs2" begin
        for col in 1:20
            @test sumabs2(view(bed, :, col)) == sum(e -> e === BEDMatrices.NA_byte ? zero(e) : abs2(e), exampledata[:, col])
            @test sumabs2(view(bed, 25:50, col)) == sum(e -> e === BEDMatrices.NA_byte ? zero(e) : abs2(e), exampledata[25:50, col])
        end
    end

    @testset "norm" begin
    end

    @testset "mapreduce" begin
    end

    @testset "moments" begin
    end
end

## To do:
# * exceptions: no file, missing fam/bim, out-of-bounds indexing, not a bed file, old bed file
# * working directory

include("../src/BEDMatrices.jl")
using BEDMatrices
using Base.Test
include("testtools.jl")

const examplepath = "./data/"
const bedfile = "example.bed"
const rawfile = "example.raw"

const singleindices = [(13, 1), (13, 2), (13, 3), (13, 4),
                       (50, 998), (50, 999), (50, 1000),
                       (47, 3), (48, 3), (49, 3), (50, 3)]
const rowarrays = [Vector{Int}(0), collect(1:50), [1, 11, 12], [1, 2], [48, 7, 49, 50]]
const colarrays = [Vector{Int}(0), collect(1:1000), [491, 13, 73], [1], [33, 31]]
const rowranges = [1:0, 1:2, 1:5, 1:50, 2:50, 3:50, 4:50]
const colranges = [1:0, 1:2, 1:101, 990:1000]
const linearranges = [1:75, 2:75, 50:51, 49:101, 40_000:50_000]
const rowlogicals = [fill(true, 50), fill(false, 50), map(x -> x > 25, 1:50)]
const collogicals = [fill(true, 1000), fill(false, 1000), map(x -> x > 777, 1:1000)]

const examplerows, examplecols, exampledata = readRAW(examplepath*rawfile)
const bed = BEDMatrix(examplepath*bedfile)
const simbeds = [simulatedBEDMatrix(4, 5), simulatedBEDMatrix(5, 5, Int, 255), simulatedBEDMatrix(6, 5, Float64),
                 simulatedBEDMatrix(7, 5, Int8, -128), simulatedBEDMatrix(1, 5, Int8, -128, false), simulatedBEDMatrix(2, 1, Int8, -128, false),
                 simulatedBEDMatrix(3, 1, UInt, 0b11, false), simulatedBEDMatrix(4, 1, Int8, -128, false),
                 simulatedBEDMatrix(13, 11, Int8, 127, true, true), simulatedBEDMatrix(13, 11, Float32, -128.0, false, true)]

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
            for j in eachindex(simdata)
                @test simbed[j] === simdata[j]
            end
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
            @test bed[rname, :] == convert(Vector{eltype(bed)}, exampledata[r, :])
            @test bed[:, cname] == convert(Vector{eltype(bed)}, exampledata[:, c])
        end
    end

    @testset "arrays" begin
        for cvec in colarrays
            cnames = [examplecols[c] for c in cvec]
            @test bed[11, cnames] == convert(Vector{eltype(bed)}, exampledata[11, cvec])
            @test bed[:, cnames] == convert(Matrix{eltype(bed)}, exampledata[:, cvec])
            @test bed[2:end, cnames] == convert(Matrix{eltype(bed)}, exampledata[2:end, cvec])

            for rvec in rowarrays
                rnames = [examplerows[r] for r in rvec]
                @test bed[rnames, cnames] == exampledata[rvec, cvec]
            end
        end
        for rvec in rowarrays
            rnames = [examplerows[r] for r in rvec]
            @test bed[rnames, 11] == convert(Vector{eltype(bed)}, exampledata[rvec, 11])
            @test bed[rnames, :] == convert(Matrix{eltype(bed)}, exampledata[rvec, :])
            @test bed[rnames, 2:end] == convert(Matrix{eltype(bed)}, exampledata[rvec, 2:end])
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

# @testset "tocontiguous" begin
#     for rrange in rowranges
#         rows = collect(rrange)
#         if isempty(rows)
#             @test BEDMatrices.tocontiguous(rows) == Vector{Int}(0)
#         else
#             @test BEDMatrices.tocontiguous(rows)[1] == rrange
#         end
#     end
# 
#     for rows in rowarrays[2:end]
#         @test mapreduce(collect, vcat, BEDMatrices.tocontiguous(rows)) == rows
#     end
# 
#     for rlogic in rowlogicals
#         if any(rlogic)
#             @test [x in mapreduce(collect, vcat, BEDMatrices.tocontiguous(rlogic)) for x in 1:50] == rlogic
#         else
#             @test BEDMatrices.tocontiguous(rlogic) == []
#         end
#     end
# end

@testset "Column Tools" begin
    @testset "typealiases" begin
        m = zeros(10, 20)

        v = view(bed, :, 1)
        w = view(bed, 1:25, 1)
        x = view(bed, [1, 2, 5], 1)

        u = view(m, :, 1)
        t = view(m, 1:4, 1)
        s = view(m, [1, 2, 5], 1)

        @test isa(v, BEDMatrices.Column) == true
        @test isa(v, BEDMatrices.ColumnUnitRange) == true
        @test isa(v, BEDMatrices.BEDColumn) == true
        @test isa(v, BEDMatrices.BEDColumnUnitRange) == true
        @test isa(v, BEDMatrices.BEDSubColumn) == true

        @test isa(w, BEDMatrices.Column) == false
        @test isa(w, BEDMatrices.ColumnUnitRange) == true
        @test isa(w, BEDMatrices.BEDColumn) == false
        @test isa(w, BEDMatrices.BEDColumnUnitRange) == true
        @test isa(w, BEDMatrices.BEDSubColumn) == true

        @test isa(x, BEDMatrices.Column) == false
        @test isa(x, BEDMatrices.ColumnUnitRange) == false
        @test isa(x, BEDMatrices.BEDColumn) == false
        @test isa(x, BEDMatrices.BEDColumnUnitRange) == false
        @test isa(x, BEDMatrices.BEDSubColumn) == true

        @test isa(u, BEDMatrices.Column) == true
        @test isa(u, BEDMatrices.ColumnUnitRange) == true
        @test isa(u, BEDMatrices.BEDColumn) == false
        @test isa(u, BEDMatrices.BEDColumnUnitRange) == false
        @test isa(u, BEDMatrices.BEDSubColumn) == false

        @test isa(t, BEDMatrices.Column) == false
        @test isa(t, BEDMatrices.ColumnUnitRange) == true
        @test isa(t, BEDMatrices.BEDColumn) == false
        @test isa(t, BEDMatrices.BEDColumnUnitRange) == false
        @test isa(t, BEDMatrices.BEDSubColumn) == false

        @test isa(s, BEDMatrices.Column) == false
        @test isa(s, BEDMatrices.ColumnUnitRange) == false
        @test isa(s, BEDMatrices.BEDColumn) == false
        @test isa(s, BEDMatrices.BEDColumnUnitRange) == false
        @test isa(s, BEDMatrices.BEDSubColumn) == false
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

        @test hasNAs(view(bed, 1:6, 7:10)) == any(e -> e === BEDMatrices.NA_byte, exampledata[1:6, 7:10])
        @test hasNAs(view(bed, 7:10, 1:20)) == any(e -> e === BEDMatrices.NA_byte, exampledata[7:10, 1:20])
    end

    @testset "countNAs" begin
        @test countNAs(bed) == count(e -> e === BEDMatrices.NA_byte, exampledata)

        for col in 1:25
            @test countNAs(view(bed, :, col)) == count(e -> e === BEDMatrices.NA_byte, exampledata[:, col])
            @test countNAs(view(bed, 1:25, col)) == count(e -> e === BEDMatrices.NA_byte, exampledata[1:25, col])
            @test countNAs(view(bed, [2, 3, 5, 7, 11, 13], col)) == count(e -> e === BEDMatrices.NA_byte, exampledata[[2, 3, 5, 7, 11, 13], col])
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
            @test sum(view(bed, [2, 3, 5, 7, 11], col)) == sum(e -> e === BEDMatrices.NA_byte ? zero(e) : e, exampledata[[2, 3, 5, 7, 11], col])

            @test sum(x -> (x-1)^10, view(bed, :, col)) == sum(e -> e === BEDMatrices.NA_byte ? zero(e) : (e - 1)^10, exampledata[:, col])
            @test sum(x -> x === NArep(bed) ? 37 : (x-1)^2, view(bed, :, col); skipna=false) == sum(e -> e === BEDMatrices.NA_byte ? 37 : (e - 1)^2, exampledata[:, col])
        end
    end

    @testset "column_sum(func, ...)" begin
        for col in 1:20
            @test column_sum(abs2, bed, col) == sum(e -> e === BEDMatrices.NA_byte ? zero(e) : abs2(e), exampledata[:, col])
            @test column_sum(x -> (x-1)^10, bed, col, 1:25) == sum(e -> e === BEDMatrices.NA_byte ? zero(e) : (e-1)^10, exampledata[1:25, col])
        end
    end

    @testset "BEDdot" begin
        v = collect(1:50)
        for col1 in 20:30
            @test BEDdot(view(bed, :, col1), v) == dot(map(e -> e === BEDMatrices.NA_byte ? zero(e) : e, exampledata[:, col1]), v)
            @test BEDdot(v, view(bed, :, col1)) == dot(map(e -> e === BEDMatrices.NA_byte ? zero(e) : e, exampledata[:, col1]), v)
            for col2 in 31:36
                @test BEDdot(view(bed, :, col1), view(bed, :, col2)) == dot(map(e -> e === BEDMatrices.NA_byte ? zero(e) : e, exampledata[:, col1]),
                                                                            map(e -> e === BEDMatrices.NA_byte ? zero(e) : e, exampledata[:, col2]))
            end
        end
    end

    @testset "column_dot" begin
        v = collect(1:50)
        for col1 in 20:30
            @test column_dot(bed, col1, v) == dot(map(e -> e === NArep(bed) ? zero(e) : e, bed[:, col1]), v)
            for col2 in 26:30
                @test column_dot(bed, col1, col2) == dot(map(e -> e === BEDMatrices.NA_byte ? zero(e) : e, exampledata[:, col1]),
                                                         map(e -> e === BEDMatrices.NA_byte ? zero(e) : e, exampledata[:, col2]))
            end
        end
    end

    @testset "column_dist" begin
        for col in 21:40
            dist = BEDMatrices.column_dist(bed, col)
            @test dist[1] == count(e -> e == 0, exampledata[:, col])
            @test dist[2] == count(e -> e == 1, exampledata[:, col])
            @test dist[3] == count(e -> e == 2, exampledata[:, col])
            @test dist[4] == count(e -> e === BEDMatrices.NA_byte, exampledata[:, col])
        end
    end

    @testset "column_dist_dot" begin
        v = collect(1:50)
        for col in 21:40
            distdot = BEDMatrices.column_dist_dot(bed, col, v)
            @test distdot[1] == count(e -> e == 0, exampledata[:, col])
            @test distdot[2] == count(e -> e == 1, exampledata[:, col])
            @test distdot[3] == count(e -> e == 2, exampledata[:, col])
            @test distdot[4] == count(e -> e === BEDMatrices.NA_byte, exampledata[:, col])
            @test distdot[5] == dot(map(e -> e === BEDMatrices.NA_byte ? zero(e) : e, exampledata[:, col]), v)
            @test distdot[6] == dot(map(e -> e === BEDMatrices.NA_byte ? one(e) : zero(e), exampledata[:, col]), v)
            @test distdot[7] == dot(map(e -> e === BEDMatrices.NA_byte ? one(e) : zero(e), exampledata[:, col]), v.^2)
        end
    end

    @testset "column_NAsup_dot" begin
        v = collect(1:50)
        for col in 1:20
            @test column_NAsup_dot(bed, col, v) == dot(map(e -> e === BEDMatrices.NA_byte ? one(e) : zero(e), exampledata[:, col]), v)
        end
    end


    @testset "sumabs2" begin
        for col in 1:20
            @test sumabs2(view(bed, :, col)) == sum(e -> e === BEDMatrices.NA_byte ? zero(e) : abs2(e), exampledata[:, col])
            @test sumabs2(view(bed, 25:50, col)) == sum(e -> e === BEDMatrices.NA_byte ? zero(e) : abs2(e), exampledata[25:50, col])
        end
    end

    @testset "column_norm" begin
        for col in 20:40
            @test column_norm(bed, col) == norm(map(e -> e === BEDMatrices.NA_byte ? zero(e) : e, exampledata[:, col]))
            @test column_norm(bed, col, :, 3) ≈ norm(map(e -> e === BEDMatrices.NA_byte ? zero(e) : e, exampledata[:, col]), 3)
            @test column_norm(bed, col, :, 4) ≈ norm(map(e -> e === BEDMatrices.NA_byte ? zero(e) : e, exampledata[:, col]), 4)
            @test column_norm(bed, col, 1:30, 4) ≈ norm(map(e -> e === BEDMatrices.NA_byte ? zero(e) : e, exampledata[1:30, col]), 4)
        end
    end
end

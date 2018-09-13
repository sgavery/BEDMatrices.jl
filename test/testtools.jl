# testtools.jl

import Random

function simulatedBEDMatrix(n::Integer, p::Integer, datatype::DataType=UInt8, navalue=BEDMatrices.NA_byte,
                            missingvalues=true, flips=false)
    Random.seed!(0)  # for reproducibility
    if missingvalues
        M = rand(0b00:0b11, n, p)
    else
        M = rand([0b00, 0b10, 0b11], n, p)
    end

    byten = cld(n, 4)
    numlastsnps = n - 4*(byten - 1)

    X = Matrix{UInt8}(undef, byten, p)
    for col in 1:p
        for r in 1:(byten - 1)
            X[r, col] = quarterstobyte(M[(4r - 3):(4r), col])
        end
        X[byten, col] = quarterstobyte([q <= numlastsnps ? M[end-numlastsnps + q, col] : 0b0 for q in 1:4])
    end

    path = "[memory]"
    rownames = [string("row_", j) for j in 1:n]
    colnames = [string("col_", j) for j in 1:p]
    if flips
        flipvec = BitVector(rand(Bool, p))
    else
        flipvec = falses(p)
    end

    bed = BEDMatrix{datatype, typeof(X)}(n, p, X, convert(datatype, navalue),
                                         path, colnames, rownames, BEDMatrices.getbytemap(convert(datatype, navalue)), flipvec)
    navalue_typed = convert(datatype, navalue)

    if flips
        data = Matrix{datatype}(undef, n, p)
        for col in 1:p
            if flipvec[col]
                for row in 1:n
                    x = M[row, col] === 0b01 ? navalue_typed : convert(datatype, BEDMatrices.rawformat(M[row, col]))
                    data[row, col] = x == navalue_typed ? x : 2 - x
                end
            else
                for row in 1:n
                    x = M[row, col] === 0b01 ? navalue_typed : convert(datatype, BEDMatrices.rawformat(M[row, col]))
                    data[row, col] = x
                end
            end
        end
    else
        data = map(x -> x === 0b01 ? navalue_typed : convert(datatype, BEDMatrices.rawformat(x)), M)
    end
    return (data, bed)
end

function quarterstobyte(v::Vector{UInt8})
    byte = 0x0
    for (qidx, q) in enumerate(v)
        byte += q << (2*(qidx - 1))
    end
    byte
end

function readRAW(filename::String)
    n = countlines(filename) - 1

    lines = readlines(filename)

    header = split(lines[1])
    p = length(header) - 6
    colnames = header[7:end]
    data = Matrix{UInt8}(undef, n, p)
    rownames = Vector{String}(undef, n)

    for (row, line) in enumerate(lines[2:end])
        fields = split(line)
        rownames[row] = fields[1]*"_"*fields[2]
        data[row, :] = map(fld -> fld == "NA" ? BEDMatrices.NA_byte : parse(UInt8, fld), fields[7:end])
    end

    return rownames, colnames, data
end

# constants.jl


#################### Helper function(s) ####################

function bytemask(num_quarters::Integer, qoffset::Integer)
    mask = 0b11_11_11_11

    for j in qoffset:(qoffset + num_quarters - 1)
        mask -= 0b11 << 2*j
    end

    mask
end

function flipsnp(x::UInt8)
    x = (~x) & 0b00_00_00_11
    (x << 1 | x >> 1) & 0b00_00_00_11
end

#################### Constants ####################


### .bed file constants
# Note: this is big endian [0x6c, 0x1b]
const plinkmagic = 0b0110_1100_0001_1011
const modes = Dict(0b0000_0001 => :SNPmajor,
                   0b0000_0000 => :SNPminor)

### BED format to Raw format mapping
const NA_byte = 0b11
const quarterstohuman = (0b10, NA_byte, 0b01, 0b00)

const flipBEDbytes = tuple([(flipsnp(snp4) << 6) | (flipsnp(snp3) << 4) | (flipsnp(snp2) << 2) | flipsnp(snp1) for
                            snp4 in 0b00:0b11 for
                            snp3 in 0b00:0b11 for
                            snp2 in 0b00:0b11 for
                            snp1 in 0b00:0b11]...)


"""
    const bytetoquarters::Tuple{Vararg{Tuple{UInt8,UInt8,UInt8,UInt8}, 256}}

A constant array storing all 256 bytes split into 4 quarters. Note
that the mapping is given by `bytetoquarters[BEDbyte + 1]` because of
julia's 1-indexing.

# Example
```julia
julia> BEDbyte = 0b11101001
0xe9
julia> BEDMatrices.bytetoquarters[BEDbyte + 1]
4-element Array{UInt8,1}:
 0x03
 0x01
 0x01
 0x00
```

This can be understood as breaking the byte into quarters, in reverse
order: `0b01`, `0b10`, `0b10`, `0b11`, and then changing from BED
format to `rawformat` with `0x03` representing missing value.

"""
const bytetoquarters = tuple([(quarterstohuman[(byte & 0b00_00_00_11) + 1],
                               quarterstohuman[(byte & 0b00_00_11_00) >> 2 + 1],
                               quarterstohuman[(byte & 0b00_11_00_00) >> 4 + 1],
                               quarterstohuman[(byte & 0b11_00_00_00) >> 6 + 1]) for
                              byte in 0x00:0xff]...)

"""
    const defaultbytemap::Matrix{Tuple{Vararg{UInt8, 4}}}

This is a matrix of `size(defaultbytemap) == (256, 2)`. The first
columne is the same as `bytetoquarters`: `defaultbytemap[:, false + 1]
== collect(bytetoquarters)`, see docstring for `bytetoquarters`. The
second column corresponds to doing a major--minor allele flip.

"""
const defaultbytemap = hcat([(quarterstohuman[(byte & 0b00_00_00_11) + 1],
                              quarterstohuman[(byte & 0b00_00_11_00) >> 2 + 1],
                              quarterstohuman[(byte & 0b00_11_00_00) >> 4 + 1],
                              quarterstohuman[(byte & 0b11_00_00_00) >> 6 + 1]) for
                             byte in 0x00:0xff],
                            [(quarterstohuman[(byte & 0b00_00_00_11) + 1],
                              quarterstohuman[(byte & 0b00_00_11_00) >> 2 + 1],
                              quarterstohuman[(byte & 0b00_11_00_00) >> 4 + 1],
                              quarterstohuman[(byte & 0b11_00_00_00) >> 6 + 1]) for
                             byte in flipBEDbytes])

### BED byte to RAW math
const onebyte = 0b10_10_10_10  # vector of ones in RAW format

const partialbytemasks = [bytemask(numq, qoff) for numq in 1:4, qoff in 0:3]  # tuple?
# const partialbytemasks = [BEDMatrices.bytemask(numq, qoff) for numq in 1:4, qoff in 0:3]

const inversebytemap = Dict(bytetoquarters[b + 1] => b for b in 0x0:0xff)
const natozeromap = [inversebytemap[map(q -> q == NA_byte ? 0b00 : q, bytetoquarters[b + 1])] for b in 0x0:0xff]  # tuple?
const na_sup_map = [inversebytemap[map(q -> q == NA_byte ? 0b01 : 0b00, bytetoquarters[b + 1])] for b in 0x0:0xff]  # tuple?

const nacountmap = [count(q -> q === NA_byte, bytetoquarters[b + 1]) for b in 0x0:0xff]  # tuple?
const distributionmap = tuple(map(b -> (count(q -> q === 0b00, bytetoquarters[b + 1]),
                                        count(q -> q === 0b01, bytetoquarters[b + 1]),
                                        count(q -> q === 0b10, bytetoquarters[b + 1]),
                                        count(q -> q === 0b11, bytetoquarters[b + 1])), 0x0:0xff)...)

const hasNAmap = tuple([c > 0 for c in nacountmap]...)

const bytebytemulttable = [dot(collect(bytetoquarters[natozeromap[b1 + 1] + 1]),
                               collect(bytetoquarters[natozeromap[b2 + 1] + 1])) for b1 in 0x0:0xff, b2 in 0x0:0xff]


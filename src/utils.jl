"""
    NON_DEGEN_BASES

Store a sorted tuple of the four standard, non-degenerate nucleotide bases (`A`, `C`, `G`, `T`).
"""
const NON_DEGEN_BASES = "ACGT" |> collect |> sort |> Tuple

"""
    DEGEN_BASES

Store a sorted tuple of the IUPAC degenerate nucleotide ambiguity codes (`B`, `D`, `H`, `K`, `M`, `N`, `R`, `S`, `V`, `W`, `Y`).
"""
const DEGEN_BASES = "MRWSYKVHDBN" |> collect |> sort |> Tuple

"""
    ALL_BASES

Store a sorted tuple of all valid nucleotide bases, including both non-degenerate and degenerate IUPAC codes.
"""
const ALL_BASES = [NON_DEGEN_BASES..., DEGEN_BASES...] |> sort |> Tuple

"""
    BASES_W_GAPS

Store a tuple of all valid nucleotide bases (non-degenerate and degenerate) along with the gap character (`-`).
"""
const BASES_W_GAPS = (ALL_BASES..., '-')

"""
    IUPAC_B2V

Map each IUPAC nucleotide character to a tuple of its constituent non-degenerate bases (e.g., `'R' => ("A", "G")`).
"""
const IUPAC_B2V = Dict(
    'A'=>Tuple("A"),  'C'=>Tuple("C"),  'G'=>Tuple("G"),  'T'=>Tuple("T"),
    'R'=>Tuple("AG"), 'Y'=>Tuple("CT"), 'S'=>Tuple("CG"), 'W'=>Tuple("AT"),
    'K'=>Tuple("GT"), 'M'=>Tuple("AC"), 'B'=>Tuple("CGT"),'D'=>Tuple("AGT"),
    'H'=>Tuple("ACT"),'V'=>Tuple("ACG"),'N'=>Tuple("ACGT")
)

"""
    IUPAC_V2B

Map a tuple of non-degenerate bases to its corresponding IUPAC ambiguity character (the inverse of [`IUPAC_B2V`](@ref)).
"""
const IUPAC_V2B = Dict(v=>k for (k,v) in IUPAC_B2V)

"""
    IUPAC_COUNTS

Map each IUPAC nucleotide character to the number of unique non-degenerate bases it represents (e.g., `'R' => 2`, `'N' => 4`).
"""
const IUPAC_COUNTS = Dict(k=>length(v) for (k,v) in IUPAC_B2V)

"""
    IUPAC_GC_CONTENT

Map each IUPAC nucleotide character to its expected GC content as a fraction (e.g., `'S' => 1.0`, `'W' => 0.0`).
"""
const IUPAC_GC_CONTENT = Dict(k=>count(in("CG"), v)/length(v) for (k,v) in IUPAC_B2V)

"""
    DNA_COMP_TABLE_DEG

Store a 256-element lookup table (`Vector{UInt8}`) for fast, branchless complementation of ASCII DNA characters, including degenerate IUPAC codes and gaps (`-`).
"""
const DNA_COMP_TABLE_DEG = let 
    _comp_dna_deg(b::UInt8)::UInt8 =
    b == UInt8('A') ? UInt8('T') : b == UInt8('a') ? UInt8('T') :
    b == UInt8('T') ? UInt8('A') : b == UInt8('t') ? UInt8('A') :
    b == UInt8('C') ? UInt8('G') : b == UInt8('c') ? UInt8('G') :
    b == UInt8('G') ? UInt8('C') : b == UInt8('g') ? UInt8('C') :
    b == UInt8('R') ? UInt8('Y') : b == UInt8('r') ? UInt8('Y') :
    b == UInt8('Y') ? UInt8('R') : b == UInt8('y') ? UInt8('R') :
    b == UInt8('S') ? UInt8('S') : b == UInt8('s') ? UInt8('S') :
    b == UInt8('W') ? UInt8('W') : b == UInt8('w') ? UInt8('W') :
    b == UInt8('K') ? UInt8('M') : b == UInt8('k') ? UInt8('M') :
    b == UInt8('M') ? UInt8('K') : b == UInt8('m') ? UInt8('K') :
    b == UInt8('B') ? UInt8('V') : b == UInt8('b') ? UInt8('V') :
    b == UInt8('D') ? UInt8('H') : b == UInt8('d') ? UInt8('H') :
    b == UInt8('H') ? UInt8('D') : b == UInt8('h') ? UInt8('D') :
    b == UInt8('V') ? UInt8('B') : b == UInt8('v') ? UInt8('B') :
    b == UInt8('N') ? UInt8('N') : b == UInt8('n') ? UInt8('N') :
    UInt8('-')
    [_comp_dna_deg(UInt8(i)) for i in 0:255]
end

"""
    IUPAC_PROBS

Map each IUPAC nucleotide character to a 4-tuple of probabilities representing the fractional likelihood of `A`, `C`, `G`, and `T`, respectively.
"""
const IUPAC_PROBS = let
    d = Dict{Char, NTuple{4, Float64}}()
    for (k, v) in IUPAC_B2V
        probs = zeros(Float64, 4)
        for base in v
            idx = findfirst(==(base), NON_DEGEN_BASES)
            probs[idx] += 1.0
        end
        if !isempty(v)
            probs ./= length(v)
        end
        d[k] = Tuple(probs)
    end
    d
end

"""
    MAX_GC_OPTIONS

Map each IUPAC nucleotide character to a tuple of its constituent bases, filtered to prioritize `G` and `C` for maximizing GC content during sampling.
"""
const MAX_GC_OPTIONS = let
    d = Dict{Char, Tuple}()
    for (k0, v0) in IUPAC_B2V
        _v = filter(in("GC"), v0)
        v = isempty(_v) ? v0 : _v
        d[k0] = v
    end
    d
end

"""
    MIN_GC_OPTIONS

Map each IUPAC nucleotide character to a tuple of its constituent bases, filtered to prioritize `A` and `T` for minimizing GC content during sampling.
"""
const MIN_GC_OPTIONS = let
    d = Dict{Char, Tuple}()
    for (k0, v0) in IUPAC_B2V
        _v = filter(!in("GC"), v0)
        v = isempty(_v) ? v0 : _v
        d[k0] = v
    end
    d
end
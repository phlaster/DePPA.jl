# DePPA.jl

**DE**generate **P**rimer **P**air **A**ssembler

`DePPA.jl` is a high-performance pure Julia package for multiple sequence alignment (MSA) analysis and PCR primer design. The package is specifically engineered to handle degenerate (IUPAC) nucleotide sequences and provides rigorous statistical calculations for the thermodynamic properties of primer pools.

## Features

- **Strict Oligonucleotide Typing:** Distinct, type-stable structures for pure, degenerate, and gapped sequences, with zero-allocation sequence slicing.
- **MSA Analysis:** Efficient parsing of FASTA alignments, calculation of position-specific metrics (depth, determinacy), and generation of consensus sequences. Includes a customizable terminal viewer.
- **Thermodynamic Primer Design:** Automated MSA scanning and multi-criteria filtering. Primers are evaluated based on the statistical distribution of $T_m$ and $\Delta G$ across their degenerate variants.
- **Lazy Thermodynamic Engine:** Core alignment and sequence manipulation requires no external dependencies. Nearest-neighbor (NN) thermodynamic calculations are seamlessly integrated via a Julia Extension linking to SeqFold.jl, loading the heavy dependencies only when required.

## Installation

To use the core alignment and sequence features, install `DePPA.jl`:

```julia
using Pkg
Pkg.add("DePPA")
```

For alignment and thermodynamic calculations, you also need to load the optional extensions:

```julia
Pkg.add("MAFFT_jll") # in-place multiple sequence alignment engine
Pkg.add("SeqFold")   # thermodynamic calculations for primer design

using DePPA, MAFFT_jll, SeqFold
```

## Quick Start

First, we define a small set of unaligned sequences in memory and write them to a temporary FASTA file.

```julia
julia> data = """>1
              GATCTGTAATGAGCGGCAGACCGACCGCGAATTAGACCTCGCCGAAGCCCTG
              GCCGCCAAGCTCAATTCGAAGCTCATTCACTTCGTGCCGCGCGAC
              >2
              CATTTGCAACGAGCGTCAGACCGACCGCGAACTCGACCTGGCCGAAGCGCTG
              GCTGCCAAACTCAATTCTAAGCTCATCCACTTCGTGCCACGC
              >3
              CATTTGTAACGAGCGTCAGACCGACCGTGAACTCGACCTCGCCGAAGCGCTG
              GCTGCCAAACTCAATTCCAAGCTCATCCACTTCGTGCCACGCGACAA
              >4
              CTGTAACGAGCGGCAGACTGACCGAGAATTAGACCTCGCTGAAGCGCTGGCC
              GCCAAGCTCAATTCGAAGCTCATTCACTTTGTGCCGCGCGACAACA
              >5
              TGTAACGAGCGGCAGACTGACCGAGAATTAGACCTCGCTGAAGCGCTGGCCG
              CCAAGCTCAATTCGAAGCTCATTCACTTTGTGCCGCGCGACAACA""";

julia> temp_file = tempname(); open(temp_file, "w") do f write(f, data) end;
```

Next, we construct an `MSA` object. Passing `mafft=true` triggers the MAFFT engine (via `MAFFT_jll`) to align the sequences in-place. The terminal output automatically displays a truncated, color-coded view with sequence depths and a consensus track.

```julia
julia> alignment = MSA(temp_file; mafft=true)
MSA with 5 sequences of length 101:
   -ATTTGTAACGAGCGGCAGACCGACCGAGAATTAGACCTCGCCGAAGCGCTGGCCGCCAAGCTCAATTCGAA…
1 >G..C.....T.................C....................C..........................…
2 >C.....C........T...........C...C.C.....G..............T.....A........T.....…
3 >C..............T...........T...C.C....................T.....A........C.....…
4 >.--C.................T....................T................................…
5 >.---.................T....................T................................…
   1        ⋅         ⋅         ⋅         ⋅         ⋅         ⋅         ⋅   75
```

Now we generate candidate primers. The `construct_primers` function scans the alignment, filters candidates based on GC content, $T_m$, and $\Delta G$ distributions, and returns a list of valid primers. We do this for both forward (`is_forward=true`) and reverse (`is_forward=false`) primers.

```julia
julia> fwd = construct_primers(alignment); first(fwd)
Constructing F... 100%|██████████| Time: 0:00:01
Forward degenerate primer with 2 deg. positions
                                             \50:66>
|==============================================================================101|

  Sequence: CTGGCYGCCAARCTCAA
  Length: 17
  Positions: 50:66
  Unique variants: 4
  Melting temperature: 56.0°C (55.7⋅56.9°C)
  Min ΔG: 1.09 kcal/mol
  GC content: 58.8%
  Description: "Degenerate consensus for 5 seq MSA"

julia> rev = construct_primers(alignment; is_forward=false); first(rev)
Constructing R... 100%|██████████| Time: 0:00:00
Reverse degenerate primer with 2 deg. positions

|==============================================================================101|
                                 <50:66\
  Sequence: TTGAGYTTGGCRGCCAG
  Length: 17
  Positions: 50:66
  Unique variants: 4
  Melting temperature: 56.0°C (55.4⋅57.3°C)
  Min ΔG: 1.14 kcal/mol
  GC content: 58.8%
  Description: "Reverse complement of Degenerate consensus for 5 seq MSA"
```

Finally, we pair the forward and reverse primers. `best_pairs` matches primers based on the desired amplicon length and $T_m$ compatibility, returning the best combinations sorted by the smallest $T_m$ difference.

```julia
julia> bp = best_pairs(fwd, rev; amplicon_len=50:51); first(bp)
PCR primer pair for 5 seq. MSA, amplicon: 18:68 (51bp)
              >_________________51bp_________________<                             
|==============================================================================101|
Forward: AGACYGACCGHGAAYTMGACCT at 18:39
Reverse: AATTGAGYTTGGCRGCCA at 51:68
Tm: 55.8±0.1 °C
```
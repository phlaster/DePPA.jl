# <div align="center"> <img src="docs/src/assets/logo.png" alt="DEPPA.jl: DEgenerate Primer Pair Assembler" width="500"></div><div align="center">Degenerate Primer Pair Assembler</div>

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://phlaster.github.io/DEPPA.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://phlaster.github.io/DEPPA.jl/dev/)
[![Build Status](https://github.com/phlaster/DEPPA.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/phlaster/DEPPA.jl/actions/workflows/CI.yml?query=branch%3Amaster)
[![Aqua](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
[![codecov](https://codecov.io/gh/phlaster/DEPPA.jl/graph/badge.svg?token=DCH8TMMXOA)](https://codecov.io/gh/phlaster/DEPPA.jl)
[![License](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
## Introduction

`DEPPA.jl` is a high-performance pure Julia package for multiple sequence alignment (MSA) analysis and PCR primer design. The package is specifically engineered to handle degenerate (IUPAC) nucleotide sequences and provides rigorous statistical calculations for the thermodynamic properties of primer pools.

## Motivation

Standard bioinformatics tools often treat degenerate positions by selecting a consensus base or calculating thermodynamics for a single "average" sequence. Biophysically, however, a degenerate primer is a mixture of distinct oligonucleotides, each with its own unique $T_m$ and $\Delta G$. 

`DEPPA.jl` addresses this by treating degenerate primers as statistical ensembles. Instead of a single point estimate, the package calculates the thermodynamic distribution of the primer pool. This allows researchers to quantitatively assess and mitigate amplification bias before entering the laboratory.

## Key Features

*   **Strict Oligonucleotide Typing:** Distinct, type-stable structures for pure, degenerate, and gapped sequences, with zero-allocation sequence slicing.
*   **MSA Analysis:** Efficient parsing of FASTA alignments, calculation of position-specific metrics (depth, determinacy), and generation of consensus sequences. Includes a customizable terminal viewer.
*   **Thermodynamic Primer Design:** Automated MSA scanning and multi-criteria filtering. Primers are evaluated based on the statistical distribution of $T_m$ and $\Delta G$ across their degenerate variants.
*   **Lazy Thermodynamic Engine:** Core alignment and sequence manipulation requires no external dependencies. Nearest-neighbor (NN) thermodynamic calculations (using parameters from SantaLucia, 2004 and Turner, 2009) are seamlessly integrated via a Julia Extension linking to [SeqFold.jl](https://github.com/phlaster/SeqFold.jl), loading the heavy dependencies only when required.

## Comparison with Existing Solutions

While several tools exist for PCR primer design, they often differ significantly in their handling of degenerate sequences, licensing models, and integration capabilities. 

Below is a factual comparison of `DEPPA.jl` against standard open-source libraries, free tools, and commercial suites.

| Package | Degenerate Primer Design | License | Ecosystem & Integration | API & Documentation |
| :--- | :--- | :--- | :--- | :--- |
| **Primer3** | Not supported natively; single sequence only. | LGPL/GPL | Standalone C library; requires wrappers. | Config-file driven; legacy docs. |
| **OpenPrimeR** | Specialized for degenerate pools. | GPL | R-only; difficult external integration. | Shiny GUI focused; secondary API. |
| **Geneious Prime** | GUI-based; native IUPAC support. | Commercial | Closed Java plugin; paid license. | GUI-centric; secondary API. |
| **CLC Genomics** | GUI-based; native IUPAC support. | Commercial | Closed ecosystem; limited scripting. | GUI-centric; manual-driven. |
| **DEPPA.jl** | Native ensemble thermodynamics; statistical distributions. | MIT | Native Julia; seamless Python interop. | Modern, type-stable API; inline docs. |

## Installation

To use the core alignment and sequence features, install `DEPPA.jl`:
```julia
julia> ]

pkg> add DEPPA
pkg> add MAFFT_jll # in-place multiple sequence alignment requires the engine
pkg> add SeqFold # enable thermodynamic calculations for primer design

julia> using DEPPA.Alignments, DEPPA.Primers
julia> using SeqFold, MAFFT_jll
julia> setMSAShowStyle!(:bw); # monochrome REPL output theme, also try `:polymorf` or `:allcolors`
```

## Basic Usage

### Quickstart

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

Next, we construct an `MSA` object. Passing `mafft=true` triggers the MAFFT engine (via `MAFFT_jll`) to align the sequences in-place. The terminal output automatically displays a truncated view with a consensus track.

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

Now we generate candidate primers. The `construct_primers` function scans the alignment, filtering candidates based on reasonable defaults for $GC$ content, $T_m$, and $\Delta G$ distributions, and returns a list of valid primers. We do this for both forward (`is_forward=true`) and reverse (`is_forward=false`) primers.

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

## Thermodynamic Engine and Extensions

`DEPPA.jl` utilizes Julia's native package extension system. The computational lifting for nearest-neighbor thermodynamics is offloaded to `SeqFold.jl`. 

When you call functions like `SeqFold.tm(primer)` or `construct_primers(...)`, the `SeqFoldExt` module intercepts these calls. It expands the degenerate primer into its non-degenerate variants, calculates the thermodynamics for each variant independently, and aggregates the results into a statistical distribution.

If `SeqFold.jl` is not installed in your environment, the core `DEPPA` modules will still load and function perfectly for sequence parsing and MSA manipulation, but thermodynamic functions will throw an informative error prompting you to install the extension.

## Architecture and Performance

Designed for modern multi-core systems, `DEPPA.jl` leverages Julia's native multithreading for computationally intensive tasks such as MSA scanning, primer candidate generation, and Monte Carlo sampling of large degenerate pools. 

The use of `OligoView` and `MSAView` ensures that slicing and subsetting operations are strictly $O(1)$ memory operations, allowing the package to handle massive metagenomic alignments without triggering garbage collection bottlenecks.

## Calling from Python

`DEPPA.jl` can be seamlessly integrated into Python bioinformatics pipelines using [`juliacall`](https://pypi.org/project/juliacall/).

Install the bridge via your preferred Python package manager:
```bash
$ uv init
$ uv add juliacall
$ uv run python
```

In your Python script or REPL:
```python
Python 3.13.13 (main, Apr 14 2026, 14:28:56) [Clang 22.1.3 ] on linux
Type "help", "copyright", "credits" or "license" for more information.
>>> from juliacall import Main as jl
>>> jl.seval("""using Pkg; Pkg.add(url="https://github.com/phlaster/DEPPA.jl"); Pkg.add("SeqFold")""")
>>> jl.seval("using DEPPA.Oligos, SeqFold")
# Define convenient wrapper
>>> jl.seval('calc_tm(seq::String) = tm(DegenOligo(seq))')
# Call as a Python function
>>> result = jl.calc_tm("AGACYGACCGHGAAYTMGACCT")
>>> print(f"Mean Tm: {result.mean}, Confidence: {result.conf}")
Mean Tm: 55.7, Confidence: (44.2, 65.6)
```

*Note: As with any Julia-Python bridge, the first function call will incur a one-time latency due to Julia's JIT compilation. Subsequent calls will execute at native speeds.*

## Citation

If you use `DEPPA.jl` in your research, please cite:

```bibtex
@misc{DEPPA.jl,
  author       = {A.D. Bezlepsky},
  title        = {{DEPPA.jl: DEgenerate Primer Pair Assembler}},
  year         = {2026},
  publisher    = {GitHub},
  journal      = {GitHub repository},
  howpublished = {\url{https://github.com/phlaster/DEPPA.jl}},
}
```
# <div align="center"> <img src="docs/src/assets/logo.png" alt="DePPA.jl: Degenerate Primer Pair Assembler" width="500"></div><div align="center">Degenerate Primer Pair Assembler</div>

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://phlaster.github.io/DePPA.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://phlaster.github.io/DePPA.jl/dev/)
[![Build Status](https://github.com/phlaster/DePPA.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/phlaster/DePPA.jl/actions/workflows/CI.yml?query=branch%3Amaster)
[![Aqua](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
[![codecov](https://codecov.io/gh/phlaster/DePPA.jl/graph/badge.svg?token=DCH8TMMXOA)](https://codecov.io/gh/phlaster/DePPA.jl)
[![License](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

## Introduction

`DePPA.jl` is a high-performance, pure Julia package for multiple sequence alignment (MSA) analysis and **degenerate PCR primer design**. It natively handles IUPAC degenerate sequences and provides rigorous statistical calculations for the thermodynamic properties of complex primer pools.

## Why DePPA.jl?

Most standard bioinformatics tools (like **Primer3**) are designed exclusively for *single, pure sequences* and cannot natively process MSAs. Commercial suites (Geneious, CLC) offer basic IUPAC support but often rely on simplistic consensus algorithms that ignore the thermodynamic reality of mixed primer pools.

`DePPA.jl` takes an **MSA as its primary input**, identifies conserved regions, and constructs degenerate primers capable of amplifying entire gene families simultaneously. 

### Ensemble Thermodynamics
To ensure successful PCR amplification across all target variants, `DePPA.jl` evaluates the statistical distribution of thermodynamic parameters for the *entire primer pool*, rather than a single "average" sequence:
*   **$T_m$ (Melting Temperature):** Calculates the distribution of melting temperatures across all variants to prevent amplification bias.
*   **$\Delta G$ (Gibbs Free Energy):** Evaluates the stability of the primer-template duplex for all variants, ensuring no sequence forms excessively stable secondary structures or binds too weakly.

The computational engine for nearest-neighbor (NN) thermodynamics is powered by extending the functionality of [SeqFold.jl](https://github.com/phlaster/SeqFold.jl). While `SeqFold` natively calculates properties for non-degenerate sequences, `DePPA.jl` expands degenerate oligos into their non-degenerate variants, computes the thermodynamics for each independently, and aggregates the results into a statistical distribution. The underlying NN parameters and complex salt corrections are based on the foundational biophysical models of SantaLucia & Hicks (2004) and Owczarzy et al. (2008).

## Comparison with Existing Solutions

While several tools exist for PCR primer design, they fundamentally differ in their handling of alignments and degenerate sequences. 

<table style="width:100%; border-collapse: collapse; text-align: left;">
  <thead>
    <tr style="background-color: #f6f8fa; border-bottom: 2px solid #dfe2e5;">
      <th style="padding: 10px; border: 1px solid #dfe2e5;">Package</th>
      <th style="padding: 10px; border: 1px solid #dfe2e5; text-align: center;">Works with MSA?</th>
      <th style="padding: 10px; border: 1px solid #dfe2e5;">Degenerate Primer Design</th>
      <th style="padding: 10px; border: 1px solid #dfe2e5;">License</th>
      <th style="padding: 10px; border: 1px solid #dfe2e5;">Ecosystem & Integration</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td style="padding: 10px; border: 1px solid #dfe2e5;"><b>Primer3</b></td>
      <td style="padding: 10px; border: 1px solid #dfe2e5; text-align: center;">No (Single seq)</td>
      <td style="padding: 10px; border: 1px solid #dfe2e5;">Not supported natively</td>
      <td style="padding: 10px; border: 1px solid #dfe2e5;">LGPL/GPL</td>
      <td style="padding: 10px; border: 1px solid #dfe2e5;">Standalone C library; requires wrappers.</td>
    </tr>
    <tr>
      <td style="padding: 10px; border: 1px solid #dfe2e5;"><b>OpenPrimeR</b></td>
      <td style="padding: 10px; border: 1px solid #dfe2e5; text-align: center;">Yes</td>
      <td style="padding: 10px; border: 1px solid #dfe2e5;">Limited (Heuristic)</td>
      <td style="padding: 10px; border: 1px solid #dfe2e5;">GPL</td>
      <td style="padding: 10px; border: 1px solid #dfe2e5;">R-only; difficult external integration.</td>
    </tr>
    <tr>
      <td style="padding: 10px; border: 1px solid #dfe2e5;"><b>Geneious Prime</b></td>
      <td style="padding: 10px; border: 1px solid #dfe2e5; text-align: center;">Yes</td>
      <td style="padding: 10px; border: 1px solid #dfe2e5;">Basic (Consensus only)</td>
      <td style="padding: 10px; border: 1px solid #dfe2e5;">Commercial</td>
      <td style="padding: 10px; border: 1px solid #dfe2e5;">Closed Java plugin; paid license.</td>
    </tr>
    <tr>
      <td style="padding: 10px; border: 1px solid #dfe2e5;"><b>CLC Genomics</b></td>
      <td style="padding: 10px; border: 1px solid #dfe2e5; text-align: center;">Yes</td>
      <td style="padding: 10px; border: 1px solid #dfe2e5;">Basic (Consensus only)</td>
      <td style="padding: 10px; border: 1px solid #dfe2e5;">Commercial</td>
      <td style="padding: 10px; border: 1px solid #dfe2e5;">Closed ecosystem; limited scripting.</td>
    </tr>
    <tr style="background-color: #eaf4ff; font-weight: bold; border: 2px solid #0969da;">
      <td style="padding: 10px; border: 1px solid #dfe2e5;"><img src="docs/src/assets/logo.png" alt="DePPA.jl" width="120"></td>
      <td style="padding: 10px; border: 1px solid #dfe2e5; text-align: center;">Yes (Native)</td>
      <td style="padding: 10px; border: 1px solid #dfe2e5;">Full Ensemble Thermodynamics</td>
      <td style="padding: 10px; border: 1px solid #dfe2e5;">MIT</td>
      <td style="padding: 10px; border: 1px solid #dfe2e5;">Native Julia; seamless Python interop.</td>
    </tr>
  </tbody>
</table>

## Installation

```julia
julia> ]
pkg> add DePPA MAFFT_jll SeqFold
```
> **`MAFFT_jll` is used for automatic alignment, and `SeqFold` serves as the backend for nearest-neighbor thermodynamic calculations.**

## REPL Experience & Quickstart

`DePPA.jl` is designed to be highly interactive directly from the Julia REPL. It features rich, color-coded terminal output for alignments and primers.

```julia
julia> using DePPA.Alignments, DePPA.Primers, DePPA.Oligos

julia> using SeqFold, MAFFT_jll

# Optional: set the REPL visualization style (:bw, :polymorf, or :allcolors)
julia> setMSAShowStyle!(:bw); 
```

### 1. Load and Align Sequences
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

### 2. Construct Degenerate Primers
The engine scans the alignment, filtering candidates based on $GC$ content, $T_m$, and $\Delta G$ distributions.

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

### 3. Pair Primers
`best_pairs` matches primers based on amplicon length and $T_m$ compatibility.

```julia
julia> bp = best_pairs(fwd, rev; amplicon_len=50:51); first(bp)
PCR primer pair for 5 seq. MSA, amplicon: 18:68 (51bp)
              >_________________51bp_________________<                             
|==============================================================================101|
Forward: AGACYGACCGHGAAYTMGACCT at 18:39
Reverse: AATTGAGYTTGGCRGCCA at 51:68
Tm: 55.8±0.1 °C
```

# Python Integration

`DePPA.jl` can be seamlessly integrated into Python bioinformatics pipelines using [`juliacall`](https://pypi.org/project/juliacall/).

```bash
$ uv add juliacall
$ uv run python
```
In Python REPL:
```python
>>> from juliacall import Main as jl
>>> jl.seval("""using Pkg; Pkg.add("DePPA"); Pkg.add("SeqFold")""")
>>> jl.seval("using DePPA.Oligos, SeqFold")
# Define convenient wrapper
>>> jl.seval('calc_tm(seq::String) = tm(DegenOligo(seq))')
>>> result = jl.calc_tm("AGACYGACCGHGAAYTMGACCT")
>>> print(f"Mean Tm: {result.mean}, Confidence: {result.conf}")
Mean Tm: 55.7, Confidence: (44.2, 65.6)
```

## Roadmap

*   **Primer Export:** Convenient export of constructed primers to standard formats (FASTA, CSV) for direct ordering.
*   **Specificity Checks:** In silico verification of primer binding against the target matrix to filter out non-specific amplification.
*   **Documentation:** Expanding the documentation with more comprehensive usage examples and tutorials.

## References

*   **SantaLucia, J., & Hicks, D. (2004)**. Thermodynamics of DNA-RNA interactions and DNA-DNA interactions. *Annual Review of Biophysics and Biomolecular Structure*, 33, 415-440.
*   **Owczarzy, R., Moreira, B. G., You, Y., Behlke, M. A., Walder, J. A., & Walder, J. (2008)**. Effects of sodium, magnesium, and spermidine on the stability of DNA duplexes. *Biochemistry*, 47(20), 5336-5353.
*   <a href="https://github.com/phlaster/SeqFold.jl" target="_blank"><img src="https://raw.githubusercontent.com/phlaster/SeqFold.jl/0ee91b0601645fba350643b9fe767dd8d89a0f90/docs/src/assets/logo.png" alt="SeqFold.jl Logo" width="120" align="middle"></a>: The underlying nearest-neighbor thermodynamic engine.

## Citation

If you use `DePPA.jl` in your research, please cite:

```bibtex
@misc{DePPA.jl,
  author       = {A.D. Bezlepsky},
  title        = {{DePPA.jl: Degenerate Primer Pair Assembler}},
  year         = {2026},
  publisher    = {GitHub},
  journal      = {GitHub repository},
  howpublished = {\url{https://github.com/phlaster/DePPA.jl}},
}
```
# API Reference

This page provides a complete list of the public API exported by `DEPPA.jl` and its submodules. 

## DEPPA

The top-level module orchestrates the package's subcomponents and provides a unified interface for degenerate primer design.

```@docs
DEPPA
```

All usefull functions are typically obtained by loading submodules. Some functions are useless without loading 3rd-party packages like `MAFFT_jll`.

## Oligos

The `Oligos` module provides a strict, type-stable hierarchy for representing nucleic acid sequences. It defines distinct concrete types—[`Oligo`](@ref), [`DegenOligo`](@ref), and [`GappedOligo`](@ref)—to handle pure, IUPAC degenerate, and gapped sequences, respectively. By subtyping `AbstractString`, these structures integrate seamlessly with Julia's standard string processing while enabling zero-allocation slicing via [`OligoView`](@ref). The module also includes utilities for expanding degenerate sequences into their non-degenerate variants, either through complete enumeration or Monte Carlo sampling.

```@autodocs
Modules = [DEPPA.Oligos]
Public = true
Order = [:type, :function]
```

## Alignments

The `Alignments` module is dedicated to the construction, visualization, and statistical analysis of Multiple Sequence Alignments (MSAs). Central to this module is the [`MSA`](@ref) type, which precomputes base frequencies and supports bootstrap resampling for robust consensus generation. Using [`MSAView`](@ref), users can subset alignments into submatrices of rows and columns with $O(1)$ memory overhead, facilitating efficient analysis of large metagenomic datasets. The module also provides functions to calculate position-specific metrics like depth and determinacy, generate consensus sequences, and filter out poorly aligned regions.

```@autodocs
Modules = [DEPPA.Alignments]
Public = true
Order = [:type, :function]
```

## Primers

The `Primers` module automates the design of degenerate PCR primers directly from an MSA. The [`construct_primers`](@ref) function performs multithreaded scanning of the alignment, evaluating candidates against strict thermodynamic and conservation filters. Unlike traditional tools that evaluate a single consensus sequence, this module treats degenerate primers as statistical ensembles, calculating distributions for melting temperature ($T_m$), free energy ($\Delta G$), and GC content across all non-degenerate variants. Finally, [`best_pairs`](@ref) matches forward and reverse primers based on amplicon length and thermodynamic compatibility.

```@autodocs
Modules = [DEPPA.Primers]
Public = true
Order = [:type, :function]
```
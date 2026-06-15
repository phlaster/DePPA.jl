# DEPPA - DEgenerate Primer Pair Assembler

[![Build Status](https://github.com/phlaster/DEPPA.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/phlaster/DEPPA.jl/actions/workflows/CI.yml?query=branch%3Amaster)
[![Aqua](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

# Quickstart

```julia
using DEPPA
using MAFFT_jll, SeqFold # for full features

file = "test/assets/TruA.fasta.gz"
aln = MSA(file; mafft=true, bootstrap=10)

fwds = construct_primers(aln)
revs = construct_primers(aln; is_forward=false)

ppairs = best_pairs(fwds, revs; amplicon_len=190:190)

open("example.txt", "w") do f
    for pp in ppairs
        show(f, MIME"text/plain"(), pp)
        println(f)
    end
end
```
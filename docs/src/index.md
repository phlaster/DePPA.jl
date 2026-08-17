# DePPA.jl

**D**egenerate **P**rimer **P**air **A**ssembler

`DePPA.jl` is a high-performance pure Julia package for multiple sequence alignment (MSA) analysis and PCR primer design. The package is specifically engineered to handle degenerate (IUPAC) nucleotide sequences and provides rigorous statistical calculations for the thermodynamic properties of primer pools.

## Features

- **Strict Oligonucleotide Typing:** Distinct, type-stable structures for pure, degenerate, and gapped sequences, with zero-allocation sequence slicing.
- **MSA Analysis:** Efficient parsing of FASTA alignments, calculation of position-specific metrics (depth, determinacy), and generation of consensus sequences. Includes a customizable terminal viewer.
- **Thermodynamic Primer Design:** Automated MSA scanning and multi-criteria filtering. Primers are evaluated based on the statistical distribution of $T_m$ and $\Delta G$ across their degenerate variants.
- **Specificity Filtering:** Evaluates candidate primers against the input MSA to reject off-target binding sites, using a probabilistic threshold (`offtarget_reject_threshold`) for both forward and reverse complement matches.
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

using DePPA.Alignments, DePPA.Primers, MAFFT_jll, SeqFold
```

## Quick Start

### Loading data
Here we define a small set of unaligned sequences in memory and write them to a temporary FASTA file. But you will use the actual fasta file from the disc.

```julia
julia> data = """>1
               GATCTGTAATGAGCGGCAGACCGACCGCGAATTAGACCTCGC
               CGAAGCCCTGGCCGCCAAGCTCAATTCGAAGCTCATTCAC
               >2
               CATTTGCAACGAGCGTCAGACCGACCGCGAACTCGACCTGGC
               CGAAGCGCTGGCTGCCAAACTCAATTCTAAGCTCATC
               >3
               CATTTGTAACGAGCGTCAGACCGACCGTGAACTCGACCTCGC
               CGAGCTGCCAAACTCAATTCCAAGCTCATCCACTT
               >4
               CTGTAACGAGCGGCAGACTGACCGAGAATTAGACCTCGCTGA
               AGCGCTGGCCGCCAAGCTCAATTCGAAGCTCATTCACTTTG
               >5
               TGTAACGAGCGGCAGACTGACCGAGAATTAGACCTCGCTGAA
               GCGCTGGCCGCCAAGCTCAATTCGAAGCTCATTCACTTAG""";

julia> temp_fasta = tempname(); open(temp_fasta, "w") do f write(f, data) end;
```

Next, we construct an `MSA` object. Passing `mafft=true` triggers the MAFFT engine (via `MAFFT_jll`) to align the sequences in-place. Your terminal output automatically displays a truncated, color-coded view with a consensus track.

```julia
julia> alignment = MSA(temp_fasta; mafft=true)
MSA with 5 sequences of length 86:
   CATCTGTAACGAGCGGCAGACCGACCGAGAATTAGACCTCGCCGAAGCGCTGGCCGCCAAGCTCAATTCGAAGCTCATTCACTT--
1 >G..C.....T.................C....................C.................................----
2 >C.....C........T...........C...C.C.....G..............T.....A........T........C-------
3 >C..............T...........T...C.C...........-------..T.....A........C........C.....--
4 >---C.................T....................T.........................................TG
5 >----.................T....................T.........................................AG
   1        ·         ·         ·         ·         ⌢         ·         ·         ·    86
```

### Main method

Now we generate candidate primers. The `construct_primers` function scans the alignment and returns a single vector containing both forward and reverse primers, filtering candidates based on GC content, $T_m$, $\Delta G$ distributions, and specificity (off-target binding within the MSA). The custom `Base.show` method provides a histogram of primer distribution across the alignment.

```julia
julia> primers = construct_primers(alignment)
Primer distribution for 5 seq MSA (L=86): 94 primers
32  ┤                                                    ▃▅▅▅▆██████▆▃                 
    │                                                  ▁▆█████████████▇▂▂▂             
    │                 ▂▂▂▂▂▄▅    ▁▁▁▁▁▁▁              ▅███████████████████▆▁           
    │              ▂▄▇████████▇█████████▇▃          ▂▇██████████████████████▄          
    │   ▁▃▅▇██████████████████████████████▇▇▇▇▇▇▇▅▃▅██████████████████████████▅▂▂▂▂    
0   ┼ ╞════════════════════════════════════════════════════════════════════════════════
    │           ▔▔▔▔▔▔█████████████████████▀▀▀▀▀▀▔▔▔████████████████████████▀▀▔▔       
    │                 ▔▔▔▔▔▀▀▀▀▀█████████▀▔         ▔████████████████████▀▔▔           
    │                           ▔▔▔▔▔▔▔▔▔            ▔█████████████████▀▔              
    │                                                 ▔▔▔▀███████████▀▔                
32  ┤                                                     ▔▔▔▔▔▔▔▔▔▔▔                  
      1                                                                              86
```

### Pairing

We pair the generated primers using `best_pairs`, which takes a single vector of mixed forward and reverse primers, matches them based on the desired amplicon length and $T_m$ compatibility, and can sort the results using the `sortby` keyword.

```julia
julia> primer_pairs = best_pairs(primers; amplicon_len=40:50, sortby=:startpos)
97 PCR primer pairs for 5 seq. MSA:
     |==============================================================86|   Tm °C
[ 1]    >-------------46bp--------------<                              56.2±1.9
[ 2]    >--------------48bp---------------<                            57.4±0.7
[ 3]     >-------------45bp-------------<                              55.8±1.4
[ 4]     >--------------47bp--------------<                            57.0±0.3
[ 5]     >--------------47bp--------------<                            58.6±1.8
[ 6]      >------------44bp-------------<                              55.2±0.9
                              ... 85 more ...                             ...
[92]                       >--------------47bp--------------<          56.4±0.6
[93]                       >--------------48bp---------------<         55.6±1.4
[94]                       >--------------48bp---------------<         57.4±0.4
[95]                       >--------------49bp---------------<         55.9±1.1
[96]                       >--------------49bp---------------<         57.6±0.6
[97]                       >---------------50bp---------------<        56.3±0.7
```

### Renaming

Before exporting, you can rename the primer pairs using `reannotated`. Let's select two distinct pairs:

```julia
julia> a1 = reannotated(primer_pairs[6], "AMPL_1")
PCR primer pair for 5 seq. MSA, amplicon: 6:49 (44bp)
      >________________44bp_________________<                                     
|==============================================================================86|
Forward: GYAAYGAGCGKCAGACYGAC at 6:25 (AMPL_1)
Reverse: SGCTTCRGCSAGGTCKARTTC at 29:49 (AMPL_1)
Tm: 55.2±0.9 °C

julia> a2 = reannotated(primer_pairs[92], "AMPL_2")
PCR primer pair for 5 seq. MSA, amplicon: 30:76 (47bp)
                           >__________________47bp___________________<            
|==============================================================================86|
Forward: AAYTMGACCTSGCYGAAGCSCTG at 30:52 (AMPL_2)
Reverse: GAGCTTVGAATTGAGYTTGGC at 56:76 (AMPL_2)
Tm: 56.4±0.6 °C
```

### Designing Primers with Adapters

You can automatically append 5' adapters (e.g., for multiplexing or sequencing platforms) to your primers using `setAdapters!`. Once set, `construct_primers` will automatically append them, recalculate $\Delta G$ at the primer's mean $T_m$, and discard candidates where the adapter worsens $\Delta G$ by more than `max_dg_drop` (default: 1.0 kcal/mol).

```julia
julia> setAdapters!("TACGGGC" => "GAACGAT"); # Set some adapters for the toy example

julia> adp_primers = construct_primers(alignment)
Primer distribution for 5 seq MSA (L=86): 53 primers
27  ┤                                                                                  
    │                                                        ▁▂▂▂▂▂▂▂▂▂▂▂▂▁            
    │                                                    ▄▇▇▇██████████████▄           
    │                                                   ▆███████████████████▆▂         
    │     ▁▁▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▄▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁     ▄███████████████████████▅▂▂▂▂    
0   ┼ ╞════════════════════════════════════════════════════════════════════════════════
    │           ▔▔▔▔▔▔▔▔▔▔▔▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▔▔▔▔█████████████████████████▀▀▔       
    │                                               ▔████████████████████▀▔▔▔          
    │                                                ▀█████████████████▀▔              
    │                                                 ▀▀▀██████████████                
27  ┤                                                    ▔██████████▀▔▔                
      1                                                                              86

julia> setAdapters!(); # to reset or pass `adapter_pair=nothing` kwarg to `construct_primers`
```

Now we pair the generated primers with adapters:

```julia
julia> adp_pairs = best_pairs(adp_primers; amplicon_len=40:50, sortby=:startpos)
16 PCR primer pairs for 5 seq. MSA:
     |===================================================================86|   Tm °C
[ 1]      >--------------44bp--------------<                                55.2±0.9
[ 2]      >---------------46bp---------------<                              56.4±0.3
[ 3]       >-------------42bp--------------<                                54.2±0.2
[ 4]       >--------------44bp---------------<                              55.4±1.4
[ 5]                     >-------------42bp-------------<                   55.0±1.1
                                 ... 6 more ...                                ...
[12]                     >---------------46bp---------------<               55.6±1.7
[13]                     >---------------47bp---------------<               55.2±1.3
[14]                     >---------------48bp----------------<              54.8±0.9
[15]                     >---------------48bp----------------<              55.9±2.0
[16]                     >----------------49bp----------------<             55.6±1.7

julia> a1a = reannotated(adp_pairs[1], "AMPL_1_ADP")
PCR primer pair for 5 seq. MSA, amplicon: 6:49 (44bp)
      >__________________44bp__________________<                                       
|===================================================================================86|
Forward: [TACGGGC]GYAAYGAGCGKCAGACYGAC at 6:25 (AMPL_1_ADP)
Reverse: [GAACGAT]SGCTTCRGCSAGGTCKARTTC at 29:49 (AMPL_1_ADP)
Tm: 55.2±0.9 °C
```

### Specificity Check

To manually ensure a primer doesn't bind elsewhere, you can use `miniblast`. It performs a fast, probabilistic sliding-window search for high-probability matches of a query sequence within the target MSA, evaluating both forward and reverse complement orientations. Note that `miniblast` automatically queries the core primer sequence (ignoring the adapter).

```julia
julia> # Let's check the binding sites of the first pair's forward primer
julia> miniblast(alignment, a1a.first, 0.45)
5-element Vector{MiniBlastHit}:
 MiniBlastHit(6:25, :forward, 0.9)
 MiniBlastHit(62:81, :reverse, 0.5349999999999999)
 MiniBlastHit(18:37, :forward, 0.49000000000000005)
 MiniBlastHit(15:34, :reverse, 0.475)
 MiniBlastHit(10:29, :forward, 0.47000000000000003)

julia> # We can also search for a short conserved motif manually
julia> miniblast(alignment, "CAT", 0.7)
11-element Vector{MiniBlastHit}:
 MiniBlastHit(76:78, :forward, 1.0)
 MiniBlastHit(9:11, :reverse, 0.7333333333333334)
```

### Exporting

Finally, you can easily export your designed primers for ordering:

```julia
julia> export_evrogen(stdout, [a1, a1a, a2], scale=0.02)
AMPL_1_F_6; GYAAYGAGCGKCAGACYGAC; 0.02
AMPL_1_R_49; SGCTTCRGCSAGGTCKARTTC; 0.02
AMPL_1_ADP_F_6; TACGGGCGYAAYGAGCGKCAGACYGAC; 0.02
AMPL_1_ADP_R_49; GAACGATSGCTTCRGCSAGGTCKARTTC; 0.02
AMPL_2_F_30; AAYTMGACCTSGCYGAAGCSCTG; 0.02
AMPL_2_R_76; GAGCTTVGAATTGAGYTTGGC; 0.02

julia> export_evrogen("primers.txt", [a1, a1a, a2])
"primers.txt"
```
module DEPPA

"""
    Package DEPPA

Nucleic acid oligomers aligning and PCR primers construction

$(isnothing(get(ENV, "CI", nothing)) ? ("\nPackage local path: $(pathof(DEPPA))") : "") 
"""
DEPPA

export Oligos, Primers, Alignments


include("Oligos.jl")
include("Alignments.jl")
include("Primers.jl")

end # module
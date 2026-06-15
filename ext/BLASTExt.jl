module BLASTExt

using NCBIBlast
using FastaIO
using ProgressMeter

using DEPPA.Alignments
using DEPPA.Oligos  
using DEPPA.Primers

import DEPPA.Primers: construct_primers


# function construct_primers(
#     msa::AbstractMSA;
#     genome::String="",
#     blast_evalue::Real=1e-10,
#     blast_identity::Real=0.85,
#     blast_coverage::Real=0.7,
#     max_primer_variants::Int=100,
#     kwargs...
# )
#     if isempty(genome)
#         # Fall back to original function
#         return construct_primers(msa; kwargs...)
#     else
#         return construct_primers_genome_aware(
#             msa, genome; 
#             blast_evalue=blast_evalue,
#             blast_identity=blast_identity, 
#             blast_coverage=blast_coverage,
#             max_primer_variants=max_primer_variants,
#             kwargs...
#         )
#     end
# end


# function construct_primers_genome_aware(
#     msa::AbstractMSA,
#     genome_fasta::String;
#     blast_evalue::Real=1e-10,
#     blast_identity::Real=0.85,
#     blast_coverage::Real=0.7,
#     max_primer_variants::Int=100,
#     kwargs...
# )
#     candidates = construct_primers(
#         msa; 
#         max_oligo_variants=max_primer_variants,
#         kwargs...
#     )
#     @info "Generated $(length(candidates)) candidate primers"
    
#     screened_primers = _screen_primers_genome(
#         candidates, genome_fasta, 
#         blast_evalue, blast_identity, blast_coverage
#     )
#     @info "Screened to $(length(screened_primers)) genome-specific primers"
    
#     return screened_primers
# end


# function _screen_primers_genome(
#     primers::Vector{T},
#     genome_fasta::String,
#     evalue::Real,
#     identity::Real,
#     coverage::Real
# ) where {T <:AbstractPrimer}
#     screened_primers = T[]
#     @showprogress desc="Screening primers against genome..." for primer in primers
#         if _has_off_target_hits(primer, genome_fasta, evalue, identity, coverage)
#             continue
#         else
#             push!(screened_primers, primer)
#         end
#     end
#     return screened_primers
# end

# function _has_off_target_hits(
#     primer::AbstractPrimer,
#     genome_fasta::String,
#     evalue::Real,
#     identity::Real,
#     coverage::Real
# )::Bool
#     sample_sequences = _generate_primer_variants(primer, 10)
#     for seq in sample_sequences
#         if _find_genome_hits(seq, genome_fasta, evalue, identity, coverage)
#             return true
#         end
#     end
#     return false  # No off-target hits found
# end

# function _generate_primer_variants(
#     primer::AbstractPrimer, 
#     n_samples::Int
# )::Vector{String}
#     variants = String[]
    
#     if n_unique_oligos(primer) <= n_samples
#         for oligo in nondegens(primer.consensus)
#             push!(variants, String(oligo))
#         end
#     else
#         rng = Random.MersenneTwister(12345)
#         for _ in 1:n_samples
#             variant = rand(rng, primer)
#             push!(variants, String(variant))
#         end
#     end
    
#     return variants
# end


# function _find_genome_hits(
#     sequence::String,
#     genome_fasta::String,
#     evalue::Real,
#     identity::Real,
#     coverage::Real
# )::Bool
#     # Strategy 1: Use BLAST+ if available
#     if _has_blast_binaries()
#         return _blast_search(sequence, genome_fasta, evalue, identity, coverage)
    
#     # Strategy 2: Simple sequence matching (fallback)
#     else
#         return _simple_sequence_search(sequence, genome_fasta, identity, coverage)
#     end
# end

# function _blast_search(
#     sequence::String,
#     genome_fasta::String,
#     evalue::Real,
#     identity::Real,
#     coverage::Real
# )::Bool
#     # Create temporary query file
#     query_file = tempname() * ".fasta"
#     open(query_file, "w") do f
#         write(f, ">query\n$sequence\n")
#     end
    
#     try
#         # Run BLASTn
#         blast_cmd = `blastn -query $query_file -subject $genome_fasta -evalue $evalue -outfmt 6`
#         blast_output = read(blast_cmd, String)
        
#         # Parse BLAST output for significant hits
#         return _parse_blast_output(blast_output, identity, coverage)
#     finally
#         # Clean up
#         isfile(query_file) && rm(query_file)
#     end
# end
# function _simple_sequence_search(
#     sequence::String,
#     genome_fasta::String,
#     identity::Real,
#     coverage::Real
# )::Bool
#     seq_len = length(sequence)
#     min_match_len = floor(Int, seq_len * coverage)
    
#     # Read genome sequences
#     FastaReader(genome_fasta) do fr
#         for (desc, genome_seq) in fr
#             # Check for approximate matches
#             if _has_approximate_match(sequence, genome_seq, identity, min_match_len)
#                 return true
#             end
#         end
#     end
    
#     return false
# end
# function _has_approximate_match(
#     query::String,
#     subject::String,
#     identity::Real,
#     min_length::Int
# )::Bool
#     q_len = length(query)
#     s_len = length(subject)
    
#     # Sliding window approach
#     for start in 1:(s_len - q_len + 1)
#         substring = subject[start:(start + q_len - 1)]
#         similarity = _calculate_similarity(query, substring)
        
#         if similarity >= identity && length(substring) >= min_length
#             return true
#         end
#     end
    
#     return false
# end
# function _calculate_similarity(seq1::String, seq2::String)::Float64
#     if length(seq1) != length(seq2)
#         return 0.0
#     end
    
#     matches = sum(c1 == c2 for (c1, c2) in zip(seq1, seq2))
#     return matches / length(seq1)
# end
# function _has_blast_binaries()::Bool
#     try
#         run(`blastn -version`)
#         return true
#     catch
#         return false
#     end
# end
# function _parse_blast_output(
#     blast_output::String,
#     identity::Real,
#     coverage::Real
# )::Bool
#     isempty(blast_output) && return false
    
#     lines = split(blast_output, '\n', keepempty=false)
    
#     for line in lines
#         # BLAST tabular format: query_id, subject_id, identity, alignment_length, etc.
#         fields = split(line, '\t')
#         if length(fields) >= 3
#             identity_pct = parse(Float64, fields[3]) / 100.0
#             alignment_len = parse(Int, fields[4])
            
#             # Check if this hit meets our criteria
#             if identity_pct >= identity && alignment_len >= coverage * length(fields[1])
#                 return true
#             end
#         end
#     end
    
#     return false
# end



end  # module

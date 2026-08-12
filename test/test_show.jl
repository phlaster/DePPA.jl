using Test
using DePPA.Alignments
using DePPA.Oligos
using DePPA.Primers
using Random

@testset "Show Methods for MSA" begin
    seqs = ["ACGTACGT", "ACGTACGT", "ACGTACGT", "ACGTACGT", "ACGTACGT"]
    msa = MSA(seqs)
    
    # Test default (polymorf)
    output = sprint((io, x) -> show(io, MIME"text/plain"(), x), msa)
    @test occursin("MSA with 5 sequences of length 8", output)
    @test occursin("ACGTACGT", output)
    @test occursin("1", output)
    
    # Test bw style
    setMSAShowStyle!(:bw)
    output_bw = sprint((io, x) -> show(io, MIME"text/plain"(), x), msa)
    @test occursin("........", output_bw)
    setMSAShowStyle!(:polymorf)
    
    # Test allcolors style
    setMSAShowStyle!(:allcolors)
    output_allcolors = sprint((io, x) -> show(io, MIME"text/plain"(), x), msa)
    @test occursin("MSA with 5 sequences of length 8", output_allcolors)
    setMSAShowStyle!(:polymorf)
    
    # Test invalid style
    @test_throws ArgumentError setMSAShowStyle!(:invalid)
    
    # Test with descriptions (triggers has_desc logic)
    gapped_seqs = [GappedOligo("ACGT", "seq1_desc"), GappedOligo("TGCA", "seq2_desc")]
    msa_desc = MSA(gapped_seqs)
    output_desc = sprint((io, x) -> show(io, MIME"text/plain"(), x), msa_desc)
    @test occursin("seq1_desc", output_desc)
    @test occursin("seq2_desc", output_desc)
    
    # Test with gaps (triggers c != '-' for dot_chars)
    msa_gap = MSA(["A---ACGT", "ACGTACGT", "ACGTACGT", "ACGTACGT", "ACGTACGT"])
    output_gap = sprint((io, x) -> show(io, MIME"text/plain"(), x), msa_gap)
    @test occursin("-", output_gap)
    
    # Test MSAView
    msa_view = msa[1:2, 2:5]
    output_view = sprint((io, x) -> show(io, MIME"text/plain"(), x), msa_view)
    @test occursin("MSAView", output_view)
    @test occursin("CGTA", output_view)
    
    # Test empty MSA
    msa_empty = MSA(String[])
    output_empty = sprint((io, x) -> show(io, MIME"text/plain"(), x), msa_empty)
    @test occursin("Empty MSA", output_empty)
    
    # Test consensus type major
    setMSAconsensusShowType!(:major)
    output_major = sprint((io, x) -> show(io, MIME"text/plain"(), x), msa)
    @test occursin("ACGTACGT", output_major)
    setMSAconsensusShowType!(:degen)
    
    # Test invalid consensus type
    @test_throws ArgumentError setMSAconsensusShowType!(:invalid)
    
    Random.seed!(123)
    
    # Test width ellipsis (truncation by width)
    long_seq = join(rand("ACGT", 200))
    seqs_long = [long_seq for _ in 1:5]
    msa_long = MSA(seqs_long)
    
    # Force a small terminal width to trigger needs_width_ellipsis
    io = IOContext(IOBuffer(), :displaysize => (24, 40))
    show(io, MIME"text/plain"(), msa_long)
    output_long = String(take!(io.io))
    @test occursin("…", output_long)
    @test occursin("200", output_long) # Should print the last position number
    
    # Test width ellipsis in bw and allcolors styles
    setMSAShowStyle!(:bw)
    io = IOContext(IOBuffer(), :displaysize => (24, 40))
    show(io, MIME"text/plain"(), msa_long)
    output_long_bw = String(take!(io.io))
    @test occursin("…", output_long_bw)
    
    setMSAShowStyle!(:allcolors)
    io = IOContext(IOBuffer(), :displaysize => (24, 40))
    show(io, MIME"text/plain"(), msa_long)
    output_long_allcolors = String(take!(io.io))
    @test occursin("…", output_long_allcolors)
    setMSAShowStyle!(:polymorf)
    
    # Test number line markers (10, 100, 1000)
    seqs_1020 = [join(rand("ACGT", 1020)) for _ in 1:2]
    msa_1020 = MSA(seqs_1020)
    io = IOContext(IOBuffer(), :displaysize => (24, 1050))
    show(io, MIME"text/plain"(), msa_1020)
    output_1020 = String(take!(io.io))
    @test occursin("1020", output_1020)
    @test occursin('·', output_1020) # 10 marker
    @test occursin('⌢', output_1020) # 50 marker
    @test occursin(':', output_1020) # 100 marker
    @test occursin('∴', output_1020) # 1000 marker
    
    # Test height ellipsis (truncation by height with descriptions)
    seqs_many = [GappedOligo(join(rand("ACGT", 10)), "seq$i") for i in 1:100]
    msa_many = MSA(seqs_many)
    io = IOContext(IOBuffer(), :displaysize => (10, 80)) # Small terminal height
    show(io, MIME"text/plain"(), msa_many)
    output_many = String(take!(io.io))
    @test occursin("seq1", output_many)
    @test occursin("...", output_many) # gap_below_seqnames dots
end

@testset "Show Methods for Primers" begin
    msa = MSA(["ACGTACGT", "ACGTACGT"])
    tm = (mean=55.0, conf=(53.0, 57.0), min=53.0, max=57.0)
    
    p = Primer(msa, 1:8, true, Oligo("ACGTTGCA", "TestPrimer"), 3, tm, -5.0, 0.5, 0.0)
    output_p = sprint(show, p)
    @test occursin("Primer(\"ACGTTGCA\", len=8, pos=1:8, forward", output_p)
    @test occursin("degen=0, variants=1", output_p)
    @test occursin("Tm=55.0°C", output_p)
    
    output_p_mime = sprint((io, x) -> show(io, MIME"text/plain"(), x), p)
    @test occursin("Forward non-degenerate primer", output_p_mime)
    @test occursin("Sequence: ACGTTGCA", output_p_mime)
    @test occursin("Length: 8", output_p_mime)
    @test occursin("Positions: 1:8", output_p_mime)
    @test occursin("Melting temperature: 55.0°C", output_p_mime)
    
    p_fwd = Primer{DegenOligo}(msa, 1:4, true, Oligo("ACGT"), 3, tm, -5.0, 0.5, 0.0)
    p_rev = Primer{DegenOligo}(msa, 5:8, false, Oligo("TGCA"), 3, tm, -5.0, 0.5, 0.0)
    pair = p_fwd => p_rev
    
    output_pair = sprint((io, x) -> show(io, MIME"text/plain"(), x), pair)
    @test occursin("PCR primer pair", output_pair)
    @test occursin("amplicon: 1:8 (8bp)", output_pair)
    @test occursin("Forward: ACGT at 1:4", output_pair)
    @test occursin("Reverse: TGCA at 5:8", output_pair)
    @test occursin("Tm: 55.0±0.0 °C", output_pair)
    
    p_rev_overlap = Primer{DegenOligo}(msa, 3:6, false, Oligo("GTAC"), 3, tm, -5.0, 0.5, 0.0)
    pair_overlap = p_fwd => p_rev_overlap
    output_overlap = sprint((io, x) -> show(io, MIME"text/plain"(), x), pair_overlap)
    @test occursin("OVERLAPPING", output_overlap)
end
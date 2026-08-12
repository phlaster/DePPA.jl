using Test
using DePPA.Alignments
using DePPA.Oligos
using DePPA.Primers

@testset "Show Methods for MSA" begin
    seqs = ["ACGTACGT", "ACGTACGT", "ACGTACGT", "ACGTACGT", "ACGTACGT"]
    msa = MSA(seqs)
    
    output = sprint((io, x) -> show(io, MIME"text/plain"(), x), msa)
    @test occursin("MSA with 5 sequences of length 8", output)
    @test occursin("ACGTACGT", output)
    @test occursin("1", output)
    
    setMSAShowStyle!(:bw)
    output_bw = sprint((io, x) -> show(io, MIME"text/plain"(), x), msa)
    @test occursin("........", output_bw)
    setMSAShowStyle!(:polymorf)
    
    msa_gap = MSA(["A---ACGT", "ACGTACGT", "ACGTACGT", "ACGTACGT", "ACGTACGT"])
    output_gap = sprint((io, x) -> show(io, MIME"text/plain"(), x), msa_gap)
    @test occursin("-", output_gap)
    
    msa_view = msa[1:2, 2:5]
    output_view = sprint((io, x) -> show(io, MIME"text/plain"(), x), msa_view)
    @test occursin("MSAView", output_view)
    @test occursin("CGTA", output_view)
    
    msa_empty = MSA(String[])
    output_empty = sprint((io, x) -> show(io, MIME"text/plain"(), x), msa_empty)
    @test occursin("Empty MSA", output_empty)
    
    setMSAconsensusShowType!(:major)
    @test DePPA.Alignments.CONSENSUS_SHOW_TYPE[] == :major
    setMSAconsensusShowType!(:degen)
    @test DePPA.Alignments.CONSENSUS_SHOW_TYPE[] == :degen
    @test_throws ArgumentError setMSAconsensusShowType!(:invalid)
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
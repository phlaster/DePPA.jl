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
    tm = (mean = 55.0, conf = (53.0, 57.0), min = 53.0, max = 57.0)

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

    p_rev_overlap = Primer{DegenOligo}(
        msa,
        3:6,
        false,
        Oligo("GTAC"),
        3,
        tm,
        -5.0,
        0.5,
        0.0,
    )
    pair_overlap = p_fwd => p_rev_overlap
    output_overlap = sprint((io, x) -> show(io, MIME"text/plain"(), x), pair_overlap)
    @test occursin("OVERLAPPING", output_overlap)
end

@testset "Show Vector of Primers" begin
    msa = MSA(["ACGTACGTACGTACGT", "ACGTACGTACGTACGT"])
    tm = (mean = 55.0, conf = (53.0, 57.0), min = 53.0, max = 57.0)

    p1 = Primer{DegenOligo}(msa, 1:4, true, Oligo("ACGT"), 3, tm, -5.0, 0.5, 0.0)
    p2 = Primer{DegenOligo}(msa, 5:8, true, Oligo("ACGT"), 3, tm, -5.0, 0.5, 0.0)
    r1 = Primer{DegenOligo}(msa, 9:12, false, Oligo("ACGT"), 3, tm, -5.0, 0.5, 0.0)
    r2 = Primer{DegenOligo}(msa, 13:16, false, Oligo("ACGT"), 3, tm, -5.0, 0.5, 0.0)

    # 1. Empty vector
    out_empty = sprint((io, x) -> show(io, MIME"text/plain"(), x), Primer{DegenOligo}[])
    @test out_empty == "Empty primer vector"

    # 2. Normal vector (forces specific terminal size for consistent testing)
    primers = [p1, p2, r1, r2]
    io = IOContext(IOBuffer(), :displaysize => (24, 80))
    show(io, MIME"text/plain"(), primers)
    out = String(take!(io.io))
    
    @test occursin("Primer distribution for 2 seq MSA (L=16): 4 primers", out)
    @test occursin("█", out)      # Check for histogram blocks
    @test occursin("╞", out)      # Check for axis start
    @test occursin("═", out)      # Check for axis line
    @test occursin("16", out)     # Check for end position label

    # 3. Different MSAs fallback
    msa2 = MSA(["TTTTTTTTTTTTTTTT", "TTTTTTTTTTTTTTTT"])
    p_diff = Primer{DegenOligo}(msa2, 1:4, true, Oligo("TTTT"), 3, tm, -5.0, 0.5, 0.0)
    out_diff = sprint((io, x) -> show(io, MIME"text/plain"(), x), [p1, p_diff])
    @test occursin("Vector of 2 primers (different MSAs):", out_diff)
    @test occursin("Primer(\"", out_diff) # Falls back to single-line show
end

@testset "Show Vector of Primer Pairs" begin
    msa = MSA(["ACGTACGTACGTACGT", "ACGTACGTACGTACGT"])
    tm1 = (mean = 55.0, conf = (53.0, 57.0), min = 53.0, max = 57.0)
    tm2 = (mean = 56.5, conf = (54.0, 58.0), min = 54.0, max = 58.0)

    p1 = Primer{DegenOligo}(msa, 1:4, true, Oligo("ACGT"), 3, tm1, -5.0, 0.5, 0.0)
    p2 = Primer{DegenOligo}(msa, 5:8, true, Oligo("ACGT"), 3, tm2, -5.0, 0.5, 0.0)
    r1 = Primer{DegenOligo}(msa, 9:12, false, Oligo("ACGT"), 3, tm1, -5.0, 0.5, 0.0)
    r2 = Primer{DegenOligo}(msa, 13:16, false, Oligo("ACGT"), 3, tm2, -5.0, 0.5, 0.0)

    # 1. Empty vector
    pair_type = Pair{Primer{DegenOligo}, Primer{DegenOligo}}
    out_empty = sprint((io, x) -> show(io, MIME"text/plain"(), x), pair_type[])
    @test out_empty == "Empty primer pair vector"

    # 2. Normal pairs
    pairs = [p1 => r1, p2 => r2]
    io = IOContext(IOBuffer(), :displaysize => (24, 80))
    show(io, MIME"text/plain"(), pairs)
    out = String(take!(io.io))
    
    @test occursin("2 PCR primer pairs for 2 seq. MSA:", out)
    @test occursin("Tm °C", out)
    # p1(1:4) => r1(9:12) -> amplicon 1:12 (12bp)
    @test occursin("12bp", out)
    # p2(5:8) => r2(13:16) -> amplicon 5:16 (12bp)
    @test occursin("55.0±0.0", out)
    @test occursin("56.5±0.0", out)

    # 3. Truncated pairs (height limit triggers "... N more ...")
    many_pairs = [p1 => r1 for _ in 1:20]
    io_small = IOContext(IOBuffer(), :displaysize => (10, 80))
    show(io_small, MIME"text/plain"(), many_pairs)
    out_trunc = String(take!(io_small.io))
    
    @test occursin("20 PCR primer pairs for 2 seq. MSA:", out_trunc)
    @test occursin("more ...", out_trunc)

    # 4. Different MSAs fallback
    msa2 = MSA(["TTTTTTTTTTTTTTTT", "TTTTTTTTTTTTTTTT"])
    p_diff = Primer{DegenOligo}(msa2, 1:4, true, Oligo("TTTT"), 3, tm1, -5.0, 0.5, 0.0)
    out_diff = sprint((io, x) -> show(io, MIME"text/plain"(), x), [p1 => r1, p_diff => r1])
    @test occursin("Vector of 2 primer pairs:", out_diff)
end

@testset "Show Edge Cases & Truncation" begin
    import DePPA.Primers: _truncate_seq, _primer_seq_str

    msa = MSA(["ACGTACGT", "ACGTACGT"])
    tm = (mean = 55.0, conf = (53.0, 57.0), min = 53.0, max = 57.0)

    # --- _truncate_seq ---
    @test _truncate_seq("ACGT", 10) == "ACGT"   # length <= w
    @test _truncate_seq("ACGT", 0) == ""        # w <= 0
    @test _truncate_seq("ACGT", -1) == ""       # w <= 0
    @test _truncate_seq("ACGT", 3) == "AC…"     # truncated

    # --- _primer_seq_str ---
    p_no_adp = Primer{DegenOligo}(msa, 1:4, true, Oligo("ACGT"), 3, tm, -5.0, 0.5, 0.0, nothing)
    p_adp = Primer{DegenOligo}(msa, 1:4, true, Oligo("ACGT"), 3, tm, -5.0, 0.5, 0.0, Oligo("GGGGG", "adp"))
    p_long_adp = Primer{DegenOligo}(msa, 1:8, true, Oligo("ACGTACGT"), 3, tm, -5.0, 0.5, 0.0, Oligo("GGGGG", "adp"))

    # 1. No adapter
    @test _primer_seq_str(p_no_adp, 10) == "ACGT"
    @test _primer_seq_str(p_no_adp, 3) == "AC…"

    # 2. Adapter, full fits
    @test _primer_seq_str(p_adp, 20) == "[GGGGG]ACGT"

    # 3. Adapter, full doesn't fit, but max_width >= length(main) + 4
    # main = "ACGT" (4), adapter = "GGGGG" (5). max_width = 8.
    # keep = 8 - 4 - 3 = 1. -> "[G…]ACGT"
    @test _primer_seq_str(p_adp, 8) == "[G…]ACGT"

    # 4. Adapter, full doesn't fit, max_width < length(main) + 4
    # main = "ACGTACGT" (8), adapter = "GGGGG" (5). max_width = 6.
    # keep_main = max(6 - 4, 1) = 2. _truncate_seq("ACGTACGT", 2) = "A…"
    # Result: "[G…]A…"
    @test _primer_seq_str(p_long_adp, 6) == "[G…]A…"

    # --- Base.show for single Primer ---
    # L == 1 case
    msa_L1 = MSA(["A", "A"])
    p_L1 = Primer(msa_L1, 1:1, true, Oligo("A", "L1"), 3, tm, -5.0, 0.5, 0.0)
    out_L1 = sprint((io, x) -> show(io, MIME"text/plain"(), x), p_L1)
    @test occursin("Forward non-degenerate primer", out_L1)
    @test occursin("Sequence: A", out_L1)

    # Reverse primer (else branch in show layout)
    p_rev = Primer{DegenOligo}(msa, 5:8, false, Oligo("ACGT", "Rev"), 3, tm, -5.0, 0.5, 0.0)
    out_rev = sprint((io, x) -> show(io, MIME"text/plain"(), x), p_rev)
    @test occursin("Reverse non-degenerate primer", out_rev)
    # Check that bar is printed BEFORE the label (which starts with '<')
    bar_idx = findfirst("|=", out_rev)
    label_idx = findfirst('<', out_rev)
    @test !isnothing(bar_idx) && !isnothing(label_idx) && bar_idx[1] < label_idx

    # --- Base.show for Pair ---
    p_fwd = Primer{DegenOligo}(msa, 1:4, true, Oligo("ACGT", "Fwd"), 3, tm, -5.0, 0.5, 0.0)

    # 1. Wrong direction fallback (!fwd.is_forward || rev.is_forward)
    out_wrong_dir = sprint((io, x) -> show(io, MIME"text/plain"(), x), p_rev => p_fwd)
    @test !occursin("PCR primer pair", out_wrong_dir)
    @test occursin("Primer(\"", out_wrong_dir) # Falls back to default Pair show

    # 2. Different MSA fallback (rev.msa !== msa)
    msa2 = MSA(["TTTTTTTT", "TTTTTTTT"])
    p_msa2 = Primer{DegenOligo}(msa2, 1:4, true, Oligo("TTTT", "Fwd2"), 3, tm, -5.0, 0.5, 0.0)
    p_rev_msa1 = Primer{DegenOligo}(msa, 5:8, false, Oligo("ACGT", "Rev1"), 3, tm, -5.0, 0.5, 0.0)
    out_diff_msa = sprint((io, x) -> show(io, MIME"text/plain"(), x), p_msa2 => p_rev_msa1)
    @test !occursin("PCR primer pair", out_diff_msa)
    @test occursin("Primer(\"", out_diff_msa)

    # 3. Exception fallback (catch e)
    # We create a mock MSA that throws an error during nseqs() calculation
    struct CrashingMSA <: AbstractMSA end
    DePPA.Alignments.nseqs(::CrashingMSA) = error("Crash MSA!")
    DePPA.Alignments.length(::CrashingMSA) = 10
    DePPA.Alignments.root(::CrashingMSA) = CrashingMSA()

    p_crash = Primer{DegenOligo}(CrashingMSA(), 1:4, true, Oligo("ACGT"), 3, tm, -5.0, 0.5, 0.0)
    p_ok_rev = Primer{DegenOligo}(CrashingMSA(), 5:8, false, Oligo("ACGT"), 3, tm, -5.0, 0.5, 0.0)
    
    # The show method should catch the error and fall back to default Base.show
    out_catch = sprint((io, x) -> show(io, MIME"text/plain"(), x), p_crash => p_ok_rev)
    @test !occursin("PCR primer pair", out_catch)
    @test occursin("Primer(\"", out_catch) # Falls back to default show
end
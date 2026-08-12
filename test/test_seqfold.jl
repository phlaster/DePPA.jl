using Test
using SeqFold
using DePPA
using DePPA.Oligos
using DePPA.Alignments
using DePPA.Primers
using Random

Random.seed!(42)

@testset "SeqFoldExt Tests" begin
    
    @testset "Complement and Reverse Complement" begin
        # Oligo
        o = Oligo("ACGTA", "test_o")
        rc_o = SeqFold.revcomp(o)
        @test rc_o isa Oligo
        @test String(rc_o) == "TACGT"
        @test occursin("Reverse complement of test_o", description(rc_o))
        
        c_o = SeqFold.complement(o)
        @test c_o isa Oligo
        @test String(c_o) == "TGCAT"
        @test occursin("Complement of test_o", description(c_o))

        # DegenOligo
        d = DegenOligo("ACNG", "test_d")
        rc_d = SeqFold.revcomp(d)
        @test rc_d isa DegenOligo
        @test String(rc_d) == "CNGT"
        
        c_d = SeqFold.complement(d)
        @test c_d isa DegenOligo
        @test String(c_d) == "TGNC"

        # GappedOligo
        g = GappedOligo("A-C-G", "test_g")
        rc_g = SeqFold.revcomp(g)
        @test rc_g isa GappedOligo
        @test String(rc_g) == "C-G-T"
        
        c_g = SeqFold.complement(g)
        @test c_g isa GappedOligo
        @test String(c_g) == "T-G-C"
    end

    @testset "GC Content" begin
        @test isnan(SeqFold.gc_content(Oligo("")))
        @test SeqFold.gc_content(Oligo("ACGT")) == 0.5
        @test SeqFold.gc_content(Oligo("GGCC")) == 1.0
        @test SeqFold.gc_content(Oligo("AATC")) == 0.25
        
        # Degenerate
        @test SeqFold.gc_content(DegenOligo("SSSS")) == 1.0  # S = C/G
        @test SeqFold.gc_content(DegenOligo("WWWW")) == 0.0  # W = A/T
        @test SeqFold.gc_content(DegenOligo("SW")) == 0.5
        
        # Gapped delegates to parent
        @test SeqFold.gc_content(GappedOligo("A-C-G")) == SeqFold.gc_content(Oligo("ACG"))
        @test isnan(SeqFold.gc_content(GappedOligo("---")))
    end

    @testset "Free Energy (dg)" begin
        @test isnan(SeqFold.dg(Oligo("")))
        @test_throws ErrorException SeqFold.dg(GappedOligo("A-C"))
        @test_throws ErrorException SeqFold.dg(DegenOligo("ACNG"), mode=:invalid_mode)

        # Non-degenerate should match SeqFold.dg(String, ...)
        d1 = SeqFold.dg(Oligo("ATGGATTTAGATAGAT"))
        d2 = SeqFold.dg("ATGGATTTAGATAGAT")
        @test d1 == d2

        # Degenerate: Average vs Worstcase
        deg = DegenOligo("ACNG", "test")
        avg_dg = SeqFold.dg(deg; mode=:average, max_samples=10)
        worst_dg = SeqFold.dg(deg; mode=:worstcase, max_samples=10)
        @test avg_dg isa Float64
        @test worst_dg isa Float64
        @test worst_dg <= avg_dg  # worstcase should be more negative or equal

        # Sampling path
        deg_large = DegenOligo("NNNNNNNN", "test_large")
        # 4^8 = 65536 variants, max_samples=100 forces sampling
        samp_dg = SeqFold.dg(deg_large; max_samples=100, mode=:average)
        @test samp_dg isa Float64
        @test !isnan(samp_dg)
        
        samp_dg_worst = SeqFold.dg(deg_large; max_samples=100, mode=:worstcase)
        @test samp_dg_worst isa Float64
        @test !isnan(samp_dg_worst)
    end

    @testset "Melting Temperature (tm)" begin
        @test_throws ErrorException SeqFold.tm(GappedOligo("A-C"))
        @test_throws ArgumentError SeqFold.tm(Oligo("ACGT"), Oligo("ACGT"), conf_int=0.0)
        @test_throws ArgumentError SeqFold.tm(Oligo("ACGT"), Oligo("ACGT"), conf_int=1.5)
        @test_throws DimensionMismatch SeqFold.tm(Oligo("ACGT"), Oligo("ACG"))

        # Non-degenerate pair
        tm1 = SeqFold.tm(Oligo("ACGTACGT"), Oligo("ACGTACGT"))
        tm2 = SeqFold.tm("ACGTACGT", "ACGTACGT")
        @test tm1.mean == tm2
        @test tm1.min == tm1.max == tm1.mean

        # Degenerate pair
        tm_deg = SeqFold.tm(DegenOligo("ACGTACGN"), DegenOligo("ACGTACGN"), max_samples=10)
        @test tm_deg.mean isa Float64
        @test tm_deg.min <= tm_deg.mean <= tm_deg.max
        @test length(tm_deg.conf) == 2

        # Self-complementary
        tm_self = SeqFold.tm(Oligo("ACGTACGT"))
        @test tm_self.mean >= 0.0

        # Sampling path
        deg_large = DegenOligo("NNNNNNNN", "test_large")
        tm_samp = SeqFold.tm(deg_large, deg_large, max_samples=10)
        @test tm_samp.mean isa Float64
    end

    @testset "Cache Methods" begin
        o = Oligo("ACGTACGT")
        d = DegenOligo("ACGTACGT")
        
        @test_throws ErrorException SeqFold.dg_cache(GappedOligo("A-C-G"), temp=37.0)
        @test_throws ErrorException SeqFold.tm_cache(GappedOligo("A-C-G"), GappedOligo("A-C-G"))
        
        @test size(SeqFold.dg_cache(o, temp=37.0)) == (8, 8)
        @test size(SeqFold.dg_cache(d, temp=37.0)) == (8, 8)
        
        @test size(SeqFold.tm_cache(o, o)) == (8, 8)
        @test size(SeqFold.tm_cache(d, d)) == (8, 8)
        @test size(SeqFold.tm_cache(o)) == (8, 8) # Self-complement cache
        
        @test size(SeqFold.gc_cache(o)) == (8, 8)
        @test size(SeqFold.gc_cache(d)) == (8, 8)
    end

    @testset "Dot Bracket" begin
        seq_str = "GGGAGGTCGTTACATCTGGGTAACACCGGTACTGATCCGGTGACCTCCC"
        o = Oligo(seq_str)
        structs = SeqFold.fold(seq_str)
        
        db_str = SeqFold.dot_bracket(seq_str, structs)
        db_o = SeqFold.dot_bracket(o, structs)
        
        @test db_o == db_str
        @test length(db_o) == length(seq_str)
    end

    @testset "Unfolded Proportion" begin
        @test isnan(Oligos.unfolded_proportion(Oligo(""), temp=37.0, max_samples=10))
        @test_throws ErrorException Oligos.unfolded_proportion(GappedOligo("A-C"), temp=37.0, max_samples=10)

        # Single variant
        p1 = Oligos.unfolded_proportion(Oligo("ACGTACGT"), temp=37.0, max_samples=10)
        @test 0.0 <= p1 <= 1.0

        # Degenerate (enumeration)
        p2 = Oligos.unfolded_proportion(DegenOligo("ACGNACGN"), temp=37.0, max_samples=100)
        @test 0.0 <= p2 <= 1.0

        # Degenerate (sampling)
        p3 = Oligos.unfolded_proportion(DegenOligo("NNNNNNNN"), temp=37.0, max_samples=10)
        @test 0.0 <= p3 <= 1.0
    end

    @testset "Primer Integration & Hooks" begin
        msa = MSA(["ACGTACGT", "ACGTACGT"])
        o = Oligo("ACGTACGT", "hook_test")
        
        # Test internal hooks
        @test Primers._ext_revcomp(o) == SeqFold.revcomp(o)
        @test Primers._ext_tm(o; max_samples=100, conf_int=0.8, conditions=:pcr) == SeqFold.tm(o; max_samples=100, conf_int=0.8, conditions=:pcr)
        @test Primers._ext_dg(o; max_samples=100, temp=37.0) == SeqFold.dg(o; max_samples=100, temp=37.0)
        @test Primers._ext_gc_content(o) == SeqFold.gc_content(o)
        
        # Test Primer field accessors
        tm_stats = (mean=55.0, conf=(53.0, 57.0), min=53.0, max=57.0)
        p = Primer(msa, 1:8, true, Oligo("ACGTTGCA", "TestPrimer"), 3, tm_stats, -5.0, 0.5, 0.0)
        
        @test SeqFold.tm(p) === tm_stats
        @test SeqFold.dg(p) == -5.0
        @test SeqFold.gc_content(p) == 0.5
    end
end
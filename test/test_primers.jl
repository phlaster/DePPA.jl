using Test
using DePPA.Primers
using DePPA.Alignments
using DePPA.Oligos
using SeqFold

@testset "Primers Tests" begin

    @testset "Primer Constructor & Extensions" begin
        msa = MSA(["ACGTACGT", "ACGTACGT"])

        # Argument validation passed down to consensus_degen
        @test_throws ArgumentError Primer(msa, 1:4; slack = -0.1)
        @test_throws ArgumentError Primer(msa, 1:4; slack = 1.0)
    end

    @testset "Primer Struct & Base Methods" begin
        msa = MSA(["ACGTACGT", "ACGTACGT"])
        tm = (mean = 55.0, conf = (53.0, 57.0), min = 53.0, max = 57.0)

        # Test with Oligo consensus
        oligo_cons = Oligo("ACGTTGCA", "TestOligoPrimer")
        p1 = Primer(msa, 1:8, true, oligo_cons, 3, tm, -5.0, 0.5, 0.0)

        @test String(p1) == "ACGTTGCA"
        @test length(p1) == 8
        @test !isempty(p1)
        @test collect(p1) == ['A', 'C', 'G', 'T', 'T', 'G', 'C', 'A']
        @test p1[1] == 'A'
        @test String(p1[2:4]) == "CGT"
        @test convert(DegenOligo, p1) isa DegenOligo
        @test String(convert(DegenOligo, p1)) == String(oligo_cons)
        @test n_unique_oligos(p1) == BigInt(1)
        @test iszero(n_deg_pos(p1))
        @test description(p1) == "TestOligoPrimer"
        @test !hasgaps(p1)
        @test nondegens(p1)[1] == oligo_cons
        @test oligo_range(p1) == 1:8

        # Test with DegenOligo consensus
        degen_cons = DegenOligo("ACGTSWCA", "TestDegenPrimer")
        p2 = Primer(msa, 1:8, true, degen_cons, 3, tm, -5.0, 0.5, 0.0)

        @test String(p2) == "ACGTSWCA"
        @test length(p2) == 8
        @test n_unique_oligos(p2) == BigInt(4)
        @test n_deg_pos(p2) == 2
        @test description(p2) == "TestDegenPrimer"
        @test !hasgaps(p2)
        @test oligo_range(p2) == 1:8

        # Test empty primer
        empty_cons = Oligo("", "EmptyPrimer")
        p_empty = Primer{DegenOligo}(msa, 1:0, true, empty_cons, 3, tm, -5.0, 0.5, 0.0)
        @test isempty(p_empty)
        @test iszero(length(p_empty))
    end

    @testset "construct_primers" begin
        msa = MSA(["ACGTACGT", "ACGTACGT"])

        # Argument validation
        @test_throws ArgumentError construct_primers(msa; slack = -0.1)
        @test_throws ArgumentError construct_primers(msa; min_msadepth = -0.1)
        @test_throws ArgumentError construct_primers(msa; length_range = -3:2)
        @test_throws ArgumentError construct_primers(msa; tail_length = -3)
        @test_throws ArgumentError construct_primers(msa; gc_range = -3:10)
        @test_throws ArgumentError construct_primers(msa; tm_range = -30:200)
        @test_throws ArgumentError construct_primers(msa; offtarget_reject_threshold = 1.2)

        seq1 = "GATCTGTAATGAGCGGCAGACCGACCGCGAATTAGACCTCGCCGAAGCCCTGGCCGCCAAGCTCAATTCGAAGCTCATTCAC----"
        seq2 = "CATTTGCAACGAGCGTCAGACCGACCGCGAACTCGACCTGGCCGAAGCGCTGGCTGCCAAACTCAATTCTAAGCTCATC-------"
        seq3 = "CATTTGTAACGAGCGTCAGACCGACCGTGAACTCGACCTCGCCGA-------GCTGCCAAACTCAATTCCAAGCTCATCCACTT--"
        seq4 = "---CTGTAACGAGCGGCAGACTGACCGAGAATTAGACCTCGCTGAAGCGCTGGCCGCCAAGCTCAATTCGAAGCTCATTCACTTTG"
        seq5 = "----TGTAACGAGCGGCAGACTGACCGAGAATTAGACCTCGCTGAAGCGCTGGCCGCCAAGCTCAATTCGAAGCTCATTCACTTAG"

        msa_1 = MSA([seq1, seq2, seq3, seq4, seq5])

        primers_relaxed = construct_primers(msa_1; offtarget_reject_threshold = 0.9)
        @test !isempty(primers_relaxed)

        primers_strict = construct_primers(msa_1; offtarget_reject_threshold = 0.1)
        @test isempty(primers_strict)
    end

    @testset "best_pairs" begin
        msa1 = MSA(["ACGTACGT", "ACGTACGT"])
        msa2 = MSA(["TGCATGCA", "TGCATGCA"])

        tm1 = (mean = 55.0, conf = (53.0, 57.0), min = 53.0, max = 57.0)
        tm2 = (mean = 56.0, conf = (54.0, 58.0), min = 54.0, max = 58.0)
        tm3 = (mean = 60.0, conf = (58.0, 62.0), min = 58.0, max = 62.0)

        # Forwards
        f1 = Primer{DegenOligo}(msa1, 1:4, true, Oligo("ACGT"), 3, tm1, -5.0, 0.5, 0.0)
        f2 = Primer{DegenOligo}(msa1, 1:3, true, Oligo("ACG"), 3, tm1, -5.0, 0.5, 0.0)

        # Reverses
        r1 = Primer{DegenOligo}(msa1, 5:8, false, Oligo("ACGT"), 3, tm1, -5.0, 0.5, 0.0)
        r2 = Primer{DegenOligo}(msa1, 3:6, false, Oligo("GTAC"), 3, tm1, -5.0, 0.5, 0.0)
        r3 = Primer{DegenOligo}(msa1, 5:8, false, Oligo("ACGT"), 3, tm3, -5.0, 0.5, 0.0)
        r4 = Primer{DegenOligo}(msa2, 5:8, false, Oligo("ACGT"), 3, tm1, -5.0, 0.5, 0.0)
        r5 = Primer{DegenOligo}(msa1, 5:8, false, Oligo("ACGT"), 3, tm2, -5.0, 0.5, 0.0)

        # Valid pair
        pairs = best_pairs([f1, r1])
        @test length(pairs) == 1
        @test pairs[1].first === f1
        @test pairs[1].second === r1

        # Overlapping primers
        @test isempty(best_pairs([f1, r2]))

        # Tm difference too large (default max_tm_diff=4.0)
        @test isempty(best_pairs([f1, r3]))

        # Tm difference acceptable
        @test length(best_pairs([f1, r3]; max_tm_diff = 6.0)) == 1

        # Amplicon length filter (f1=1:4, r1=5:8 -> amplicon = 8 - 1 + 1 = 8)
        @test isempty(best_pairs([f1, r1]; amplicon_len = 1:7))
        @test length(best_pairs([f1, r1]; amplicon_len = 8:10)) == 1

        # Invalid MSA
        @test_throws ArgumentError best_pairs([f1, r4])

        # No forwards / no reverses -> empty result
        @test isempty(best_pairs([r1]))   # only reverses
        @test isempty(best_pairs([f1]))   # only forwards

        # Empty list
        @test isempty(best_pairs(Primer{DegenOligo}[]))

        # Sort by Tm diff
        pairs_sorted = best_pairs([f1, r1, r5]; sortby = :tm_diff)
        @test length(pairs_sorted) == 2
        @test pairs_sorted[1].second === r1 # diff 0.0
        @test pairs_sorted[2].second === r5 # diff 1.0
    end

    @testset "Non-specific binding check" begin
        # MSA with exact off-target duplications
        msa_exact = MSA(["ACGTACGTGGGGACGTACGT", "ACGTACGTGGGGACGTACGT"])
        # Target 1:8 (ACGTACGT). Off-target at 13:20 (ACGTACGT).
        primer_exact = "ACGTACGT"

        @test DePPA.Primers._has_nonspecific_match(
            primer_exact,
            msa_exact,
            1:8;
            min_identity = 0.8,
        )
        @test DePPA.Primers._has_nonspecific_match(
            primer_exact,
            msa_exact,
            1:8;
            min_identity = 1.0,
        )

        # MSA with partial off-target (75% identity: 6/8 matches)
        # ACGTACGT vs ACGGAGGT -> 6/8 match = 0.75
        msa_partial = MSA(["ACGTACGTGGGGACGGAGGT", "ACGTACGTGGGGACGGAGGT"])
        @test !DePPA.Primers._has_nonspecific_match(
            primer_exact,
            msa_partial,
            1:8;
            min_identity = 0.8,
        )
        @test DePPA.Primers._has_nonspecific_match(
            primer_exact,
            msa_partial,
            1:8;
            min_identity = 0.7,
        )

        # MSA with reverse complement off-target
        # Target: AAAAAAAA (1:8)
        # RC of AAAAAAAA is TTTTTTTT. Located at 13:20.
        msa_rc = MSA(["AAAAAAAATTTTCCCCTTTTTTTT", "AAAAAAAATTTTCCCCTTTTTTTT"])
        primer_a = "AAAAAAAA"
        @test DePPA.Primers._has_nonspecific_match(
            primer_a,
            msa_rc,
            1:8;
            min_identity = 0.8,
        )

        # Clean MSA (no off-targets)
        msa_clean = MSA(["AAAAAAAATTTTCCCCGGGG", "AAAAAAAATTTTCCCCGGGG"])
        @test !DePPA.Primers._has_nonspecific_match(
            primer_a,
            msa_clean,
            1:8;
            min_identity = 0.8,
        )

        # Make sure it doesn't flag the target itself as off-target
        @test !DePPA.Primers._has_nonspecific_match(
            primer_a,
            msa_clean,
            1:8;
            min_identity = 1.0,
        )
    end
end

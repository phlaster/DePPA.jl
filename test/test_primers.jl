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

    @testset "Adapter Management" begin
        # Ensure we start with a clean state
        setAdapters!()
        @test isnothing(getAdapters())

        # Set via strings
        setAdapters!("AAA" => "TTT")
        @test getAdapters() == (Oligo("AAA", "Custom Fwd Adapter") => Oligo("TTT", "Custom Rev Adapter"))

        # Set via Oligo
        setAdapters!(Oligo("CCC") => Oligo("GGG"))
        @test getAdapters() == (Oligo("CCC") => Oligo("GGG"))

        # Reset
        setAdapters!()
        @test isnothing(getAdapters())
    end

    @testset "Primer Constructor & Hooks" begin
        msa = MSA(["ACGTACGT", "ACGTACGT"])

        # 1. Basic constructor with default kwargs
        p1 = Primer(msa, 1:8)
        @test p1 isa Primer{DegenOligo}
        @test isnothing(p1.adapter)
        @test String(p1) == "ACGTACGT"

        # 2. Constructor with explicit adapter_pair
        p2 = Primer(msa, 1:8; adapter_pair = Oligo("GGG") => Oligo("CCC"))
        @test p2.adapter == Oligo("GGG")
        @test String(p2) == "GGGACGTACGT" # Check that adapter was added to the 5' end

        # 3. Constructor with custom description (descr)
        p3 = Primer(msa, 1:8; descr = "MyCustomPrimer")
        @test description(p3) == "MyCustomPrimer"

        # 4. Error when gaps are present in consensus
        # If both sequences have gaps in the interval, sum(p) will be 0.0 < 0.5 -> '-'
        msa_gap = MSA(["A---ACGT", "T---TGCA"])
        @test_throws ArgumentError Primer(msa_gap, 1:4)
    end

    @testset "Reannotation" begin
        msa = MSA(["ACGTACGT", "ACGTACGT"])
        tm = (mean = 55.0, conf = (53.0, 57.0), min = 53.0, max = 57.0)

        p1 = Primer{DegenOligo}(msa, 1:4, true, Oligo("ACGT", "OldFwd"), 3, tm, -5.0, 0.5, 0.0)
        p2 = Primer{DegenOligo}(msa, 5:8, false, Oligo("TGCA", "OldRev"), 3, tm, -5.0, 0.5, 0.0)

        # Reannotation of a single primer (immutability)
        p1_new = reannotated(p1, "NewFwd")
        @test p1_new.pos == p1.pos
        @test p1_new.tm == p1.tm
        @test description(p1_new) == "NewFwd"
        @test description(p1) == "OldFwd" # Original unchanged

        # Reannotation of a pair
        pair = p1 => p2
        pair_new = reannotated(pair, "PairAnnot")
        @test description(pair_new.first) == "PairAnnot"
        @test description(pair_new.second) == "PairAnnot"
    end

    @testset "miniblast" begin
        # MSA with forward and reverse (reverse-complement) matches
        msa = MSA([
            "ATTGGTTCCCCCCCCCCCCCCAACCAAT",
            "ATTGGTACCCCCCCCCCCCCCTACCAAT"
        ])

        # 1. Exact match (should find 2 hits: forward and reverse-complement)
        hits = miniblast(msa, "ATTGGTT")
        @test length(hits) == 2
        @test :forward in [h.strand for h in hits]
        @test :reverse in [h.strand for h in hits]

        # 2. High threshold (no partial matches)
        hits_none = miniblast(msa, "ATTGGTT", 0.99)
        @test isempty(hits_none)

        # 3. Low threshold
        hits_partial = miniblast(msa, "ATTGGTT", 0.8)
        @test !isempty(hits_partial)
        @test all(h.identity ≈ 0.9285714285714286 for h in hits_partial)

        # 4. Query as a Primer object
        p = Primer(msa, 1:8; descr="query")
        hits_p = miniblast(msa, p)
        @test length(hits_p) == 2

        # 5. Edge cases
        @test isempty(miniblast(msa, "")) # Empty query
        @test isempty(miniblast(msa, "ACGTACGTACGTACGTACGTACGT")) # Too long query
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

    @testset "construct_primers Advanced" begin
        # Use MSA from previous tests
        seq1 = "GATCTGTAATGAGCGGCAGACCGACCGCGAATTAGACCTCGCCGAAGCCCTGGCCGCCAAGCTCAATTCGAAGCTCATTCAC----"
        seq2 = "CATTTGCAACGAGCGTCAGACCGACCGCGAACTCGACCTGGCCGAAGCGCTGGCTGCCAAACTCAATTCTAAGCTCATC-------"
        seq3 = "CATTTGTAACGAGCGTCAGACCGACCGTGAACTCGACCTCGCCGA-------GCTGCCAAACTCAATTCCAAGCTCATCCACTT--"
        seq4 = "---CTGTAACGAGCGGCAGACTGACCGAGAATTAGACCTCGCTGAAGCGCTGGCCGCCAAGCTCAATTCGAAGCTCATTCACTTTG"
        seq5 = "----TGTAACGAGCGGCAGACTGACCGAGAATTAGACCTCGCTGAAGCGCTGGCCGCCAAGCTCAATTCGAAGCTCATTCACTTAG"
        msa_1 = MSA([seq1, seq2, seq3, seq4, seq5])

        @testset "Negative MSA & Adapters" begin
            # 1. negative_msa validation branch
            @test_throws ArgumentError construct_primers(msa_1; negative_msa=[(msa_1, 1.5)])
            
            # 2. _evaluate filtering branch: if !isempty(negative_msa)
            # Use PERFECT MSA (no gaps or degeneracies).
            # Any primer generated from it will 100% match it.
            msa_clean = MSA(["ACGTACGTACGTACGTACGT", "ACGTACGTACGTACGTACGT"])
            
            # Pass it as negative with threshold 1.0 (strict 100% match)
            primers_neg = construct_primers(msa_clean; negative_msa=[(msa_clean, 1.0)])
            @test isempty(primers_neg) # All primers should be rejected

            # 3. if !isnothing(adapter_pair) branch - rejection by max_dg_drop
            adp = Oligo("GCGCGC", "FwdAdp") => Oligo("GCGCGC", "RevAdp")
            
            # max_dg_drop = -100.0 means that the condition (final_dg - dg_val) < 100.0
            # will trigger for any primer (since the ΔG difference is actually units).
            # Thus we test the rejection mechanics itself, independent of thermodynamics.
            primers_adp_strict = construct_primers(msa_1; adapter_pair=adp, offtarget_reject_threshold=0.9, max_dg_drop=-100.0)
            @test isempty(primers_adp_strict)
            
            # max_dg_drop = 100.0 means that the condition (final_dg - dg_val) < -100.0
            # will never trigger, so all primers with adapter will pass the filter.
            primers_adp = construct_primers(msa_1; adapter_pair=adp, offtarget_reject_threshold=0.9, max_dg_drop=100.0)
            @test !isempty(primers_adp)
            @test all(!isnothing(p.adapter) for p in primers_adp)
            @test all(p.is_forward ? p.adapter == adp.first : p.adapter == adp.second for p in primers_adp)
        end
    end

    @testset "best_pairs Advanced" begin
        msa_sort = MSA(["ACGTACGTACGTACGTACGT", "ACGTACGTACGTACGTACGT"])
        tm1 = (mean = 55.0, conf = (53.0, 57.0), min = 53.0, max = 57.0)
        tm2 = (mean = 56.0, conf = (54.0, 58.0), min = 54.0, max = 58.0)
        
        # 4 primers: 2 forward, 2 reverse
        f1 = Primer{DegenOligo}(msa_sort, 1:4, true, Oligo("ACGT"), 3, tm1, -5.0, 0.5, 0.0)
        f2 = Primer{DegenOligo}(msa_sort, 5:8, true, Oligo("ACGT"), 3, tm2, -5.0, 0.5, 0.0)
        r1 = Primer{DegenOligo}(msa_sort, 13:16, false, Oligo("ACGT"), 3, tm1, -5.0, 0.5, 0.0)
        r2 = Primer{DegenOligo}(msa_sort, 17:20, false, Oligo("ACGT"), 3, tm2, -5.0, 0.5, 0.0)
        primers_all = [f1, f2, r1, r2]

        @testset "Sorting" begin
            # Sort by :startpos (f1=1, f2=5)
            pairs_startpos = best_pairs(primers_all; sortby=:startpos)
            @test first(pairs_startpos[1].first.pos) < first(pairs_startpos[end].first.pos)
            
            # Sort by :length (f2->r1 = 12bp, f1->r2 = 20bp)
            pairs_length = best_pairs(primers_all; sortby=:length)
            len1 = last(pairs_length[1].second.pos) - first(pairs_length[1].first.pos) + 1
            len_end = last(pairs_length[end].second.pos) - first(pairs_length[end].first.pos) + 1
            @test len1 < len_end
        end

        @testset "Nested Pairs (best_pairs)" begin
            # Flanking pair: 3:6 and 15:18
            f_flank = Primer{DegenOligo}(msa_sort, 3:6, true, Oligo("ACGT"), 3, tm1, -5.0, 0.5, 0.0)
            r_flank = Primer{DegenOligo}(msa_sort, 15:18, false, Oligo("ACGT"), 3, tm1, -5.0, 0.5, 0.0)
            flank_pair = f_flank => r_flank

            # 1. MSA validation for flanking pair
            msa_other = MSA(["TTTTTTTT", "TTTTTTTT"])
            f_other = Primer{DegenOligo}(msa_other, 1:4, true, Oligo("TTTT"), 3, tm1, -5.0, 0.5, 0.0)
            r_other = Primer{DegenOligo}(msa_other, 5:8, false, Oligo("AAAA"), 3, tm1, -5.0, 0.5, 0.0)
            @test_throws ArgumentError best_pairs(primers_all; nested_pair=(f_other => r_other, 0))
            
            # 2. Validation: amp_start > amp_stop
            f_oob = Primer{DegenOligo}(msa_sort, 15:18, true, Oligo("ACGT"), 3, tm1, -5.0, 0.5, 0.0)
            r_oob = Primer{DegenOligo}(msa_sort, 1:4, false, Oligo("ACGT"), 3, tm1, -5.0, 0.5, 0.0)
            @test_throws ArgumentError best_pairs(primers_all; nested_pair=(f_oob => r_oob, 0))

            # 3. Validation: offset < 0 out of bounds (nesting region goes negative)
            @test_throws ArgumentError best_pairs(primers_all; nested_pair=(flank_pair, -10))
            # 4. Validation: offset < 0 amplicon_len too large
            @test_throws ArgumentError best_pairs(primers_all; nested_pair=(flank_pair, -1), amplicon_len=20:30)

            # 5. Validation: offset > 0 amplicon_len too small (min amplicon larger than theoretical)
            @test_throws ArgumentError best_pairs(primers_all; nested_pair=(flank_pair, 5), amplicon_len=1:10)

            # 6. Filtering: offset < 0 (Inner pairs)
            # bounds: 3-(-1)=4, 18+(-1)=17.
            # f1(1:4) excluded (first=1 < 4). f2(5:8) included.
            # r1(13:16) included. r2(17:20) excluded (last=20 > 17).
            inner_pairs = best_pairs(primers_all; nested_pair=(flank_pair, -1))
            @test length(inner_pairs) == 1
            @test inner_pairs[1].first === f2
            @test inner_pairs[1].second === r1

            # 7. Filtering: offset > 0 (Outer pairs)
            msa_outer = MSA([repeat("ACGT", 10), repeat("ACGT", 10)]) # length 40
            f1_out = Primer{DegenOligo}(msa_outer, 1:4, true, Oligo("ACGT"), 3, tm1, -5.0, 0.5, 0.0)
            f2_out = Primer{DegenOligo}(msa_outer, 16:19, true, Oligo("ACGT"), 3, tm2, -5.0, 0.5, 0.0)
            r1_out = Primer{DegenOligo}(msa_outer, 22:25, false, Oligo("ACGT"), 3, tm1, -5.0, 0.5, 0.0)
            r2_out = Primer{DegenOligo}(msa_outer, 37:40, false, Oligo("ACGT"), 3, tm2, -5.0, 0.5, 0.0)
            primers_outer = [f1_out, f2_out, r1_out, r2_out]
            
            f_flank_out = Primer{DegenOligo}(msa_outer, 10:13, true, Oligo("ACGT"), 3, tm1, -5.0, 0.5, 0.0)
            r_flank_out = Primer{DegenOligo}(msa_outer, 28:31, false, Oligo("ACGT"), 3, tm1, -5.0, 0.5, 0.0)
            flank_pair_out = f_flank_out => r_flank_out
            
            # bounds: 10-5=5, 28+5=33.
            # f1(1:4) included (last=4 < 5). f2(16:19) excluded (last=19 >= 5).
            # r1(22:25) excluded (first=22 <= 33). r2(37:40) included (first=37 > 33).
            outer_pairs = best_pairs(primers_outer; nested_pair=(flank_pair_out, 5))
            @test length(outer_pairs) == 1
            @test outer_pairs[1].first === f1_out
            @test outer_pairs[1].second === r2_out
        end
    end
end

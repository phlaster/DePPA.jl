using Test
using DEPPA.Alignments
using DEPPA.Oligos
using Random

@testset "Alignments Tests" begin
    
    @testset "MSA Constructors" begin
        # Basic construction from strings
        seqs = ["ACGT", "ACGT"]
        msa = MSA(seqs)
        @test msa isa MSA
        
        # Empty MSA
        msa_empty = MSA(String[])
        @test msa_empty isa MSA
        @test isempty(msa_empty.seqs)
        
        # Lowercase handling
        msa_lower = MSA(["acgt", "tgca"])
        @test String(getsequence(msa_lower, 1)) == "ACGT"
        
        # Length mismatch
        @test_throws ArgumentError MSA(["ACGT", "ACG"])
        
        # Negative bootstrap
        @test_throws ArgumentError MSA(["ACGT"]; bootstrap=-1)
        
        # Bootstrap and seed
        msa_boot = MSA(["ACGT", "TGCA", "GGGG"]; bootstrap=5, seed=123)
        @test bval(msa_boot) == 5
        @test size(get_base_count(msa_boot)) == (4, 4)
    end

    @testset "Fasta I/O" begin
        tmpfile = tempname() * ".fasta"
        try
            # Write test fasta
            open(tmpfile, "w") do f
                println(f, ">seq1")
                println(f, "ACGT")
                println(f, ">seq2")
                println(f, "TGCA")
                println(f, ">seq3")
                println(f, "GGGG")
            end
            
            msa = MSA(tmpfile)
            @test nseqs(msa) == 3
            @test width(msa) == 4
            
            # Predicate filtering
            msa_f = MSA(desc -> desc == "seq1", tmpfile)
            @test nseqs(msa_f) == 1
            @test String(getsequence(msa_f, 1)) == "ACGT"
            
            # Invalid characters
            open(tmpfile, "w") do f
                println(f, ">bad")
                println(f, "ACGT*")
            end
            @test_throws ArgumentError MSA(tmpfile)
            
            # mafft=true without MAFFT loaded
            open(tmpfile, "w") do f
                println(f, ">seq1")
                println(f, "ACGT")
            end
            @test_throws ErrorException MSA(tmpfile; mafft=true)
            
        finally
            rm(tmpfile, force=true)
        end
    end

    @testset "Basic Properties" begin
        msa = MSA(["ACGT", "TGCA"])
        @test nseqs(msa) == 2
        @test width(msa) == 4
        @test height(msa) == 2
        @test length(msa) == 4
        @test size(msa) == (2, 4)
        @test size(msa, 1) == 2
        @test size(msa, 2) == 4
        @test axes(msa, 1) == Base.OneTo(2)
        @test axes(msa, 2) == Base.OneTo(4)
        @test lastindex(msa) == 4
        @test lastindex(msa, 1) == 2
        @test ndims(msa) == 2
        
        # Empty MSA properties
        msa_empty = MSA(String[])
        @test nseqs(msa_empty) == 0
        @test width(msa_empty) == 0
        @test height(msa_empty) == 0
        @test length(msa_empty) == 0
        @test size(msa_empty) == (0, 0)
    end

    @testset "Sequences & Base Counts" begin
        msa = MSA(["ACGT", "TGCA"])
        
        # getsequence
        @test getsequence(msa, 1) isa GappedOligo
        @test String(getsequence(msa, 1)) == "ACGT"
        @test getsequence(msa, 1, 1) == 'A'
        @test getsequence(msa, 2, 4) == 'A'
        
        # get_base_count
        # A=[1,0,0,0], T=[0,0,0,1] -> avg=[0.5, 0, 0, 0.5]
        @test get_base_count(msa, 1) ≈ [0.5, 0.0, 0.0, 0.5]
        @test get_base_count(msa, 1:2) isa AbstractMatrix{Float64}
        @test get_base_count(msa) isa AbstractMatrix{Float64}
        
        # Gaps
        msa_gap = MSA(["A---", "ACGT"])
        # col 2: -, C -> [0, 0.5, 0, 0]
        @test get_base_count(msa_gap, 2) ≈ [0.0, 0.5, 0.0, 0.0]
        
        # Degenerate bases
        msa_deg = MSA(["ACGT", "RNWS"])
        # col 1: A vs R(0.5, 0, 0.5, 0) -> [0.75, 0, 0.25, 0]
        @test get_base_count(msa_deg, 1) ≈ [0.75, 0.0, 0.25, 0.0]
    end

    @testset "Depth & Determinacy" begin
        msa = MSA(["ACGT", "ACGT"])
        @test msadepth(msa, 1) ≈ 1.0
        @test msadepth(msa, 1:4) ≈ [1.0, 1.0, 1.0, 1.0]
        @test msadepth(msa) ≈ [1.0, 1.0, 1.0, 1.0]
        
        @test msadet(msa, 1) ≈ 1.0
        @test msadet(msa, 1:4) ≈ [1.0, 1.0, 1.0, 1.0]
        @test msadet(msa) ≈ [1.0, 1.0, 1.0, 1.0]

        msa_gap = MSA(["A---", "ACGT"])
        @test msadepth(msa_gap, 1) ≈ 1.0
        @test msadepth(msa_gap, 2) ≈ 0.5
        @test msadet(msa_gap, 1) ≈ 1.0
        @test msadet(msa_gap, 2) ≈ 1.0
        @test msadepth(msa_gap) ≈ [1.0, 0.5, 0.5, 0.5]
        @test msadet(msa_gap) ≈ [1.0, 1.0, 1.0, 1.0]

        msa_mixed = MSA(["AC", "AC", "GT"])
        @test msadepth(msa_mixed, 1) ≈ 1.0
        @test msadet(msa_mixed, 1) ≈ 2/3
        
        # Gap-only column
        msa_gap_only = MSA(["---", "---"])
        @test msadepth(msa_gap_only, 1) ≈ 0.0
        @test msadet(msa_gap_only, 1) ≈ 0.0
        
        # Empty
        msa_empty = MSA(String[])
        @test msadepth(msa_empty) == Float64[]
        @test msadet(msa_empty) == Float64[]
    end

    @testset "Consensus" begin
        msa = MSA(["ACGT", "ACGT"])
        @test consensus_major(msa, 1) == 'A'
        @test consensus_major(msa, 1:4) isa GappedOligo
        @test String(consensus_major(msa, 1:4)) == "ACGT"
        @test String(consensus_major(msa)) == "ACGT"
        
        @test consensus_degen(msa, 1) == 'A'
        @test String(consensus_degen(msa)) == "ACGT"
        
        msa_mixed = MSA(["AC", "AG"])
        @test consensus_major(msa_mixed, 2) == 'C' # C and G both 0.5, argmax gives first -> C
        @test consensus_degen(msa_mixed, 2) == 'S' # S = C/G
        @test String(consensus_degen(msa_mixed)) == "AS"
        
        # slack parameter
        @test consensus_degen(msa_mixed, 2; slack=0.0) == 'S'
        @test consensus_degen(msa_mixed, 2; slack=0.6) == '-'
        
        # Gap-only column
        msa_gap = MSA(["---", "A--"])
        @test consensus_major(msa_gap, 1) == 'A'
        @test consensus_major(msa_gap, 2) == '-'
        @test consensus_degen(msa_gap, 2) == '-'
    end

    @testset "Dry MSA" begin
        msa = MSA(["ACGT", "ACGT"])
        dry = dry_msa(msa)
        @test nseqs(dry) == 2
        @test String(getsequence(dry, 1)) == "ACGT"
        
        msa_gap = MSA(["A---", "ACGT", "----"])
        dry = dry_msa(msa_gap; gap_content=0.8)
        # Row 1: prop=0.75 < 0.8 -> kept
        # Row 2: prop=0.0 < 0.8 -> kept
        # Row 3: prop=1.0 < 0.8 -> dropped
        @test nseqs(dry) == 2
        
        # Drop all-gap columns
        msa_gap_cols = MSA(["-A", "-C"])
        dry2 = dry_msa(msa_gap_cols)
        @test nseqs(dry2) == 2
        @test String(getsequence(dry2, 1)) == "A"
        @test width(dry2) == 1
        
        # Completely empty gap MSA
        msa_all_gaps = MSA(["----", "----"])
        dry3 = dry_msa(msa_all_gaps)
        @test dry3 isa MSAView
        @test nseqs(dry3) == 2
        @test width(dry3) == 0
    end

    @testset "Nucleotide Diversity" begin
        # Identical
        msa_id = MSA(["ACGT", "ACGT"])
        @test nucleotide_diversity(msa_id) ≈ 0.0

        # All different
        msa_diff = MSA(["ACGT", "TGCA"])
        @test nucleotide_diversity(msa_diff) ≈ 1.0

        # With gaps, ignore_gaps=true
        msa_gap = MSA(["ACGT", "A-GT"])
        @test DEPPA.Alignments._pairwise_distance(msa_gap, 1, 2; ignore_gaps=true) ≈ 0.0
        @test nucleotide_diversity(msa_gap; ignore_gaps=true) ≈ 0.0

        # With gaps, ignore_gaps=false
        @test DEPPA.Alignments._pairwise_distance(msa_gap, 1, 2; ignore_gaps=false) ≈ 0.25
        @test nucleotide_diversity(msa_gap; ignore_gaps=false) ≈ 0.25

        # Degenerate bases
        msa_deg = MSA(["ACGT", "RNWS"])
        # A vs R -> diff 0.5, C vs N -> diff 0.75, G vs W -> diff 1.0, T vs S -> diff 1.0
        # Total diff = 3.25. Sites = 4. dist = 3.25/4 = 0.8125
        @test nucleotide_diversity(msa_deg) ≈ 0.8125
        
        # Identical degenerate bases short-circuit to exact match in the code
        msa_deg_same = MSA(["R", "R"])
        @test nucleotide_diversity(msa_deg_same) ≈ 0.0

        # Different degenerate bases with overlapping probabilities
        msa_deg_overlap = MSA(["R", "S"]) # R=(A/G), S=(C/G). Match prob = 0.25, diff = 0.75
        @test nucleotide_diversity(msa_deg_overlap) ≈ 0.75

        # Empty / Single
        @test nucleotide_diversity(MSA(String[])) == 0.0
        @test nucleotide_diversity(MSA(["ACGT"])) == 0.0

        # Large MSA sampling
        seqs_large = [Random.randstring("ACGT", 10) for _ in 1:250]
        msa_large = MSA(seqs_large)
        @test nucleotide_diversity(msa_large; max_pairs=100) isa Float64
    end

    @testset "MSAView & Indexing" begin
        msa = MSA(["ACGT", "TGCA"])
        
        # Basic slicing
        view = msa[1:2, 2:3]
        @test view isa MSAView
        @test root(view) === msa
        @test bval(view) == bval(msa)
        @test nseqs(view) == 2
        @test width(view) == 2
        @test height(view) == 2
        @test size(view) == (2, 2)
        
        # Sequences from view
        @test getsequence(view, 1) isa OligoView
        @test String(getsequence(view, 1)) == "CG"
        @test String(getsequence(view, 2)) == "GC"
        @test getsequence(view, 1, 1) == 'C'
        @test getsequence(view, 2, 2) == 'C'
        
        # Base counts from full height view
        @test get_base_count(view, 1) == get_base_count(msa, 2)
        @test get_base_count(view, 1:2) == get_base_count(msa, 2:3)
        @test get_base_count(view) == get_base_count(msa, 2:3)
        
        # Base counts from partial height view (should throw)
        view_row = msa[1:1, 1:4]
        @test_throws ErrorException get_base_count(view_row, 1)
        
        # Depth and determinacy on view
        @test msadepth(view, 1) == msadepth(msa, 2)
        @test msadet(view, 1) == msadet(msa, 2)
        
        # Consensus on view
        cons = consensus_major(view)
        @test length(cons) == 2
        
        # Colon indexing
        @test msa[1, :] isa GappedOligo
        @test String(msa[1, :]) == "ACGT"
        @test msa[:, 1] isa MSAView
        @test msa[:, :] isa MSAView
        
        # Single element
        @test msa[1, 1] == 'A'
        @test msa[2, 4] == 'A'
        
        # View of a view
        view2 = view[1:2, 1:1]
        @test root(view2) === msa
        @test String(getsequence(view2, 1)) == "C"
    end
end
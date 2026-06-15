using DEPPA
using Test
using Aqua
using JET
using Random
using MAFFT_jll

@testset "DEPPA.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(DEPPA)
    end
    @testset "Code linting (JET.jl)" begin
        JET.test_package(DEPPA; target_defined_modules = true)
    end
    
    @testset "Oligos" begin
        @testset "Oligo" begin
            # Construction tests
            @test Oligo("ACGT") isa Oligo
            @test Oligo("ACGT", "test description") isa Oligo
            @test Oligo("ACGT", "") isa Oligo
            @test Oligo("ACGT", nothing) isa Oligo
            @test Oligo(['A', 'C', 'G', 'T']) isa Oligo
            @test Oligo(['A', 'C', 'G', 'T'], "test description") isa Oligo
            @test Oligo("ACGT", 123) isa Oligo
            @test_throws ErrorException Oligo("ACGTX")
            @test Oligo() == Oligo("")
            @test Oligo("", "empty") == Oligo("", "empty")
            
            # String interface tests
            oligo = Oligo("ACGT", "test")
            @test String(oligo) == "ACGT"
            @test length(oligo) == 4
            @test isempty(Oligo("")) == true
            @test isempty(oligo) == false
            @test collect(oligo) == ['A', 'C', 'G', 'T']
            @test oligo[1] == 'A'
            @test oligo[2:3] == "CG"
            @test occursin("CG", oligo)
            @test first(oligo) == 'A'
            @test last(oligo) == 'T'
            @test oligo[4:-1:1] == "TGCA"
            
            # Description tests
            @test description(oligo) == "test"
            @test description(Oligo("ACGT", "")) == ""
            
            # Conversion tests
            deg_oligo = DegenerateOligo("ACGT")
            @test Oligo(deg_oligo) isa Oligo
            @test_throws InexactError Oligo(DegenerateOligo("ACGN"))
            @test convert(Oligo, deg_oligo) isa Oligo
            
            # Concatenation tests
            oligo2 = Oligo("TGCA")
            @test (oligo * oligo2) == Oligo("ACGTTGCA", "concat")
            @test description((oligo * oligo2)) == "concat"
            @test (Oligo("A") * Oligo("C") * Oligo("G") * Oligo("T")) == oligo
            
            # Equality tests
            @test Oligo("ACGT") == Oligo("ACGT")
            @test Oligo("ACGT") != Oligo("TGCA")
            @test Oligo("ACGT") == "ACGT"
            @test "ACGT" == Oligo("ACGT")
            @test Oligo("ACGT") != "TGCA"
            
            # Case handling
            @test Oligo("acgt") == Oligo("ACGT")
            
            # Empty sequence tests
            empty_oligo = Oligo("", "empty")
            @test isempty(empty_oligo) == true
            @test length(empty_oligo) == 0
            @test String(empty_oligo) == ""
            @test_throws BoundsError empty_oligo[1]

            @testset "Display" begin
                oligo = Oligo("AGTC", "descr")
                oligo_nodesc = Oligo("AGTC", "")
                oligo_long = Oligo("A"^25, "")
                # Short display
                @test sprint(show, oligo) == "Oligo(\"AGTC\", len=4, desc=\"descr\")"
                @test sprint(show, oligo_nodesc) == "Oligo(\"AGTC\", len=4)"
                @test sprint(show, oligo_long) == "Oligo(\"AAAAAAAAAAAAAAAAA...\", len=25)"
                # Full display
                @test sprint(show, MIME"text/plain"(), oligo) == "Oligo\n  Sequence: AGTC\n  Length: 4\n  Description: \"descr\"\n"
                @test sprint(show, MIME"text/plain"(), oligo_nodesc) == "Oligo\n  Sequence: AGTC\n  Length: 4\n  Description: (none)\n"
                @test sprint(show, MIME"text/plain"(), oligo_long) == "Oligo\n  Sequence: AAAAAAAAAAAAAAAAAAAAAAAAA\n  Length: 25\n  Description: (none)\n"
            end
        end
        
        @testset "DegenerateOligo" begin
            # Construction tests
            @test DegenerateOligo("ACGT") isa DegenerateOligo
            @test DegenerateOligo("ACGN") isa DegenerateOligo
            @test DegenerateOligo("ACGN", "test description") isa DegenerateOligo
            @test DegenerateOligo("ACGN", "") isa DegenerateOligo
            @test DegenerateOligo("ACGN", nothing) isa DegenerateOligo
            @test DegenerateOligo(['A', 'C', 'G', 'N']) isa DegenerateOligo
            @test DegenerateOligo(['A', 'C', 'G', 'N'], "test description") isa DegenerateOligo
            @test DegenerateOligo("ACGT", 123) isa DegenerateOligo
            @test_throws ErrorException DegenerateOligo("ACGTX")
            @test DegenerateOligo() == DegenerateOligo("")
            @test DegenerateOligo("", "empty") == DegenerateOligo("", 0, 1, "empty")
            
            # String interface tests
            deg_oligo = DegenerateOligo("ACGN", "test")
            @test String(deg_oligo) == "ACGN"
            @test length(deg_oligo) == 4
            @test isempty(DegenerateOligo(""))
            @test !isempty(deg_oligo)
            @test collect(deg_oligo) == ['A', 'C', 'G', 'N']
            @test deg_oligo[1] == 'A'
            @test deg_oligo[2:3] == "CG"
            @test occursin("CG", deg_oligo)
            @test first(deg_oligo) == 'A'
            @test last(deg_oligo) == 'N'
            @test deg_oligo[4:-1:1] == "NGCA"
            
            # Properties tests
            @test n_deg_pos(deg_oligo) == 1
            @test n_unique_oligos(deg_oligo) == 4
            @test n_deg_pos(DegenerateOligo("ACGT")) == 0
            @test n_unique_oligos(DegenerateOligo("ACGT")) == 1
            @test n_deg_pos(DegenerateOligo("NNNN")) == 4
            @test n_unique_oligos(DegenerateOligo("NNNN")) == 256
            @test n_deg_pos(DegenerateOligo("RY")) == 2
            @test n_unique_oligos(DegenerateOligo("RY")) == 4
            @test n_deg_pos(DegenerateOligo("BVDH")) == 4
            @test n_unique_oligos(DegenerateOligo("BVDH")) == 81
            
            # Description tests
            @test description(deg_oligo) == "test"
            @test description(DegenerateOligo("ACGN", "")) == ""
            
            # Conversion tests
            oligo = Oligo("ACGT")
            @test DegenerateOligo(oligo) isa DegenerateOligo
            @test convert(DegenerateOligo, oligo) isa DegenerateOligo
            @test String(DegenerateOligo(oligo)) == "ACGT"
            @test n_deg_pos(DegenerateOligo(oligo)) == 0
            @test n_unique_oligos(DegenerateOligo(oligo)) == 1
            
            # Concatenation tests
            deg_oligo2 = DegenerateOligo("NNGT")
            @test (deg_oligo * deg_oligo2) isa DegenerateOligo
            @test String((deg_oligo * deg_oligo2)) == "ACGNNNGT"
            @test n_deg_pos((deg_oligo * deg_oligo2)) == 3
            @test n_unique_oligos((deg_oligo * deg_oligo2)) == 4^3
            
            @test (oligo * deg_oligo) isa DegenerateOligo
            @test String((oligo * deg_oligo)) == "ACGTACGN"
            @test n_deg_pos((oligo * deg_oligo)) == 1
            @test n_unique_oligos((oligo * deg_oligo)) == 4
            
            @test (deg_oligo * oligo) isa DegenerateOligo
            @test String((deg_oligo * oligo)) == "ACGNACGT"
            @test n_deg_pos((deg_oligo * oligo)) == 1
            @test n_unique_oligos((deg_oligo * oligo)) == 4
            
            # NonDegenIterator tests
            @test length(nondegens(deg_oligo)) == 4
            variants = collect(nondegens(deg_oligo))
            @test length(variants) == 4
            @test Set(String.(variants)) == Set(["ACGA", "ACGC", "ACGG", "ACGT"])
            
            # Complex degenerate sequence
            complex_deg = DegenerateOligo("RYSWKMBDHVN")
            @test n_deg_pos(complex_deg) == 11
            @test n_unique_oligos(complex_deg) == 20736
            
            # NonDegenIterator for complex sequence
            @test length(nondegens(complex_deg)) == 20736
            complex_variants = collect(Iterators.take(nondegens(complex_deg), 5))
            @test all(x -> length(x) == 11, complex_variants)
            
            # Empty sequence
            empty_deg = DegenerateOligo("", "empty")
            @test isempty(empty_deg) == true
            @test length(empty_deg) == 0
            @test String(empty_deg) == ""
            @test n_deg_pos(empty_deg) == 0
            @test n_unique_oligos(empty_deg) == 1
            @test length(nondegens(empty_deg)) == 1

            @testset "Display" begin
                doligo = DegenerateOligo("AGNTC", "descr")
                doligo_nodesc = DegenerateOligo("AGNTC", "")
                doligo_long = DegenerateOligo("N"^25, "")
                # Short display
                @test sprint(show, doligo) == "DegenerateOligo(\"AGNTC\", len=5, n_deg=1, vars=4, desc=\"descr\")"
                @test sprint(show, doligo_nodesc) == "DegenerateOligo(\"AGNTC\", len=5, n_deg=1, vars=4)"
                @test sprint(show, doligo_long) == "DegenerateOligo(\"NNNNNNNNNNNNNNNNN...\", len=25, n_deg=25, vars=>10k)"
                # Full display
                @test sprint(show, MIME"text/plain"(), doligo) == "DegenerateOligo\n  Sequence: AGNTC\n  Length: 5\n  Degenerate positions: 1\n  Unique variants: 4\n  Description: \"descr\"\n"
                @test sprint(show, MIME"text/plain"(), doligo_nodesc) == "DegenerateOligo\n  Sequence: AGNTC\n  Length: 5\n  Degenerate positions: 1\n  Unique variants: 4\n  Description: (none)\n"
                @test sprint(show, MIME"text/plain"(), doligo_long) == "DegenerateOligo\n  Sequence: NNNNNNNNNNNNNNNNNNNNNNNNN\n  Length: 25\n  Degenerate positions: 25\n  Unique variants: 1125899906842624\n  Description: (none)\n"
            end
        end
        
        @testset "OligoView" begin
            # Basic view creation
            oligo = Oligo("ACGTACGT", "test")
            deg_oligo = DegenerateOligo("ACGNACGN", "deg test")
            
            view1 = oligo[2:5]
            @test view1 isa OligoView{Oligo}
            @test String(view1) == "CGTA"
            @test length(view1) == 4
            @test collect(view1) == ['C', 'G', 'T', 'A']
            @test description(view1) == "test"
            
            view2 = deg_oligo[2:5]
            @test view2 isa OligoView{DegenerateOligo}
            @test String(view2) == "CGNA"
            @test length(view2) == 4
            @test collect(view2) == ['C', 'G', 'N', 'A']
            @test description(view2) == "deg test"
            
            # Edge cases
            @test isempty(oligo[1:0])
            @test_throws BoundsError oligo[10:11]
            
            # String interface
            @test view1[1] == 'C'
            @test view1[2:3] == "GT"
            @test view1[end] == 'A'
            @test view1[3:-1:1] == "TGC"
            @test occursin("GT", view1)
            
            # Properties
            @test n_deg_pos(view1) == 0
            @test n_unique_oligos(view1) == 1
            @test n_deg_pos(view2) == 1
            @test n_unique_oligos(view2) == 4
            
            # Concatenation
            @test (view1 * Oligo("TG")) isa Oligo
            @test String(view1 * Oligo("TG")) == "CGTATG"
            
            @test (view2 * Oligo("TG")) isa DegenerateOligo
            @test String(view2 * Oligo("TG")) == "CGNATG"
            @test n_deg_pos(view2 * Oligo("TG")) == 1
            @test n_unique_oligos(view2 * Oligo("TG")) == 4
            
            @test (Oligo("AT") * view1) isa Oligo
            @test String(Oligo("AT") * view1) == "ATCGTA"
            
            @test (Oligo("AT") * view2) isa DegenerateOligo
            @test String(Oligo("AT") * view2) == "ATCGNA"
            @test n_deg_pos(Oligo("AT") * view2) == 1
            @test n_unique_oligos(Oligo("AT") * view2) == 4
            
            # Conversion
            @test convert(Oligo, view1) isa Oligo
            @test String(convert(Oligo, view1)) == "CGTA"
            @test_throws InexactError convert(Oligo, view2)
            
            # NonDegenIterator for views
            @test length(nondegens(view2)) == 4
            view_variants = collect(nondegens(view2))
            @test length(view_variants) == 4
            @test Set(String.(view_variants)) == Set(["CGAA", "CGCA", "CGGA", "CGTA"])
            
            # Empty view
            empty_view = oligo[1:0]
            @test isempty(empty_view)
            @test length(empty_view) == 0
            @test collect(empty_view) == Char[]
            @test n_deg_pos(empty_view) == 0
            @test n_unique_oligos(empty_view) == 1
            
            # View of a view
            subview = view1[2:3]
            @test String(subview) == "GT"
            @test description(subview) == "test"
            @test subview isa OligoView{Oligo}
            @test n_deg_pos(subview) == 0
            @test n_unique_oligos(subview) == 1
            
            # View of degenerate view
            deg_subview = view2[2:3]
            @test String(deg_subview) == "GN"
            @test description(deg_subview) == "deg test"
            @test deg_subview isa OligoView{DegenerateOligo}
            @test n_deg_pos(deg_subview) == 1
            @test n_unique_oligos(deg_subview) == 4

            @testset "Display" begin
                doligo_view = DegenerateOligo("AGNTC", "descr")[2:3]
                doligo_view_nodesc = DegenerateOligo("AGNTC", "")[2:4]
                doligo_view_long = DegenerateOligo("N"^25, "")[10:22]
                # Short display
                @test sprint(show, doligo_view) == "OligoView(\"GN\", len=2, range=2:3, desc=\"descr\")"
                @test sprint(show, doligo_view_nodesc) == "OligoView(\"GNT\", len=3, range=2:4)"
                @test sprint(show, doligo_view_long) == "OligoView(\"NNNNNNNNNNNNN\", len=13, range=10:22)"
                # Full display
                @test sprint(show, MIME"text/plain"(), doligo_view) == "OligoView{DegenerateOligo}\n  Viewed sequence: GN\n  Length: 2\n  Range: 2:3\n  Parent description: descr\n"
                @test sprint(show, MIME"text/plain"(), doligo_view_nodesc) == "OligoView{DegenerateOligo}\n  Viewed sequence: GNT\n  Length: 3\n  Range: 2:4\n  Parent description: \n"
                @test sprint(show, MIME"text/plain"(), doligo_view_long) == "OligoView{DegenerateOligo}\n  Viewed sequence: NNNNNNNNNNNNN\n  Length: 13\n  Range: 10:22\n  Parent description: \n"
            end
        end
        
        @testset "GappedOligo" begin
            # Construction tests with Oligo parent
            seq = "ACGTACGT"
            gapped_seq = "AC--GT-ACGT"
            oligo = Oligo(seq, "test")
            go = GappedOligo(oligo, [3=>2, 5=>1])
            @test go isa GappedOligo{Oligo}
            @test String(go) == gapped_seq
            @test length(go) == 11
            @test parent(go) === oligo
            @test go.gaps == [3=>2, 5=>1]
            @test description(go) == description(oligo)
            @test hasgaps(go)
            @test n_deg_pos(go) == 0
            @test n_unique_oligos(go) == 1
            @test collect(go) == collect(gapped_seq)
            
            # Construction from string
            go_str = GappedOligo(gapped_seq, "test")
            @test go_str isa GappedOligo{Oligo}
            @test String(go_str) == gapped_seq
            @test parent(go_str) == Oligo(seq, "test")
            @test go_str.gaps == [3=>2, 5=>1]
            @test length(go_str) == 11
            
            # Construction with DegenerateOligo parent
            deg_seq = "ACGNRT"
            gapped_deg_seq = "AC--GN-RT"
            deg_oligo = DegenerateOligo(deg_seq, "deg test")
            go_deg = GappedOligo(deg_oligo, [3=>2, 5=>1])
            @test go_deg isa GappedOligo{DegenerateOligo}
            @test String(go_deg) == gapped_deg_seq
            @test length(go_deg) == 9
            @test parent(go_deg) === deg_oligo
            @test go_deg.gaps == [3=>2, 5=>1]
            @test n_deg_pos(go_deg) == 2
            @test n_unique_oligos(go_deg) == 8
            @test collect(go_deg) == collect(gapped_deg_seq)
            
            # Construction from degenerate string
            go_deg_str = GappedOligo(gapped_deg_seq, "deg test")
            @test go_deg_str isa GappedOligo{DegenerateOligo}
            @test String(go_deg_str) == gapped_deg_seq
            @test parent(go_deg_str) == DegenerateOligo(deg_seq, "deg test")
            @test go_deg_str.gaps == [3=>2, 5=>1]
            
            # Invalid gap positions
            @test_throws ArgumentError GappedOligo(oligo, [0=>1]) # Start < 1
            @test_throws ArgumentError GappedOligo(oligo, [10=>1]) # Start > parent len + 1
            @test_throws ArgumentError GappedOligo(oligo, [3=>0]) # Non-positive gap length
            @test_throws ArgumentError GappedOligo(oligo, [3=>2, 3=>1]) # Overlapping gaps
            
            # Edge cases: empty and no gaps
            empty_go = GappedOligo(Oligo(""), Pair{Int, Int}[])
            @test String(empty_go) == ""
            @test length(empty_go) == 0
            @test isempty(empty_go) == true
            @test hasgaps(empty_go) == false
            no_gap_go = GappedOligo(oligo, Pair{Int, Int}[])
            @test String(no_gap_go) == "ACGTACGT"
            @test length(no_gap_go) == 8
            @test hasgaps(no_gap_go) == false
            
            # Indexing tests
            @test go[1] == 'A'
            @test go[3] == '-'
            @test go[5] == 'G'
            @test go[11] == 'T'
            @test_throws BoundsError go[12]
            @test go[1:3] isa OligoView{GappedOligo{Oligo}}
            @test String(go[1:3]) == "AC-"
            @test String(go[3:4]) == "--"
            @test String(go[5:7]) == "GT-"
            @test String(go[8:11]) == "ACGT"
            @test String(go[11:11]) == "T"
            
            # Slicing with DegenerateOligo
            go_view = go_deg[2:7]
            @test go_view isa OligoView{GappedOligo{DegenerateOligo}}
            @test String(go_view) == "C--GN-"
            @test n_deg_pos(go_view) == 1
            @test n_unique_oligos(go_view) == 4
            @test collect(nondegens(go_view)) == GappedOligo.(["C--GA-", "C--GC-", "C--GG-", "C--GT-"])
            
            # Gaps at start/end
            go_start = GappedOligo(oligo, [1=>2])
            @test String(go_start) == "--ACGTACGT"
            @test length(go_start) == 10
            @test go_start[1] == '-' && go_start[3] == 'A'
            go_end = GappedOligo(oligo, [9=>2])
            @test String(go_end) == "ACGTACGT--"
            @test length(go_end) == 10
            @test go_end[8] == 'T' && go_end[9] == '-'
            
            # Complex slicing
            nested_view = go[2:8][2:5]
            @test nested_view isa OligoView{GappedOligo{Oligo}}
            @test nested_view == "--GT"
            @test go[4:5] == "-G"
            @test go[7:7] == "-"
            
            # Concatenation edge cases
            @test_throws ErrorException go * go
            go_view = go[1:3]
            @test_throws ErrorException go_view * Oligo("CG")
            
            # SeqFold methods
            @test SeqFold.revcomp(go) == "ACGT-AC--GT"
            @test SeqFold.revcomp(go) isa GappedOligo{Oligo}
            @test SeqFold.revcomp(go).gaps == [length(go)-a-1 => b for (a,b) in reverse(go.gaps)]
            @test SeqFold.complement(go) == "TG--CA-TGCA"
            @test SeqFold.complement(go).gaps == go.gaps
            @test SeqFold.gc_content(go) ≈ SeqFold.gc_content(oligo)
            @test_throws ErrorException SeqFold.fold(go)
            @test_throws ErrorException SeqFold.dg(go)
            @test_throws ErrorException SeqFold.tm(go)
            
            # Iteration
            @test collect(go) == collect(gapped_seq)
            @test collect(go_deg) == collect(gapped_deg_seq)
            empty_iter = GappedOligo(Oligo(""), Pair{Int, Int}[])
            @test iterate(empty_iter) === nothing
            
            # Nondegens for degenerate parent
            @test collect(nondegens(go)) == [go]
            @test collect(nondegens(go_deg)) == GappedOligo.(["AC--GA-AT", "AC--GA-GT", "AC--GC-AT", "AC--GC-GT", "AC--GG-AT", "AC--GG-GT", "AC--GT-AT", "AC--GT-GT"])
        end

        @testset "Extra" begin
            @testset "GappedOligo with DegenerateOligo" begin
                deg_oligo = DegenerateOligo("ACGNRT", "deg gap test")
                gaps = [3=>2, 5=>1]
                go = GappedOligo(deg_oligo, gaps)
                @test go isa GappedOligo{DegenerateOligo}
                @test String(go) == "AC--GN-RT"  # 6 bases + 3 gaps = 9 length
                @test parent(go) == deg_oligo
                @test go.gaps == gaps
                @test length(go) == 9
                @test n_deg_pos(go) == 2  # N and R from parent
                @test n_unique_oligos(go) == 8  # N=4, R=2
                @test collect(go) == ['A','C','-','-','G','N','-','R','T']
                @test String(go[2:7]) == "C--GN-"
                @test String(go[8:end]) == "RT"
                @test go[3] == '-'
                @test go[6] == 'N'
                go_view = go[3:6]
                @test go_view isa OligoView{GappedOligo{DegenerateOligo}}
                @test String(go_view) == "--GN"
                @test n_deg_pos(go_view) == 1
                @test n_unique_oligos(go_view) == 4
                @test collect(nondegens(go_view)) == GappedOligo.(["--GA", "--GC", "--GG", "--GT"])

                deg_oligo = DegenerateOligo("ACGNRT", "deg gap test")
                # Gaps at start
                go_start = GappedOligo(deg_oligo, [1=>2])
                @test String(go_start) == "--ACGNRT"
                @test length(go_start) == 8
                @test go_start[1] == '-' && go_start[2] == '-' && go_start[3] == 'A'
                # Gaps at end
                go_end = GappedOligo(deg_oligo, [7=>2])
                @test String(go_end) == "ACGNRT--"
                @test length(go_end) == 8
                @test go_end[6] == 'T' && go_end[7] == '-' && go_end[8] == '-'
                # Multiple gaps with degenerate bases
                go_multi = GappedOligo(deg_oligo, [2=>1, 3=>1, 4=>1])
                @test String(go_multi) == "A-C-G-NRT"
                @test length(go_multi) == 9
                @test collect(go_multi) == ['A', '-', 'C', '-', 'G', '-', 'N', 'R', 'T']
                # Slice including partial gap
                @test String(go_multi[2:5]) == "-C-G"
                # Gap-only slice
                go = GappedOligo(deg_oligo, [3=>2, 6=>1])
                @test String(go[3:4]) == "--"
            end

            @testset "Random Sampling" begin
                deg_oligo = DegenerateOligo("ACRN", "rand test")
                rng = Random.MersenneTwister(42)
                rand_oligo = rand(rng, deg_oligo)
                @test rand_oligo isa Oligo
                @test String(rand_oligo) in String.(nondegens(deg_oligo))
                @test description(rand_oligo) == "rand test"
                deg_view = deg_oligo[2:4]
                rand_view = rand(rng, deg_view)
                @test rand_view isa Oligo
                @test String(rand_view) in ["CAA", "CAG", "CGA", "CGG"]
                @test description(rand_view) == "rand test"
                # Test multiple samples to ensure coverage
                samples = Set(String(rand(rng, deg_oligo)) for _ in 1:100)
                @test length(samples) > 1  # Likely to hit multiple variants
            end

            @testset "Complex Slicing for GappedOligo" begin
                oligo = Oligo("ACGTACGT", "slice test")
                go = GappedOligo(oligo, [3=>2, 5=>1])
                @test String(go) == "AC--GT-ACGT"
                @test String(go[3:4]) == "--"  # Only gaps
                @test String(go[4:5]) == "-G"  # Partial gap
                @test String(go[6:8]) == "T-A"  # Gap in middle
                @test String(go[1:2]) == "AC"  # Before gaps
                @test String(go[end:end]) == "T"  # Single position
                nested_view = go[2:8][2:5]
                @test nested_view isa OligoView{GappedOligo{Oligo}}
                @test String(nested_view) == "--GT"  # Nested slice
                deg_oligo = DegenerateOligo("ACGNRT", "deg slice")
                go_deg = GappedOligo(deg_oligo, [3=>2])
                @test String(go_deg[2:5]) == "C--G"
                @test String(go_deg[3:4]) == "--"
            end

            @testset "Concatenation Edge Cases" begin
                go = GappedOligo(Oligo("ACGT", "test"), [3=>1])
                @test_throws ErrorException go * go  # Gapped concatenation not supported
                go_view = go[1:3]
                @test go_view isa OligoView{GappedOligo{Oligo}}
                @test_throws ErrorException go_view * Oligo("CG")
                deg_oligo = DegenerateOligo("ACN")
                go_deg = GappedOligo(deg_oligo, [2=>1])
                go_deg_view = go_deg[1:3]
                @test_throws ErrorException go_deg_view * Oligo("CG")
            end

            @testset "SeqFold Methods" begin
                deg_oligo = DegenerateOligo("ACGN", "tm test")
                tm_result = SeqFold.tm(deg_oligo, conditions=:pcr)
                @test tm_result.mean isa Float64
                @test tm_result.conf[1] <= tm_result.mean <= tm_result.conf[2]
                @test all(x in ["ACGA", "ACGC", "ACGG", "ACGT"] for x in nondegens(deg_oligo))
                deg_view = deg_oligo[2:4]
                tm_view = SeqFold.tm(deg_view, conditions=:pcr)
                @test tm_view.mean isa Float64
                @test tm_view.conf[1] <= tm_view.mean <= tm_view.conf[2]
                @test SeqFold.gc_content(deg_oligo) ≈ (0.0 + 1.0 + 1.0 + 0.5) / 4  # A,C,G,N weights
                @test SeqFold.gc_content(deg_view) ≈ (1.0 + 1.0 + 0.5) / 3  # C,G,N
                no_gap_go = GappedOligo(Oligo("ACGT"), Pair{Int, Int}[])
                @test SeqFold.tm(no_gap_go, conditions=:pcr).mean == SeqFold.tm("ACGT", conditions=:pcr)
            end

            @testset "Boundary Indexing and Slicing" begin
                # Single-position slices
                oligo = Oligo("ACGT", "single")
                @test oligo[1:1] == "A"
                @test oligo[end:end] == "T"
                deg_oligo = DegenerateOligo("ACGN", "deg single")
                @test deg_oligo[1:1] == "A"
                @test deg_oligo[end:end] == "N"
                go = GappedOligo(Oligo("ACGT"), [2=>1])
                @test go[2:2] == "-"  # Single gap position
                @test go[3:3] == "C"  # Single non-gap position

                # Negative step ranges
                @test oligo[4:-1:1] == "TGCA"
                @test deg_oligo[4:-1:1] == "NGCA"
                @test go[5:-1:1] == "TGC-A"  # Includes gap
            end

            @testset "Degenerate Base Edge Cases" begin
                # All degenerate bases
                all_deg = DegenerateOligo("NNNNNN", "all deg")
                @test n_deg_pos(all_deg) == 6
                @test n_unique_oligos(all_deg) == 4^6  # 4096
                @test length(collect(nondegens(all_deg))) == 4^6
                @test all(length(ol) == 6 for ol in Iterators.take(nondegens(all_deg), 10))

                # Mixed degenerate/non-degenerate
                mixed_deg = DegenerateOligo("ACGRYN", "mixed")
                @test n_deg_pos(mixed_deg) == 3  # R, Y, N
                @test n_unique_oligos(mixed_deg) == 2 * 2 * 4  # 16
                @test Set(String.(collect(nondegens(mixed_deg)))) ⊆ Set([String(['A','C','G',d,e,f]) for d in "AG" for e in "CT" for f in "ACGT"])
            end

            @testset "Error Handling for Invalid Inputs" begin
                # Invalid degenerate codes
                @test_throws ErrorException DegenerateOligo("ACGTZ")  # Z not in ALL_BASES
                @test_throws ErrorException DegenerateOligo("ACG#")  # Non-IUPAC character
                @test_throws ArgumentError GappedOligo(Oligo("ACGT"), [2=>-1])  # Negative gap length
            end

            @testset "SeqFold Interoperability" begin
                # revcomp for OligoView
                oligo = Oligo("ACGT", "revcomp test")
                ov = oligo[2:4]
                @test SeqFold.revcomp(ov) == "ACG"
                @test SeqFold.revcomp(ov) isa OligoView{Oligo}
                deg_oligo = DegenerateOligo("ACGN", "deg revcomp")
                deg_ov = deg_oligo[2:4]
                @test SeqFold.revcomp(deg_ov) == "NCG"
                @test SeqFold.revcomp(deg_ov) isa OligoView{DegenerateOligo}

                # complement for degenerate sequences
                @test SeqFold.complement(deg_oligo) == "TGCN"
                @test SeqFold.complement(deg_oligo) isa DegenerateOligo
                @test n_deg_pos(SeqFold.complement(deg_oligo)) == 1
                @test n_unique_oligos(SeqFold.complement(deg_oligo)) == 4
            end

            @testset "Long Sequences and Many Gaps" begin
                # Long sequence
                long_oligo = Oligo("A"^100, "long")
                @test length(long_oligo) == 100
                @test String(long_oligo[50:60]) == "A"^11
                @test SeqFold.gc_content(long_oligo) ≈ 0.0

                # Many gaps
                gaps = [i=>1 for i in 1:10:91]  # Gaps every 10 positions
                go_long = GappedOligo(long_oligo, gaps)
                @test length(go_long) == 110  # 100 bases + 10 gaps
                @test go_long[20:25] == "AAA-AA"
                @test n_deg_pos(go_long) == 0
                @test n_unique_oligos(go_long) == 1
            end

            @testset "Conversion Edge Cases" begin
                # OligoView to Oligo/DegenerateOligo
                ov = Oligo("ACGT")[2:3]
                @test convert(Oligo, ov) == Oligo("CG")
                deg_ov = DegenerateOligo("ACGN")[2:4]
                @test convert(DegenerateOligo, deg_ov) == DegenerateOligo("CGN")
                @test_throws InexactError convert(Oligo, DegenerateOligo("NN")[1:2])

                # GappedOligo to string and back
                go = GappedOligo(Oligo("ACGT"), [2=>1])
                go_str = String(go)
                @test GappedOligo(go_str) == go
            end
            
        end
    end
end

import DePPA.Oligos: NON_DEGEN_BASES, DEGEN_BASES, ALL_BASES, IUPAC_COUNTS

const NUM_RANDOM_TESTS = 10

rseq(len, bases) = join(rand(bases, len))
rdesc() = join(rand('1':'z', rand(0:150)))
rolg(::Type{Oligo}, len) =  Oligo(rseq(len, NON_DEGEN_BASES), rdesc())
rolg(::Type{DegenOligo}, len) =  DegenOligo(rseq(len, ALL_BASES), rdesc())
function rolg(::Type{GappedOligo}, len)
    len == 0 && return GappedOligo("", rdesc())
    seq_chars = collect(ALL_BASES)
    append!(seq_chars, fill('-', 5)) 
    seq = rseq(len, seq_chars)
    while count(==('-'), seq) == len
        seq = rseq(len, seq_chars)
    end
    return GappedOligo(seq, rdesc())
end
@testset "Types" begin
    @test Oligo <: AbstractOligo <: AbstractString
    @test DegenOligo <: AbstractDegen <: AbstractOligo
    @test GappedOligo <: AbstractGapped <: AbstractDegen <: AbstractOligo

    @test OligoView <: AbstractOligo
    @test OligoView{Union{Oligo, DegenOligo, GappedOligo}} <: AbstractOligo
    @test_throws TypeError OligoView{Int}
    
    @test NonDegenIterator{Union{Oligo, DegenOligo, GappedOligo}} <: NonDegenIterator
    @test_throws TypeError NonDegenIterator{Int}
end

@testset "Oligo Construction" begin
    @test Oligo() == Oligo("", "") == Oligo("")
    @test Oligo() !== Oligo("", "123")
    
    for _ in 1:NUM_RANDOM_TESTS
        len = rand(0:50)
        seq = rseq(len, NON_DEGEN_BASES)
        descr = rdesc()
        oligo = Oligo(seq, descr)
        @test oligo isa Oligo
        @test String(oligo) == uppercase(seq)
        @test description(oligo) == descr
        
        @test Oligo(DegenOligo(seq)) == oligo
        @test Oligo(GappedOligo(seq)) == oligo
    end
    
    for bad_char in setdiff('A':'Z', NON_DEGEN_BASES)
        @test_throws ErrorException Oligo("ACGT$bad_char")
    end
end

@testset "DegenOligo Construction" begin
    @test DegenOligo() == DegenOligo("")
    
    for _ in 1:NUM_RANDOM_TESTS
        len = rand(0:50)
        seq = rseq(len, ALL_BASES)
        descr = rdesc()
        deg_oligo = DegenOligo(seq, descr)
        @test deg_oligo isa DegenOligo
        @test String(deg_oligo) == uppercase(seq)
        @test description(deg_oligo) == descr
        
        if !any(c -> c in DEGEN_BASES, seq)
            nondeg = Oligo(seq, descr)
            @test DegenOligo(nondeg) == deg_oligo
        end
    end
    
    for bad_char in setdiff('A':'Z', ALL_BASES)
        @test_throws ErrorException DegenOligo("ACGTN$bad_char")
    end
end

@testset "GappedOligo Construction" begin
    @test GappedOligo() == GappedOligo("")
    
    for seq in ["-", "--A--", "A--", "--A", "A---T", "-G-C-", "ACGT", ""]
        go = GappedOligo(seq)
        @test String(go) == seq
        @test String(parent(go)) == filter(!=('-'), seq)
    end
    
    @test GappedOligo("A---T").gaps == [2 => 3]
    @test GappedOligo("-G-C-").gaps == [1 => 1, 2 => 1, 3 => 1]
    @test GappedOligo("ACGT").gaps == Pair{Int,Int}[]
    
    for _ in 1:NUM_RANDOM_TESTS
        len = rand(1:50)
        go = rolg(GappedOligo, len)
        seq_str = String(go)
        
        @test length(go) == len
        @test String(go) == seq_str
        @test String(parent(go)) == filter(!=('-'), seq_str)
        @test hasgaps(go) == any(==('-'), seq_str)
    end
end

@testset "Indexing and Views" begin
    for T in (Oligo, DegenOligo, GappedOligo)
        for _ in 1:NUM_RANDOM_TESTS
            len = rand(0:50)
            oligo = rolg(T, len)
            
            @test_throws BoundsError oligo[0]

            if len > 0
                idx = rand(1:len)
                @test oligo[idx] == String(oligo)[idx]
                @test oligo[idx] isa Char

                full_view = oligo[1:len]
                @test full_view isa OligoView{T}
                @test String(full_view) == String(oligo)
                @test length(full_view) == len
            end
            
            if len >= 2
                start = rand(1:len-1)
                stop = rand(start:len)
                view = oligo[start:stop]
                @test view isa OligoView
                @test String(view) == String(oligo)[start:stop]
                @test length(view) == stop - start + 1
                @test oligo_range(view) == start:stop
            end
            
            if len >= 4
                outer = oligo[2:len-1]
                inner = outer[1:length(outer)-1]
                @test String(inner) == String(oligo)[2:len-2]
            end
        end
    end
    
    gapped = GappedOligo("A-C-G-T")
    @test gapped[1] == 'A'
    @test gapped[2] == '-'
    @test gapped[3] == 'C'
    @test gapped[4] == '-'
    @test gapped[5] == 'G'
    @test gapped[6] == '-'
    @test gapped[7] == 'T'
    
    empty_oligo = Oligo("")
    @test_throws BoundsError empty_oligo[1]
    @test_throws BoundsError empty_oligo[1:1]
    non_empty = Oligo("A")
    @test_throws BoundsError non_empty[2]
    @test_throws BoundsError non_empty[0]
    @test_throws BoundsError non_empty[0:1]
end

@testset "Degeneracy Properties" begin
    nondeg = Oligo("ACGT")
    @test n_deg_pos(nondeg) == 0
    @test n_unique_oligos(nondeg) == 1
    @test !hasgaps(nondeg)
    
    deg = DegenOligo("ACGTNRYS")
    @test n_deg_pos(deg) == count(c -> c in DEGEN_BASES, "NRYS")
    @test n_unique_oligos(deg) == prod(IUPAC_COUNTS[c] for c in "NRYS")
    
    gapped_deg = GappedOligo("A-N-C-")
    @test n_deg_pos(gapped_deg) == 1
    @test n_unique_oligos(gapped_deg) == 4
    @test hasgaps(gapped_deg)
    
    view = deg[3:6]
    @test n_deg_pos(view) == count(c -> c in DEGEN_BASES, "GTNR")
    @test n_unique_oligos(view) == prod(IUPAC_COUNTS[c] for c in "NR")
end

@testset "Non-Degenerate Iteration" begin
    @test isempty(collect(nondegens(Oligo(""))))
    
    nondeg = Oligo("ACGT", "test")
    nondeg_iter = collect(nondegens(nondeg))
    @test only(nondeg_iter) == nondeg
    @test description(only(nondeg_iter)) == "test"
    
    deg = DegenOligo("ACN", "test")
    nondeg_iter = collect(nondegens(deg))
    @test length(nondeg_iter) == 4
    @test Set(String.(nondeg_iter)) == Set(["ACA", "ACC", "ACG", "ACT"])
    @test all(description(o) == "Non-degen sample from: test" for o in nondeg_iter)
    
    gapped_deg = GappedOligo("A-N-", "test")
    @test_throws ErrorException nondegens(gapped_deg)
    
    gapped_nogaps_deg = GappedOligo("AANC", "test")
    nondeg_iter = collect(nondegens(gapped_nogaps_deg))
    @test length(nondeg_iter) == 4
    expected = ["AAAC", "AAGC", "AACC", "AATC"]
    @test Set(String.(nondeg_iter)) == Set(expected)
    
    for _ in 1:NUM_RANDOM_TESTS
        len = rand(1:8)
        deg = rolg(DegenOligo, len)
        nd_iter = nondegens(deg)
        @test length(nd_iter) == n_unique_oligos(deg)
        @test length(collect(nd_iter)) == n_unique_oligos(deg)
    end
end

@testset "Sampling Functions" begin
    for T in (Oligo, DegenOligo, GappedOligo)
        for _ in 1:NUM_RANDOM_TESTS
            len = rand(1:50)
            oligo = rolg(T, len)
            sampled = sampleChar(oligo)
            @test sampled isa Char
            @test sampled in String(oligo)
        end
        @test_throws ArgumentError sampleChar(rolg(T, 0))
    end
    
    for T in (Oligo, DegenOligo, GappedOligo)
        for _ in 1:NUM_RANDOM_TESTS
            len = rand(5:50)
            oligo = rolg(T, len)
            view_len = rand(1:len)
            view = sampleView(oligo, view_len)
            @test view isa OligoView
            @test length(view) == view_len
            @test occursin(String(view), String(oligo))
        end
        empty_oligo = rolg(T, 0)
        @test_throws ArgumentError sampleView(empty_oligo, 1)
        non_empty = rolg(T, 3)
        @test_throws ArgumentError sampleView(non_empty, 4)
    end
    
    for _ in 1:NUM_RANDOM_TESTS
        nondeg = rolg(Oligo, rand(0:50))
        @test sampleNondeg(nondeg) == nondeg
        
        deg = rolg(DegenOligo, rand(1:50))
        sampled = sampleNondeg(deg)
        @test sampled isa DegenOligo
        @test length(sampled) == length(deg)
        @test all(c in NON_DEGEN_BASES for c in String(sampled))
        
        gapped_deg = rolg(GappedOligo, rand(1:50))
        sampled = sampleNondeg(gapped_deg)
        @test sampled isa GappedOligo
        @test length(sampled) == length(gapped_deg)
        @test all(c in NON_DEGEN_BASES || c == '-' for c in String(sampled))
        @test n_deg_pos(sampled) == 0
    end
end

@testset "Conversion and Promotion" begin
    @test convert(Oligo, DegenOligo("ACGT")) == Oligo("ACGT")
    @test convert(DegenOligo, Oligo("ACGT")) == DegenOligo("ACGT")
    @test convert(GappedOligo, Oligo("ACGT")) == GappedOligo("ACGT")
    
    @test promote_type(Oligo, DegenOligo) == DegenOligo
    @test promote_type(Oligo, GappedOligo) == GappedOligo
    @test promote_type(DegenOligo, GappedOligo) == GappedOligo
    @test promote_type(Oligo, String) == Oligo
    @test promote_type(DegenOligo, SubString) == DegenOligo
    
    @test Oligo("A") == "A"
    @test "A" == Oligo("A")
    @test DegenOligo("N") != Oligo("A")
    gapped = GappedOligo("A-C")
    @test gapped == GappedOligo("A-C")
    @test gapped != GappedOligo("AC-")
end

@testset "Iteration and Base Functions" begin
    for T in (Oligo, DegenOligo, GappedOligo)
        for _ in 1:NUM_RANDOM_TESTS
            len = rand(0:50)
            oligo = rolg(T, len)
            iterated = collect(oligo)
            @test length(iterated) == len
            @test String(oligo) == join(iterated)
        end
    end
    
    gapped = GappedOligo("--A--")
    @test collect(gapped) == ['-','-','A','-','-']
    
    for T in (Oligo, DegenOligo, GappedOligo)
        for _ in 1:NUM_RANDOM_TESTS
            len = rand(0:50)
            oligo = rolg(T, len)
            @test isempty(oligo) == (len == 0)
            @test length(oligo) == len
            @test lastindex(oligo) == len
            @test ncodeunits(oligo) == len
        end
    end
end
using Test
using DePPA.Alignments
using DePPA.Oligos
using DePPA.Primers

@testset "Evrogen Export" begin
    msa = MSA(["ACGTACGT", "ACGTACGT"])
    tm = (mean=55.0, conf=(53.0, 57.0), min=53.0, max=57.0)
    
    # Setup primers with different descriptions
    p1 = Primer{DegenOligo}(msa, 1:4, true, Oligo("ACGT", "MyForward"), 3, tm, -5.0, 0.5, 0.0)
    p2 = Primer{DegenOligo}(msa, 5:8, false, Oligo("TGCA", "MyReverse"), 3, tm, -5.0, 0.5, 0.0)
    p3 = Primer{DegenOligo}(msa, 1:4, true, Oligo("ACGT", ""), 3, tm, -5.0, 0.5, 0.0) # Empty description
    # Description with forbidden characters (should be replaced with spaces)
    p4 = Primer{DegenOligo}(msa, 5:8, false, Oligo("TGCA", "Bad;Desc\t"), 3, tm, -5.0, 0.5, 0.0)
    
    pair1 = p1 => p2
    
    # --- IO Export Tests ---
    @testset "IO Export" begin
        # Single primer
        io = IOBuffer()
        export_evrogen(io, p1)
        content = String(take!(io))
        @test content == "MyForward_F_1; ACGT; 0.04\n"
        
        # Vector of primers with custom scale (Float)
        io = IOBuffer()
        export_evrogen(io, [p1, p2]; scale=1.0)
        content = String(take!(io))
        @test content == "MyForward_F_1; ACGT; 1.0\nMyReverse_R_8; TGCA; 1.0\n"
        
        # Single pair with scale as a string
        io = IOBuffer()
        export_evrogen(io, pair1; scale="0.2")
        content = String(take!(io))
        @test content == "MyForward_F_1; ACGT; 0.2\nMyReverse_R_8; TGCA; 0.2\n"
        
        # Vector of pairs
        io = IOBuffer()
        export_evrogen(io, [pair1, pair1])
        content = String(take!(io))
        @test occursin("MyForward_F_1; ACGT; 0.04", content)
        @test occursin("MyReverse_R_8; TGCA; 0.04", content)
        @test length(split(content, '\n')) == 5 # 4 lines + empty trailing line
        
        # Check name generation: empty description and sanitization
        io = IOBuffer()
        export_evrogen(io, [p3, p4])
        content = String(take!(io))
        # p3: desc empty -> DePPA_F_1_1
        @test occursin("DePPA_F_1_1; ACGT; 0.04", content)
        # p4: desc "Bad;Desc\t" -> "Bad Desc ", pos=5, idx=2
        @test occursin("Bad Desc _R_8; TGCA; 0.04", content)
    end
    
    # --- File Export Tests ---
    @testset "File Export" begin
        tmpfile = tempname() * ".txt"
        try
            # Vector of primers
            ret = export_evrogen(tmpfile, [p1, p2])
            @test ret == tmpfile
            lines = readlines(tmpfile)
            @test length(lines) == 2
            @test lines[1] == "MyForward_F_1; ACGT; 0.04"
            @test lines[2] == "MyReverse_R_8; TGCA; 0.04"
            
            # Single primer
            export_evrogen(tmpfile, p1; scale=0.2)
            lines = readlines(tmpfile)
            @test length(lines) == 1
            @test lines[1] == "MyForward_F_1; ACGT; 0.2"
            
            # Vector of pairs
            export_evrogen(tmpfile, [pair1])
            lines = readlines(tmpfile)
            @test length(lines) == 2
            @test lines[1] == "MyForward_F_1; ACGT; 0.04"
            @test lines[2] == "MyReverse_R_8; TGCA; 0.04"
            
            # Single pair
            export_evrogen(tmpfile, pair1)
            lines = readlines(tmpfile)
            @test length(lines) == 2
            @test lines[1] == "MyForward_F_1; ACGT; 0.04"
            @test lines[2] == "MyReverse_R_8; TGCA; 0.04"
        finally
            rm(tmpfile, force=true)
        end
    end
end
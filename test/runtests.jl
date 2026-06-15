using Test
using Aqua
using JET
using Random

using DEPPA
using DEPPA.Oligos
using DEPPA.Alignments
using DEPPA.Primers

Random.seed!(42)

@testset verbose=true failfast=true "DEPPA.jl"  begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(DEPPA)
    end
    @testset "Code linting (JET.jl)" begin
        JET.test_package(DEPPA; target_modules=(DEPPA,))
    end
    
    @testset "Oligos" include("test_oligos.jl")
    @testset "Alignments" include("test_alignments.jl")

    @testset "Primers" include("test_primers.jl")
    # @testset "SeqFold methods" include("test_seqfold.jl")


end
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
    # Passing
    @testset "Code quality (Aqua.jl)" Aqua.test_all(DEPPA)
    @testset "Code linting (JET.jl)" JET.test_package(DEPPA; target_modules=(DEPPA,))
    @testset "Oligos" include("test_oligos.jl")
    @testset "Alignments" include("test_alignments.jl")
    @testset "Primers" include("test_primers.jl")

    # New
    @testset "Real pipeline" include("test_pipeline.jl")


    # @testset "SeqFold methods" include("test_seqfold.jl")


end
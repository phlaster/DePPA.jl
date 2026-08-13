using Test
using Aqua
using JET
using Random

using DePPA
using DePPA.Oligos
using DePPA.Alignments
using DePPA.Primers

Random.seed!(42)

@testset verbose = true failfast = true "DePPA.jl" begin
    @testset "Code quality (Aqua.jl)" Aqua.test_all(DePPA)
    @testset "Code linting (JET.jl)" JET.test_package(DePPA; target_modules = (DePPA,))
    @testset "Oligos" include("test_oligos.jl")
    @testset "Alignments" include("test_alignments.jl")
    @testset "Primers" include("test_primers.jl")
    @testset "Real pipeline" include("test_pipeline.jl")
    @testset "Show methods" include("test_show.jl")
    @testset "SeqFold methods" include("test_seqfold.jl")
    @testset "Export Primers" include("test_export.jl")
end

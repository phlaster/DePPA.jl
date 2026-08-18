using DePPA.Alignments, DePPA.Primers, MAFFT_jll, SeqFold, Random

Random.seed!(123)

file = "assets/TruA.fasta.gz"

aln = MSA(file; mafft = true, bootstrap = 10);
primers = construct_primers(aln)
ppairs = best_pairs(primers; amplicon_len = 190:190, sortby=:tm)

tempfile = tempname()
reffile = "assets/TruA_primers.txt"

open(tempfile, "w") do f
    for pp in ppairs
        show(f, MIME"text/plain"(), pp)
        println(f)
    end
end

@test readlines(tempfile) == readlines(reffile)

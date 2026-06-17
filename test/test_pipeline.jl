using DePPA, MAFFT_jll, SeqFold

Random.seed!(123)

file = "assets/TruA.fasta.gz"

aln = MSA(file; mafft=true, bootstrap=10);
fwds = construct_primers(aln)
revs = construct_primers(aln; is_forward=false)

ppairs = best_pairs(fwds, revs; amplicon_len=190:190)

tempfile = tempname()
reffile = "assets/TruA_primers.txt"

open(tempfile, "w") do f
    for pp in ppairs
        show(f, MIME"text/plain"(), pp)
        println(f)
    end
end

@test readlines(tempfile) == readlines(reffile)

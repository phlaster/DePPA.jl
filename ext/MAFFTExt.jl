module MAFFTExt

using MAFFT_jll
using FastaIO

import DePPA.Alignments: _align!

function _align!(fasta_content::Vector{NTuple{2, String}})
    buffer = IOBuffer()
    writefasta(buffer, fasta_content; check_description=false)
    frombuffer = take!(buffer)

    proc = open(`$(MAFFT_jll.mafft()) --quiet -`, "r+")
    # proc = open(pipeline(`$(MAFFT_jll.mafft()) -`, stderr=devnull), "r+")
    write(proc, frombuffer)
    close(proc.in)
    aligned_output = read(proc, String)

    fasta_content .= readfasta(IOBuffer(aligned_output))
    return
end

end  # module

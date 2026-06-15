function _show_primer_common(io::IO, primer::AbstractPrimer)
    println(io, "  Sequence: ", primer.consensus)
    println(io, "  Length: ", length(primer.consensus))
    println(io, "  Positions: ", primer.pos)
    println(io, "  Unique variants: ", n_unique_oligos(primer))
    println(io, "  Melting temperature: ", round(primer.tm.mean, digits=1), "°C (",
            round(primer.tm.conf[1], digits=1), "⋅",
            round(primer.tm.conf[2], digits=1), "°C)")
    println(io, "  Min ΔG: ", round(primer.dg, digits=2), " kcal/mol")
    println(io, "  GC content: ", round(primer.gc * 100, digits=1), "%")
    print(io, "  Description: \"", description(primer), "\"")
end

function Base.show(io::IO, primer::AbstractPrimer)
    max_width = 15
    seq = primer.consensus
    seq_display = length(seq) > max_width ? seq[1:max_width-3] * "..." : seq

    dir_str = primer.is_forward ? "forward" : "reverse"
    print(io, "Primer(\"", seq_display, "\", len=$(length(primer.consensus)), ",
          "pos=$(primer.pos.start):$(primer.pos.stop), $dir_str")

    print(io, ", degen=$(n_deg_pos(primer)), variants=$(n_unique_oligos(primer) > 10000 ? ">10k" : n_unique_oligos(primer)))")
    print(io, ", Tm=$(round(primer.tm.mean, digits=1))°C, ",
          "ΔG=$(round(primer.dg, digits=2))kcal/mol, ",
          "GC=$(round(primer.gc * 100, digits=1))%)")
end

function Base.show(io::IO, ::MIME"text/plain", primer::AbstractPrimer)
    dir_str = primer.is_forward ? "Forward" : "Reverse"
    ndeg = n_deg_pos(primer)
    deg_status = ndeg > 0 ?
        "degenerate primer with $ndeg deg. positions" :
        "non-degenerate primer"
    println(io, "$dir_str $deg_status")

    L = length(primer.msa)
    s, _..., e = primer.pos
    term_width = displaysize(io)[2] - 1
    bar_start = "|"
    bar_end = string(L, '|')
    eq_len = max(0, term_width - length(bar_start) - length(bar_end))
    bar = bar_start * repeat('=', eq_len) * bar_end

    if L == 1
        map_start_col = length(bar_start) + 1
        map_end_col = map_start_col + eq_len - 1
        col_s = col_e = map_start_col
    else
        map_start_col = length(bar_start) + 1
        map_end_col = map_start_col + eq_len - 1
        scale = (eq_len - 1) / (L - 1.0)
        col_s = map_start_col + round(Int, (s - 1) * scale)
        col_e = map_start_col + round(Int, (e - 1) * scale)
    end

    three_prime_col, label = primer.is_forward ?
        (col_e, string('\\', s, ':', e, '>')) :
        (col_s, string('<', s, ':', e, '\\'))

    label_len = length(label)
    indent = max(0, three_prime_col - label_len)
    label_line = " "^indent * label

    if primer.is_forward
        println(io, label_line)
        println(io, bar)
        println(io)
    else
        println(io)
        println(io, bar)
        println(io, label_line)
    end

    _show_primer_common(io, primer)
end

function Base.show(io::IO, ::MIME"text/plain", pp::Pair{<:AbstractPrimer})
    fwd, rev = pp.first, pp.second

    # Validate that this is a proper forward/reverse pair
    if !fwd.is_forward || rev.is_forward
        invoke(show, Tuple{IO, MIME"text/plain", Pair}, io, MIME"text/plain"(), pp)
        return
    end

    try
        msa = fwd.msa
        if rev.msa !== msa
            invoke(show, Tuple{IO, MIME"text/plain", Pair}, io, MIME"text/plain"(), pp)
            return
        end

        N = nseqs(msa)
        L = length(msa)
        amp_start = fwd.pos.start
        amp_end = rev.pos.stop
        amp_len = amp_end - amp_start + 1
        overlap = fwd.pos.stop >= rev.pos.start ? "!!! OVERLAPPING !!! " : ""
        header = "$(overlap)PCR primer pair for $N seq. MSA, amplicon: $amp_start:$amp_end ($(amp_len)bp)"
        println(io, header)

        term_width = displaysize(io)[2] - 1
        bar_start = "|"
        bar_end = string(L, '|')
        eq_len = max(0, term_width - length(bar_start) - length(bar_end))
        bar = bar_start * repeat('=', eq_len) * bar_end
        scale = L == 1 ? 0.0 : (eq_len - 1) / (L - 1.0)
        map_start_col = length(bar_start) + 1
        col_amp_start = map_start_col + (L > 1 ? round(Int, (amp_start - 1) * scale) : 0)
        col_amp_end = map_start_col + (L > 1 ? round(Int, (amp_end - 1) * scale) : 0)
        arrow_line = [' ' for _ in 1:term_width]
        if col_amp_start <= term_width && col_amp_start >= 1
            arrow_line[col_amp_start] = '>'
        end
        if col_amp_end <= term_width && col_amp_end >= 1 && col_amp_end != col_amp_start
            arrow_line[col_amp_end] = '<'
        end
        for i in (col_amp_start + 1):(col_amp_end - 1)
            if 1 <= i <= term_width
                arrow_line[i] = '_'
            end
        end
        label = fwd.pos.stop >= rev.pos.start ? "" : "$(amp_len)bp"
        label_len = length(label)
        inner_len = col_amp_end - col_amp_start - 1
        if inner_len >= label_len + 2
            mid_col = col_amp_start + div(col_amp_end - col_amp_start + 1, 2)
            label_start = mid_col - div(label_len, 2)
            if label_start >= col_amp_start + 2 && label_start + label_len - 1 <= col_amp_end - 2
                for (j, c) in enumerate(label)
                    pos = label_start + j - 1
                    if 1 <= pos <= term_width
                        arrow_line[pos] = c
                    end
                end
            end
        end
        println(io, join(arrow_line))
        println(io, bar)
        println(io, "Forward: ", fwd.consensus, " at ", fwd.pos.start, ":", fwd.pos.stop)
        println(io, "Reverse: ", rev.consensus, " at ", rev.pos.start, ":", rev.pos.stop)
        mean_tm = (fwd.tm.mean + rev.tm.mean) / 2
        delta_tm = abs(fwd.tm.mean - rev.tm.mean) / 2
        print(io, "Tm: ", round(mean_tm, digits=1), "±", round(delta_tm, digits=1), " °C")
    catch e
        invoke(show, Tuple{IO, MIME"text/plain", Pair}, io, MIME"text/plain"(), pp)
    end
end

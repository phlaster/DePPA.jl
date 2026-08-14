# fatou-ignore-file undefined-name

function _show_primer_common(io::IO, primer::AbstractPrimer)
    println(io, "  Sequence: ", String(primer))
    println(io, "  Length: ", length(primer))
    println(io, "  Positions: ", primer.pos)
    println(io, "  Unique variants: ", n_unique_oligos(primer))
    println(
        io,
        "  Melting temperature: ",
        round(primer.tm.mean, digits = 1),
        "°C (",
        round(primer.tm.conf[1], digits = 1),
        "⋅",
        round(primer.tm.conf[2], digits = 1),
        "°C)",
    )
    println(io, "  Min ΔG: ", round(primer.dg, digits = 2), " kcal/mol")
    println(io, "  GC content: ", round(primer.gc * 100, digits = 1), "%")
    print(io, "  Description: \"", description(primer), "\"")
end

function Base.show(io::IO, primer::AbstractPrimer)
    max_width = 15
    seq = String(primer)
    seq_display = length(seq) > max_width ? seq[1:max_width - 3] * "..." : seq

    dir_str = primer.is_forward ? "forward" : "reverse"
    print(
        io,
        "Primer(\"",
        seq_display,
        "\", len=$(length(primer.consensus)), ",
        "pos=$(primer.pos.start):$(primer.pos.stop), $dir_str",
    )

    print(
        io,
        ", degen=$(n_deg_pos(primer)), variants=$(n_unique_oligos(primer) > 10000 ? ">10k" : n_unique_oligos(primer)))",
    )
    print(
        io,
        ", Tm=$(round(primer.tm.mean, digits = 1))°C, ",
        "ΔG=$(round(primer.dg, digits = 2))kcal/mol, ",
        "GC=$(round(primer.gc * 100, digits = 1))%)",
    )
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
        col_s = col_e = map_start_col
    else
        map_start_col = length(bar_start) + 1
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
            if label_start >= col_amp_start + 2 &&
                    label_start + label_len - 1 <= col_amp_end - 2
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
        f_desc = description(fwd)
        r_desc = description(rev)

        f_str = isempty(f_desc) ?
            "" :
            " (" * (length(f_desc) > 12 ? first(f_desc, 12) * "…" : f_desc) * ")"
        r_str = isempty(r_desc) ?
            "" :
            " (" * (length(r_desc) > 12 ? first(r_desc, 12) * "…" : r_desc) * ")"

        println(
            io,
            "Forward: ",
            String(fwd),
            " at ",
            fwd.pos.start,
            ":",
            fwd.pos.stop,
            f_str,
        )
        println(
            io,
            "Reverse: ",
            String(rev),
            " at ",
            rev.pos.start,
            ":",
            rev.pos.stop,
            r_str,
        )
        mean_tm = (fwd.tm.mean + rev.tm.mean) / 2
        delta_tm = abs(fwd.tm.mean - rev.tm.mean) / 2
        print(
            io,
            "Tm: ",
            round(mean_tm, digits = 1),
            "±",
            round(delta_tm, digits = 1),
            " °C",
        )
    catch e
        invoke(show, Tuple{IO, MIME"text/plain", Pair}, io, MIME"text/plain"(), pp)
    end
end

function Base.show(io::IO, ::MIME"text/plain", primers::AbstractVector{<:AbstractPrimer})
    isempty(primers) && (print(io, "Empty primer vector"); return)
    
    anymsa = root(first(primers).msa)
    
    if !all(root(p.msa) == anymsa for p in primers)
        println(io, "Vector of $(length(primers)) primers (different MSAs):")
        for p in primers
            print(io, "  ")
            show(io, p)
            println(io)
        end
        return
    end
        
    L = length(anymsa)
    fwd_counts = zeros(Int, L)
    rev_counts = zeros(Int, L)
    
    for p in primers
        if p.is_forward
            fwd_counts[p.pos.start:p.pos.stop] .+= 1
        else
            rev_counts[p.pos.start:p.pos.stop] .+= 1
        end
    end
    
    term_width = displaysize(io)[2] - 1
    max_val = max(maximum(fwd_counts), maximum(rev_counts), 1)
    max_val_str = string(max_val)
    
    scale_w = max(length(max_val_str), 3) + 3
    W = max(10, term_width - scale_w)
    
    fwd_bins = zeros(Int, W)
    rev_bins = zeros(Int, W)
    scale = L == 1 ? 0.0 : (W - 1) / (L - 1)
    
    for pos in 1:L
        col = floor(Int, (pos - 1) * scale) + 1
        fwd_bins[col] += fwd_counts[pos]
        rev_bins[col] += rev_counts[pos]
    end
    
    H = 4
    
    function get_lower_block(frac)
        sub = floor(Int, frac * 8)
        return " ▁▂▃▄▅▆▇█"[sub + 1]
    end
    function get_upper_block(frac)
        if frac >= 0.5
            return '▀'
        elseif frac > 0.0
            return '▔'
        else
            return ' '
        end
    end

    println(io, "Primer distribution for $(nseqs(anymsa)) seq MSA (L=$L): $(length(primers)) primers")
    
    for i in 1:H
        level = H - i + 1
        frac_threshold = level / H
        next_frac_threshold = (level - 1) / H
        
        if i == 1
            print(io, rpad(string(max_val), scale_w - 3), " ┤ ")
        else
            print(io, rpad("", scale_w - 3), " │ ")
        end
        
        for col in 1:W
            val = fwd_bins[col] / max_val
            if val >= frac_threshold
                print(io, '█')
            elseif val > next_frac_threshold
                print(io, get_lower_block((val - next_frac_threshold) / (frac_threshold - next_frac_threshold)))
            else
                print(io, ' ')
            end
        end
        println(io)
    end
    
    print(io, rpad("0", scale_w - 3), " ┼ ═")
    print(io, repeat('═', W))
    println(io)
    
    for i in 1:H
        level = H - i + 1
        frac_threshold = level / H
        next_frac_threshold = (level - 1) / H
        
        if i == H
            print(io, rpad(string(max_val), scale_w - 3), " ┘ ")
        else
            print(io, rpad("", scale_w - 3), " │ ")
        end
        
        for col in 1:W
            val = rev_bins[col] / max_val
            if val >= frac_threshold
                print(io, '█')
            elseif val > next_frac_threshold
                print(io, get_upper_block((val - next_frac_threshold) / (frac_threshold - next_frac_threshold)))
            else
                print(io, ' ')
            end
        end
        println(io)
    end
    
    print(io, rpad("", scale_w - 3), "   ")
    str_L = string(L)
    str_1 = "1"
    if W >= length(str_L) + length(str_1) + 2
        print(io, str_1)
        print(io, repeat(' ', W - length(str_L) - length(str_1)))
        print(io, str_L)
    else
        print(io, repeat(' ', W))
    end
    println(io)
end
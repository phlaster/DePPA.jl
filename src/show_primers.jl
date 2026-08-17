# fatou-ignore-file undefined-name

function _truncate_seq(s::AbstractString, w::Integer)
    length(s) <= w && return s
    w <= 0 && return ""
    return first(s, w - 1) * "…"
end

function _primer_seq_str(primer::AbstractPrimer, max_width::Integer=typemax(Int) >> 2)
    main = String(primer.consensus)
    if isnothing(primer.adapter) || isempty(primer.adapter)
        return _truncate_seq(main, max_width)
    end
    adapter = String(primer.adapter)

    full = "[$adapter]$main"
    length(full) <= max_width && return full

    if max_width >= length(main) + 4
        keep = max_width - length(main) - 3
        return "[" * first(adapter, keep) * "…]" * main
    end

    keep_main = max(max_width - 4, 1)
    return "[" * first(adapter, 1) * "…]" * _truncate_seq(main, keep_main)
end

function _show_primer_common(io::IO, primer::AbstractPrimer)
    seq_label = "  Sequence: "
    seq_w = displaysize(io)[2] - length(seq_label)
    println(io, seq_label, _primer_seq_str(primer, seq_w))
    println(io, "  Length: ", length(primer))
    println(io, "  Positions: ", primer.pos)
    println(io, "  Unique variants: ", n_unique_oligos(primer))
    println(
        io,
        "  Melting temperature: ",
        round(primer.tm.mean, digits=1),
        "°C (",
        round(primer.tm.conf[1], digits=1),
        "⋅",
        round(primer.tm.conf[2], digits=1),
        "°C)",
    )
    println(io, "  Min ΔG: ", round(primer.dg, digits=2), " kcal/mol")
    println(io, "  GC content: ", round(primer.gc * 100, digits=1), "%")
    print(io, "  Description: \"", description(primer), "\"")
end

function Base.show(io::IO, primer::AbstractPrimer)
    max_width = 15
    seq_display = _primer_seq_str(primer, max_width)

    dir_str = primer.is_forward ? "forward" : "reverse"
    print(
        io,
        "Primer(\"",
        seq_display,
        "\", len=$(length(primer.consensus)), ",
        "pos=$(first(primer.pos)):$(last(primer.pos)), $dir_str",
    )

    print(
        io,
        ", degen=$(n_deg_pos(primer)), variants=$(n_unique_oligos(primer) > 10000 ? ">10k" : n_unique_oligos(primer)), ",
    )
    print(
        io,
        "Tm=$(round(primer.tm.mean, digits = 1))°C, ",
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
    s, e = first(primer.pos), last(primer.pos)
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

function Base.show(io::IO, ::MIME"text/plain", pp::Pair{<:AbstractPrimer,<:AbstractPrimer})
    fwd, rev = pp.first, pp.second

    if !fwd.is_forward || rev.is_forward
        invoke(show, Tuple{IO,MIME"text/plain",Pair}, io, MIME"text/plain"(), pp)
        return
    end

    try
        msa = fwd.msa
        if rev.msa !== msa
            invoke(show, Tuple{IO,MIME"text/plain",Pair}, io, MIME"text/plain"(), pp)
            return
        end

        N = nseqs(msa)
        L = length(msa)
        amp_start = first(fwd.pos)
        amp_end = last(rev.pos)
        amp_len = amp_end - amp_start + 1
        overlap = last(fwd.pos) >= first(rev.pos) ? "!!! OVERLAPPING !!! " : ""
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
        for i in (col_amp_start+1):(col_amp_end-1)
            if 1 <= i <= term_width
                arrow_line[i] = '_'
            end
        end
        label = last(fwd.pos) >= first(rev.pos) ? "" : "$(amp_len)bp"
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

        f_tail = " at $(first(fwd.pos)):$(last(fwd.pos))" * f_str
        r_tail = " at $(first(rev.pos)):$(last(rev.pos))" * r_str
        f_seq = _primer_seq_str(fwd, max(5, term_width - length("Forward: ") - length(f_tail)))
        r_seq = _primer_seq_str(rev, max(5, term_width - length("Reverse: ") - length(r_tail)))
        println(io, "Forward: ", f_seq, f_tail)
        println(io, "Reverse: ", r_seq, r_tail)

        mean_tm = (fwd.tm.mean + rev.tm.mean) / 2
        delta_tm = abs(fwd.tm.mean - rev.tm.mean) / 2
        print(
            io,
            "Tm: ",
            round(mean_tm, digits=1),
            "±",
            round(delta_tm, digits=1),
            " °C",
        )
    catch e
        invoke(show, Tuple{IO,MIME"text/plain",Pair}, io, MIME"text/plain"(), pp)
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
            fwd_counts[first(p.pos):last(p.pos)] .+= 1
        else
            rev_counts[first(p.pos):last(p.pos)] .+= 1
        end
    end

    terminal_height, terminal_width = displaysize(io)
    H = clamp(div(terminal_height - 6, 2), 3, 12)

    max_val = max(maximum(fwd_counts), maximum(rev_counts), 1)
    scale_w = max(length(string(max_val)), 3) + 3
    W = max(10, terminal_width - 1 - scale_w)

    fwd_bins = zeros(Int, W)
    rev_bins = zeros(Int, W)
    scale = L == 1 ? 0.0 : (W - 1) / (L - 1)
    col_of(pos) = floor(Int, (pos - 1) * scale) + 1
    for pos in 1:L
        c = col_of(pos)
        fwd_bins[c] = max(fwd_bins[c], fwd_counts[pos])
        rev_bins[c] = max(rev_bins[c], rev_counts[pos])
    end
    max_val = max(maximum(fwd_bins), maximum(rev_bins), 1)

    lower_blocks = (' ', '▁', '▂', '▃', '▄', '▅', '▆', '▇', '█')
    get_lower_block(frac) = lower_blocks[clamp(floor(Int, frac*8), 0, 7)+1]

    function get_upper_block(frac)
        frac >= 0.5 && return '▀'
        frac > 0.0 && return '▔'
        return ' '
    end

    axis = fill('═', W)
    for mark in 100:100:L
        axis[col_of(mark)] = '╪'
    end
    axis[1] = '╞'

    nearest_hundred = round(Int, L / 100) * 100
    cell_width = (L - 1) / (W - 1)
    if nearest_hundred > 0 && abs(L - nearest_hundred) <= 2 * cell_width
        axis[W] = '╡'
    end

    prefix(label, tick) = rpad(label, scale_w - 3) * tick
    show_half = iseven(H) && iseven(max_val)
    half_label = string(div(max_val, 2))
    half_row_fwd = div(H, 2) + 1
    half_row_rev = div(H, 2)

    println(io, "Primer distribution for $(nseqs(anymsa)) seq MSA (L=$L): $(length(primers)) primers")

    for i in 1:H
        if i == 1
            print(io, prefix(string(max_val), " ┤ "))
        elseif show_half && i == half_row_fwd
            print(io, prefix(half_label, " ┤ "))
        else
            print(io, prefix("", " │ "))
        end
        for col in 1:W
            h = fwd_bins[col] / max_val * H
            lo, hi = H - i, H - i + 1
            if h >= hi
                print(io, '█')
            elseif h > lo
                print(io, get_lower_block(h - lo))
            else
                print(io, ' ')
            end
        end
        println(io)
    end

    print(io, prefix("0", " ┼ "))
    print(io, String(axis))
    println(io)

    for i in 1:H
        if i == H
            print(io, prefix(string(max_val), " ┤ "))
        elseif show_half && i == half_row_rev
            print(io, prefix(half_label, " ┤ "))
        else
            print(io, prefix("", " │ "))
        end
        for col in 1:W
            h = rev_bins[col] / max_val * H
            lo, hi = i - 1, i
            if h >= hi
                print(io, '█')
            elseif h > lo
                print(io, get_upper_block(h - lo))
            else
                print(io, ' ')
            end
        end
        println(io)
    end

    print(io, prefix("", "   "))
    str_L = string(L)
    if W >= length(str_L) + 2
        print(io, '1', repeat(' ', W - length(str_L) - 1), str_L)
    else
        print(io, repeat(' ', W))
    end

    return nothing
end

function Base.show(io::IO, ::MIME"text/plain", pairs::AbstractVector{<:Pair{<:AbstractPrimer, <:AbstractPrimer}})
    isempty(pairs) && (print(io, "Empty primer pair vector"); return)

    npairs = length(pairs)
    anymsa = root(first(pairs).first.msa)

    same_msa = all(
        p.first.is_forward && !p.second.is_forward &&
            root(p.first.msa) == anymsa && root(p.second.msa) == anymsa
        for p in pairs
    )

    if !same_msa
        println(io, "Vector of $npairs primer pairs:")
        for p in pairs
            print(io, "  ")
            show(io, p)
            println(io)
        end
        return
    end

    L = length(anymsa)
    terminal_height, terminal_width = displaysize(io)

    idx_w = length(string(npairs))
    lab_w = idx_w + 2

    max_pair_lines = max(4, terminal_height - 5)
    if npairs <= max_pair_lines
        top_idx = 1:npairs
        bottom_idx = 1:0
        hidden = 0
    else
        n_half = div(max_pair_lines - 1, 2)
        top_idx = 1:n_half
        bottom_idx = (npairs - n_half + 1):npairs
        hidden = npairs - 2 * n_half
    end

    function tm_str(i)
        fwd, rev = pairs[i]
        mean_tm = (fwd.tm.mean + rev.tm.mean) / 2
        delta_tm = abs(fwd.tm.mean - rev.tm.mean) / 2
        return string(round(mean_tm, digits = 1), "±", round(delta_tm, digits = 1))
    end

    tm_w = length("Tm °C")
    for i in Iterators.flatten((top_idx, bottom_idx))
        tm_w = max(tm_w, length(tm_str(i)))
    end

    map_w = max(10, terminal_width - lab_w - 1 - 4 - tm_w)

    bar_start = "|"
    bar_end = string(L, '|')
    eq_len = max(0, map_w - length(bar_start) - length(bar_end))
    bar = bar_start * repeat('=', eq_len) * bar_end
    scale = L == 1 ? 0.0 : (eq_len - 1) / (L - 1.0)

    function pair_line(i)
        fwd, rev = pairs[i]
        amp_start = first(fwd.pos)
        amp_end = last(rev.pos)
        amp_len = amp_end - amp_start + 1

        chars = fill(' ', map_w)
        col_s = 2 + (L > 1 ? round(Int, (amp_start - 1) * scale) : 0)
        col_e = 2 + (L > 1 ? round(Int, (amp_end - 1) * scale) : 0)
        col_s = clamp(col_s, 1, map_w)
        col_e = clamp(col_e, 1, map_w)

        chars[col_s] = '>'
        if col_e != col_s
            chars[col_e] = '<'
            for c in (col_s + 1):(col_e - 1)
                chars[c] = '-'
            end
            label = "$(amp_len)bp"
            inner_len = col_e - col_s - 1
            if inner_len >= length(label) + 2
                label_start = col_s + 1 + div(inner_len - length(label), 2)
                for (j, ch) in enumerate(label)
                    chars[label_start + j - 1] = ch
                end
            end
        end

        return "[" * lpad(i, idx_w) * "] " * String(chars) * lpad(tm_str(i), tm_w)
    end

    lines = String[pair_line(i) for i in top_idx]
    if hidden > 0
        msg = "... $hidden more ..."
        mchars = fill(' ', map_w)
        start = max(1, 1 + div(map_w - length(msg), 2))
        for (j, ch) in enumerate(msg)
            c = start + j - 1
            c <= map_w && (mchars[c] = ch)
        end
        push!(lines, " "^(lab_w + 1) * String(mchars) * lpad("...", tm_w-2))
    end
    append!(lines, (pair_line(i) for i in bottom_idx))

    println(io, "$npairs PCR primer pairs for $(nseqs(anymsa)) seq. MSA:")
    println(io, " "^(lab_w + 1), bar, lpad("Tm °C", tm_w))
    print(io, join(lines, '\n'))

    return nothing
end
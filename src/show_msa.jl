const BASE_COLORS = Dict{Char, Symbol}(
    'A' => :green,
    'C' => :blue,
    'G' => :yellow,
    'T' => :red,
    'M' => :light_cyan,   # A/C
    'R' => :light_green,  # A/G
    'W' => :light_green,  # A/T
    'S' => :light_blue,   # C/G
    'Y' => :light_magenta,# C/T
    'K' => :light_red,    # G/T
    'V' => :cyan,         # A/C/G
    'H' => :magenta,      # A/C/T
    'D' => :green,        # A/G/T
    'B' => :blue,         # C/G/T
    'N' => :white,
    '-' => :normal
)

const histbars = ['▁', '▂', '▃', '▄', '▅', '▆', '▇', '█']

const MSA_SHOW_STYLE = Ref(:polymorf)

"""
    setMSAShowStyle!(style::Symbol)

Sets the global display style for the MSA viewer. 

Valid options are:
- `:bw` (black and white, consensus with `.` for matches)
- `:polymorf` (colored polymorphic characters, `.` for consensus matches)
- `:allcolors` (fully colored sequences and depth histogram bars, default)
"""
function setMSAShowStyle!(style::Symbol)
    if style in (:bw, :polymorf, :allcolors)
        MSA_SHOW_STYLE[] = style
    else
        throw(ArgumentError(
            "Invalid MSA show style: :$style. " *
            "Valid options are :bw, :polymorf, or :allcolors."
        ))
    end
end

function Base.show(io::IO, msa::AbstractMSA)
    n_sequences = nseqs(msa)
    if n_sequences == 0
        print(io, "Empty ", typeof(msa))
        return
    end
    seq_length = length(msa)
    print(io, typeof(msa), " with $n_sequences sequences of length $seq_length")
    seq_length == 0 && return

    terminal_height, terminal_width = displaysize(io)

    max_display_height = max(5, terminal_height - 9)
    n_display_seqs = min(n_sequences, max_display_height)

    processed_descs = Vector{String}(undef, n_display_seqs)
    has_desc = false
    for i in 1:n_display_seqs
        desc = description(getsequence(msa, i))
        processed = replace(desc, '\n' => ' ', '\t' => ' ')
        processed_descs[i] = processed
        if !isempty(processed)
            has_desc = true
        end
    end

    desc_width = 0
    padded_descs = String[]
    if has_desc
        max_possible_desc = min(50, floor(Int, 0.2 * terminal_width))
        max_desc_len = maximum(length, processed_descs)
        desc_width = min(max_possible_desc, max_desc_len)
        padded_descs = Vector{String}(undef, n_display_seqs)
        for i in 1:n_display_seqs
            trunc = processed_descs[i][1:min(end, desc_width)]
            padded_descs[i] = rpad(trunc, desc_width)
        end
        desc_width += 2  # for separator
    end

    max_display_width = max(20, terminal_width - desc_width - 7)

    needs_width_ellipsis = seq_length > max_display_width
    max_seq_chars = needs_width_ellipsis ? max_display_width - 3 : max_display_width
    displayed_cols_range = 1:min(seq_length, max_seq_chars)

    abs_cols = if msa isa MSAView
        msa.cols.start .+ (displayed_cols_range .- 1)
    else
        displayed_cols_range
    end

    println(io, ":")
    has_desc && print(io, " "^desc_width)
    
    style = MSA_SHOW_STYLE[]
    
    # Precompute consensus characters for :bw and :polymorf styles
    consensus_chars = Vector{Char}(undef, length(displayed_cols_range))
    if style != :allcolors
        for dj in 1:length(displayed_cols_range)
            counts = Dict{Char, Int}()
            for i in 1:n_sequences
                c = getsequence(msa, i, dj)
                counts[c] = get(counts, c, 0) + 1
            end
            
            max_count = -1
            maj_c = '-'
            for (c, count) in counts
                if count > max_count
                    max_count = count
                    maj_c = c
                end
            end
            consensus_chars[dj] = maj_c
        end
    end

    # Print the header line (histbars or consensus)
    if style == :allcolors
        msa_all_rows = _returnrows(msa)
        for dj in 1:length(displayed_cols_range)
            pos_counts = vec(get_base_count(msa_all_rows, dj))
            major_nuc = "ACGT"[argmax(pos_counts)]
            depth = only(msadepth(msa_all_rows, dj))
            bar_index = clamp(floor(Int, depth * 8) + 1, 1, 8)
            bar_char = histbars[bar_index]
            color = get(BASE_COLORS, major_nuc, :normal)
            printstyled(io, bar_char; color=color)
        end
    else
        for dj in 1:length(displayed_cols_range)
            c = consensus_chars[dj]
            if style == :polymorf
                color = get(BASE_COLORS, c, :normal)
                printstyled(io, c; color=color, reverse=true)
            else
                print(io, c)
            end
        end
    end

    if needs_width_ellipsis
        if style == :allcolors || style == :polymorf
            printstyled(io, "…"; color=:light_black, reverse=true)
        else
            print(io, "…")
        end
    end

    println(io)

    # Print the sequences
    for i in 1:n_display_seqs
        if has_desc
            print(io, padded_descs[i], " >")
        end
        for dj in 1:length(displayed_cols_range)
            c = getsequence(msa, i, dj)
            if style == :allcolors
                col = get(BASE_COLORS, c, :normal)
                printstyled(io, c; color=col, reverse=true)
            elseif style == :bw
                if c == consensus_chars[dj]
                    print(io, '.')
                else
                    print(io, c)
                end
            else # :polymorf
                if c == consensus_chars[dj]
                    print(io, '.')
                else
                    col = get(BASE_COLORS, c, :normal)
                    printstyled(io, c; color=col)
                end
            end
        end
        
        if needs_width_ellipsis
            if style == :allcolors || style == :polymorf
                printstyled(io, "…"; color=:light_black, reverse=true)
            else
                print(io, "…")
            end
        end
        println(io)
    end

    # Print the gap below sequence names if truncated
    gap_below_seqnames = fill(' ', desc_width)
    if n_display_seqs < n_sequences && length(gap_below_seqnames) ≥ 3
        gap_below_seqnames[1:3] .= '.'
    end
    has_desc && print(io, join(gap_below_seqnames))

    # Print the position number line
    num_line_chars = fill(' ', length(displayed_cols_range))
    
    start_num_str = string(first(abs_cols))
    start_num_len = length(start_num_str)
    if length(num_line_chars) >= start_num_len
        @inbounds @simd for i in 1:start_num_len
            if i <= length(num_line_chars)
                num_line_chars[i] = start_num_str[i]
            end
        end
    end
    
    end_num_str = string(last(abs_cols))
    end_num_len = length(end_num_str)
    if length(num_line_chars) >= end_num_len
        @inbounds @simd for i in 1:end_num_len
            idx = length(num_line_chars) - end_num_len + i
            if idx >= 1
                num_line_chars[idx] = end_num_str[i]
            end
        end
    end
    
    @inbounds for (rel_idx, abs_pos) in enumerate(abs_cols)
        if abs_pos % 10 == 0 && num_line_chars[rel_idx] == ' '
            num_line_chars[rel_idx] = '⋅'
        end
        if abs_pos % 100 == 0 && num_line_chars[rel_idx] == '⋅'
            num_line_chars[rel_idx] = '*'
        end
        if abs_pos % 1000 == 0 && num_line_chars[rel_idx] == '*'
            num_line_chars[rel_idx] = '#'
        end
    end
    
    needs_width_ellipsis && print(io, " ")
    print(io, String(num_line_chars))
end
_truncate_seq(seq::AbstractString, max_width::Int=20) = length(seq) > max_width ?
    seq[1:max(0, max_width-3)] * "..." : seq

function _show_header(io::IO, oligo::AbstractOligo)
    println(io, typeof(oligo))
    println(io, "  Sequence: ", String(oligo))
    println(io, "  Length: ", length(oligo))
end

function _show_common_fields(io::IO, oligo::AbstractOligo)
    print(io, "  Description: ")
    if isempty(description(oligo))
        print(io, "(none)")
    else
        print(io, "\"", description(oligo), "\"")
    end
end

function Base.show(io::IO, oligo::Oligo)
    seq_display = _truncate_seq(String(oligo))
    print(io, "Oligo(\"", seq_display, "\", len=", length(oligo))
    if !isempty(description(oligo))
        print(io, ", desc=\"", description(oligo), "\"")
    end
    print(io, ")")
end

function Base.show(io::IO, ::MIME"text/plain", oligo::Oligo)
    _show_header(io, oligo)
    _show_common_fields(io, oligo)
end

function Base.show(io::IO, deg::DegenOligo)
    seq_display = _truncate_seq(String(deg))
    print(io, "DegenOligo(\"", seq_display, "\", len=", length(deg))
    print(io, ", n_deg=", n_deg_pos(deg))
    vars_str = n_unique_oligos(deg) > 10000 ? ">10k" : string(n_unique_oligos(deg))
    print(io, ", vars=", vars_str)
    if !isempty(description(deg))
        print(io, ", desc=\"", description(deg), "\"")
    end
    print(io, ")")
end

function Base.show(io::IO, ::MIME"text/plain", deg::DegenOligo)
    _show_header(io, deg)
    println(io, "  Degenerate positions: ", n_deg_pos(deg))
    println(io, "  Unique variants: ", n_unique_oligos(deg))
    _show_common_fields(io, deg)
end

function Base.show(io::IO, go::GappedOligo)
    seq_display = _truncate_seq(String(go))
    print(io, "GappedOligo(\"", seq_display, "\", len=", length(go))
    print(io, ", gaps=", length(go.gaps))
    if !isempty(description(go))
        print(io, ", desc=\"", description(go), "\"")
    end
    print(io, ")")
end

function Base.show(io::IO, ::MIME"text/plain", go::GappedOligo)
    println(io, typeof(go))
    println(io, "  Gapped sequence: ", String(go))
    println(io, "  Length (with gaps): ", length(go))
    println(io, "  Gaps: ", length(go.gaps))
    _show_common_fields(io, go)
end

function Base.show(io::IO, ov::OligoView)
    seq_display = _truncate_seq(String(ov))
    print(io, "OligoView(\"", seq_display, "\", len=", length(ov))
    print(io, ", range=", ov.range)
    if !isempty(description(ov))
        print(io, ", desc=\"", description(ov), "\"")
    end
    print(io, ")")
end

function Base.show(io::IO, ::MIME"text/plain", ov::OligoView)
    println(io, typeof(ov))
    println(io, "  Viewed sequence: ", String(ov))
    println(io, "  Length: ", length(ov))
    println(io, "  Range: ", ov.range)
    println(io, "  Parent description: ", description(ov))
end
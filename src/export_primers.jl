# fatou-ignore-file undefined-name

"""
    export_evrogen(io::IO, primers; scale=0.04)
    export_evrogen(io::IO, pairs; scale=0.04)
    export_evrogen(filename::AbstractString, primers; scale=0.04)
    export_evrogen(filename::AbstractString, pairs; scale=0.04)

Export primers to a text stream or file formatted for the Evrogen DNA synthesis order form (Form I).

The format used is: `Name; Sequence; Scale` (e.g., `Primer_F_18; AGACYGACCGHGAAYTMGACCT; 0.04`).
IUPAC ambiguity codes are preserved, as required by Evrogen.

# Arguments
- `io::IO`: An output stream (e.g., `stdout` or a buffer).
- `filename::AbstractString`: Path to the output text file.
- `primers`: A single primer or a vector of primers.
- `pairs`: A single primer pair or a vector of primer pairs (e.g., output from `best_pairs`).
- `scale`: Synthesis scale (e.g., `0.04`, `0.2`, `1.0`). Defaults to `0.04`.

# Returns
- For `IO` methods: returns `nothing`;
- For file methods: returns the `filename`.
"""
function export_evrogen(
    io::IO,
    primers::AbstractVector{<:AbstractPrimer};
    scale::Union{Real, AbstractString} = 0.04,
)
    return _write_evrogen!(io, primers; scale = scale)
end

function export_evrogen(
    io::IO,
    primer::AbstractPrimer;
    scale::Union{Real, AbstractString} = 0.04,
)
    return _write_evrogen!(io, [primer]; scale = scale)
end

function export_evrogen(
    io::IO,
    pairs::AbstractVector{<:Pair};
    scale::Union{Real, AbstractString} = 0.04,
)
    return _write_evrogen!(io, pairs; scale = scale)
end

function export_evrogen(
    io::IO,
    pair::Pair{<:AbstractPrimer, <:AbstractPrimer};
    scale::Union{Real, AbstractString} = 0.04,
)
    return _write_evrogen!(io, [pair]; scale = scale)
end

function export_evrogen(
    filename::AbstractString,
    primers::AbstractVector{<:AbstractPrimer};
    scale::Union{Real, AbstractString} = 0.04,
)
    open(filename, "w") do io
        _write_evrogen!(io, primers; scale = scale)
    end
    return filename
end

function export_evrogen(
    filename::AbstractString,
    primer::AbstractPrimer;
    scale::Union{Real, AbstractString} = 0.04,
)
    open(filename, "w") do io
        _write_evrogen!(io, [primer]; scale = scale)
    end
    return filename
end

function export_evrogen(
    filename::AbstractString,
    pairs::AbstractVector{<:Pair};
    scale::Union{Real, AbstractString} = 0.04,
)
    open(filename, "w") do io
        _write_evrogen!(io, pairs; scale = scale)
    end
    return filename
end

function export_evrogen(
    filename::AbstractString,
    pair::Pair{<:AbstractPrimer, <:AbstractPrimer};
    scale::Union{Real, AbstractString} = 0.04,
)
    open(filename, "w") do io
        _write_evrogen!(io, [pair]; scale = scale)
    end
    return filename
end

function _write_evrogen!(
    io::IO,
    primers::AbstractVector{<:AbstractPrimer};
    scale::Union{Real, AbstractString} = 0.04,
)
    for (i, p) in enumerate(primers)
        name = _evrogen_name(p, i)
        seq = String(p)
        println(io, "$name; $seq; $scale")
    end
    return nothing
end

function _write_evrogen!(
    io::IO,
    pairs::AbstractVector{<:Pair};
    scale::Union{Real, AbstractString} = 0.04,
)
    idx = 1
    for pair in pairs
        name_f = _evrogen_name(pair.first, idx)
        seq_f = String(pair.first)
        println(io, "$name_f; $seq_f; $scale")

        name_r = _evrogen_name(pair.second, idx)
        seq_r = String(pair.second)
        println(io, "$name_r; $seq_r; $scale")

        idx += 1
    end
    return nothing
end

function _evrogen_name(p::AbstractPrimer, idx::Int)
    base_desc = replace(description(p), r"[;\t\n\r]" => " ")
    dir = p.is_forward ? "F" : "R"

    pos = p.is_forward ? first(p.pos) : last(p.pos)

    if isempty(base_desc)
        return "DePPA_$(dir)_$(pos)_$(idx)"
    else
        return "$(base_desc)_$(dir)_$(pos)"
    end
end

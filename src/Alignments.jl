module Alignments
include("utils.jl")

using ..Oligos
using ProgressMeter
using FastaIO
using Random

export AbstractMSA, MSA, MSAView
export nseqs, width, height, getsequence, get_base_count
export msadepth, msadet, root, bval
export consensus_major, consensus_degen, dry_msa, nucleotide_diversity
export setMSAShowStyle!

function _bootstrap_base_counts(
    seqs::Vector{<:AbstractString}, 
    bootstrap::Int;
    progress_label::String, 
    barlen::Int
)
    n = length(seqs)
    L = length(first(seqs))
    all(length(s) == L for s in seqs) || throw(ArgumentError("All sequences must have the same length"))
    base_count = zeros(4, L)
    if bootstrap > 0
        @showprogress desc=progress_label barlen=barlen for b in 1:bootstrap
            boot_rows = rand(1:n, n)
            boot_seqs = seqs[boot_rows]
            
            Threads.@threads for j in 1:L
                pos_counts = zeros(4)
                for s in boot_seqs
                    c = uppercase(s[j])
                    probs = get(IUPAC_PROBS, c, (.0,.0,.0,.0))
                    pos_counts .+= probs
                end
                current_freq = pos_counts / n
                base_count[:, j] += (current_freq - base_count[:, j]) / b
            end
        end
    else
        Threads.@threads for j in 1:L
            counts = zeros(4)
            for s in seqs
                c = uppercase(s[j])
                probs = get(IUPAC_PROBS, c, (.0,.0,.0,.0))
                counts .+= probs
            end
            base_count[:, j] = counts / n
        end
    end
    return base_count
end

"""
    AbstractMSA

Abstract supertype for Multiple Sequence Alignments.
"""
abstract type AbstractMSA end

"""
    MSA <: AbstractMSA

A concrete Multiple Sequence Alignment type. Stores sequences and precomputed base frequencies.
"""
struct MSA <: AbstractMSA
    seqs::Vector{<:AbstractGapped}
    base_count::Matrix{Float64}
    bootstrap::Int

    """
        MSA(seqs::Vector{<:AbstractString}; bootstrap::Int=0, seed=nothing)

    Construct an MSA from a vector of equal-length strings.

    If `bootstrap > 0`, compute base frequencies using bootstrap resampling.
    """
    function MSA(seqs::Vector{<:AbstractString}; bootstrap::Int=0, seed=nothing)
        bootstrap >= 0 || throw(ArgumentError("bootstrap must be non-negative"))
        isnothing(seed) || Random.seed!(seed)

        isempty(seqs) && return new(GappedOligo[], zeros(4, 0), bootstrap)

        base_count = _bootstrap_base_counts(seqs, bootstrap;
            progress_label="Bootstrap, $bootstrap it.",
            barlen=19
        )
        gapped_seqs = GappedOligo.(seqs)
        return new(gapped_seqs, base_count, bootstrap)
    end
end

"""
    MSAView <: AbstractMSA

A lightweight view into an [`MSA`](@ref), representing a submatrix of rows and columns.
"""
struct MSAView <: AbstractMSA
    parent::AbstractMSA
    rows::UnitRange{Int}
    cols::UnitRange{Int}
end


_returnrows(m::AbstractMSA) = m
_returnrows(m::MSAView) = MSAView(m.parent, 1:height(m.parent), m.cols)

"""
    root(msa::AbstractMSA) -> MSA

Return the underlying root `MSA` object, resolving any `MSAView` layers.

See also [`bval`](@ref), [`MSAView`](@ref).
"""
root(msa::MSA) = msa
root(msav::MSAView) = root(msav.parent)

"""
    bval(msa::AbstractMSA) -> Int

Return the number of bootstrap iterations used to compute base frequencies.

See also [`root`](@ref), [`MSA`](@ref).
"""
bval(msa::MSA) = msa.bootstrap
bval(msav::MSAView) = root(msav).bootstrap

function _align!(args...; kwargs...)
    # This is overloaded in ext/MAFFTExt.jl to load MAFFT_jll artifact dynamically
    error("Alignment requires MAFFT_jll artifact to be loaded.\n" *
          "In order to align your FASTAs, please `]add MAFFT_jll` to your project\n" *
          "and load it with `using MAFFT_jll` before calling MSA with `mafft=true`.")
end

"""
    MSA(predicate::Function, fasta::AbstractString; mafft::Bool=false, bootstrap::Int=0, seed=nothing)

Construct an MSA from a FASTA file.

# Arguments
- `predicate::Function`: A function that takes the sequence description and returns `true` to include the sequence.
- `fasta::AbstractString`: Path to the FASTA file.
- `mafft::Bool=false`: If `true`, align the sequences using MAFFT (requires `MAFFT_jll` package to be loaded).
- `bootstrap::Int=0`: Number of bootstrap iterations for base frequencies.
- `seed=nothing`: Random seed for reproducibility.
"""
function MSA(predicate::Function, fasta::AbstractString; mafft::Bool=false, bootstrap::Int=0, seed=nothing)
    fasta_content = Tuple{String, String}[]
    FastaReader(fasta) do fr
        counter = 0
        for (desc, seq) in fr
            predicate(desc) ? (counter += 1) : continue
            desc = isempty(desc) ? "seq$counter" : desc
            push!(fasta_content, (desc, seq))
        end
    end

    allowed_chars = mafft ? (NON_DEGEN_BASES..., 'N') : (NON_DEGEN_BASES..., '-', 'N')
    for (desc, seq) in fasta_content
        upper_seq = uppercase(seq)
        if !all(in(allowed_chars), upper_seq)
            invalid_chars = setdiff(unique(upper_seq), collect(allowed_chars))
            throw(ArgumentError(
                "Sequence '$desc' contains invalid characters: $(join(invalid_chars, ", ")).\\n" *
                "Only $(join(collect(allowed_chars), ", ")) allowed when `mafft=$mafft`"
            ))
        end
    end

    mafft && _align!(fasta_content)

    gapped_oligos = [GappedOligo(seq, desc) for (desc, seq) in fasta_content]
    return MSA(gapped_oligos; bootstrap=bootstrap, seed=seed)
end

"""
    MSA(fasta::AbstractString; kwargs...)

Construct an `MSA` from a FASTA file, including all sequences.

This is a convenience method that delegates to the predicate-based constructor with a filter that always returns `true`.

# Arguments
- `fasta::AbstractString`: Path to the FASTA file.
- `kwargs...`: Keyword arguments passed to the underlying `MSA(predicate, fasta; kwargs...)` constructor (e.g., `mafft::Bool`, `bootstrap::Int`, `seed`).

# Returns
- `MSA`: A new `MSA` object containing all sequences from the file.

"""
MSA(fasta::AbstractString; kwargs...) = MSA(x->true, fasta; kwargs...)

"""
    MSA(msav::MSAView; bootstrap::Int=0, seed=nothing)

Construct a new concrete `MSA` by materializing the sliced rows and columns from an `MSAView` into a standalone alignment.

# Arguments
- `msav::MSAView`: The MSA view to materialize.
- `bootstrap::Int=0`: Number of bootstrap iterations for computing base frequencies.
- `seed=nothing`: Random seed for reproducibility during bootstrap resampling.

# Returns
- `MSA`: A new, independent `MSA` object containing the sliced sequences.

See also [`MSAView`](@ref).
"""
function MSA(msav::MSAView; bootstrap::Int=0, seed=nothing)
    seqs = [getsequence(msav, i) for i in 1:nseqs(msav)]
    return MSA(seqs; bootstrap=bootstrap, seed=seed)
end

"""
    getindex(msa::AbstractMSA, rows, cols)

Get a submatrix or element from an MSA using multi-dimensional indexing.

Supports:
- `msa[row, col]`: single element
- `msa[row_range, col_range]`: submatrix view
- `msa[row, :]`: entire row
- `msa[:, col_range]`: entire column range
"""
function Base.getindex(msa::AbstractMSA, rows::UnitRange{Int}, cols::UnitRange{Int})
    root_msa = root(msa)
    abs_rows = isa(msa, MSA) ? rows : (msa.rows.start + rows.start - 1):(msa.rows.start + rows.stop - 1)
    abs_cols = isa(msa, MSA) ? cols : (msa.cols.start + cols.start - 1):(msa.cols.start + cols.stop - 1)
    return MSAView(root_msa, abs_rows, abs_cols)
end
Base.getindex(msa::AbstractMSA, row::Int, col::Int) = getsequence(msa, row, col)
Base.getindex(msa::AbstractMSA, row::Int, cols::UnitRange{Int}) = msa[row:row, cols]
Base.getindex(msa::AbstractMSA, row::Int, ::Colon) = getsequence(msa, row)
Base.getindex(msa::AbstractMSA, rows::UnitRange{Int}, col::Int) = msa[rows, col:col]
Base.getindex(msa::AbstractMSA, ::Colon, col::Int) = msa[:, col:col]
Base.getindex(msa::AbstractMSA, ::Colon, cols::UnitRange{Int}) = msa[1:nseqs(msa), cols]
Base.getindex(msa::AbstractMSA, rows::UnitRange{Int}, ::Colon) = msa[rows, 1:length(msa)]
Base.getindex(msa::AbstractMSA, ::Colon, ::Colon) = msa[1:nseqs(msa), 1:length(msa)]

function Base.checkbounds(msa::AbstractMSA, rows::Colon, cols::UnitRange{<:Integer})
    if ! (1 <= cols.start <= cols.stop <= length(msa))
        throw(BoundsError(msa, (:, cols)))
    end
end
function Base.checkbounds(msa::AbstractMSA, rows::UnitRange{<:Integer}, cols::Colon)
    if ! (1 <= rows.start <= rows.stop <= nseqs(msa))
        throw(BoundsError(msa, (rows, :)))
    end
end
function Base.checkbounds(msa::AbstractMSA, rows::UnitRange{<:Integer}, cols::UnitRange{<:Integer})
    if ! (1 <= rows.start <= rows.stop <= nseqs(msa))
        throw(BoundsError(msa, (rows, cols)))
    end
    if ! (1 <= cols.start <= cols.stop <= length(msa))
        throw(BoundsError(msa, (rows, cols)))
    end
end
function Base.checkbounds(::AbstractMSA, ::Colon, ::Colon) end

"""
    nseqs(msa::AbstractMSA) -> Int

Return the number of sequences (rows) in the alignment.

See also [`height`](@ref), [`width`](@ref), [`length`](@ref).
"""
nseqs(msa::MSA) = length(msa.seqs)
nseqs(v::MSAView) = length(v.rows)

Base.length(msa::MSA) = size(msa.base_count, 2)
Base.length(v::MSAView) = length(v.cols)

"""
    width(msa::AbstractMSA) -> Int

Return the number of columns (alignment length) in the MSA.

See also [`length`](@ref), [`nseqs`](@ref), [`height`](@ref).
"""
width(msa::AbstractMSA) = length(msa)

"""
    height(msa::AbstractMSA) -> Int

Return the number of sequences (rows) in the MSA. Alias for [`nseqs`](@ref).

See also [`nseqs`](@ref), [`width`](@ref).
"""
height(msa::AbstractMSA) = nseqs(msa)

Base.ndims(::AbstractMSA) = 2

Base.size(msa::AbstractMSA) = (nseqs(msa), length(msa))
function Base.size(msa::AbstractMSA, dim::Int)
    if dim == 1
        nseqs(msa)
    elseif dim == 2
        length(msa)
    else
        throw(ArgumentError("AbstractMSA only has dimensions 1 and 2"))
    end
end

Base.axes(msa::AbstractMSA, dim::Int) = Base.OneTo(size(msa, dim))

Base.lastindex(msa::AbstractMSA, dim::Int) = size(msa, dim)
Base.lastindex(msa::AbstractMSA) = lastindex(msa, ndims(msa))

"""
    getsequence(msa::AbstractMSA, row::Int)
    getsequence(msa::AbstractMSA, row::Int, col::Int)

Get a sequence or an individual position from an MSA.

# Arguments
- `msa::AbstractMSA`: The MSA.
- `row::Int`: Sequence index (1-based).
- `col::Int`: Position index (1-based, optional).

# Returns
- For `getsequence(msa, row)`: The full sequence ([`GappedOligo`](@ref)).
- For `getsequence(msa, row, col)`: Single character at position.
"""
getsequence(msa::MSA, row::Int) = msa.seqs[row]

function getsequence(v::MSAView, row::Int)
    abs_row = v.rows.start + row - 1
    parent_seq = getsequence(root(v), abs_row)
    return parent_seq[v.cols]
end
getsequence(msa::AbstractMSA, row::Int, col::Int) = getsequence(msa, row)[col]


_is_full_height(msa::MSA) = true
_is_full_height(msav::MSAView) = msav.rows == 1:nseqs(root(msav))

"""
    get_base_count(msa::AbstractMSA, pos::Int)
    get_base_count(msa::AbstractMSA, interval::UnitRange{Int})
    get_base_count(msa::AbstractMSA)

Get base frequency counts from an MSA.

# Arguments
- `msa::AbstractMSA`: The MSA.
- `pos::Int`: Single position (1-based).
- `interval::UnitRange{Int}`: Range of positions.
- If no position or interval is provided, returns counts for all positions.

# Returns
- A vector of 4 floats (A, C, G, T probabilities) for a single position.
- A matrix view for multiple positions.
"""
get_base_count(msa::MSA, pos::Int) = @view msa.base_count[:, pos]
get_base_count(msa::MSA, interval::UnitRange{Int}) = @view msa.base_count[:, interval]
get_base_count(msa::MSA) = msa.base_count
function get_base_count(msav::MSAView, pos::Int)
    if !_is_full_height(msav)
        throw(ErrorException("get_base_count not supported for views that slice rows (height), as it requires recomputation"))
    end
    abs_pos = msav.cols.start + pos - 1
    return @view root(msav).base_count[:, abs_pos]
end
function get_base_count(msav::MSAView, interval::UnitRange{Int})
    if !_is_full_height(msav)
        throw(ErrorException("get_base_count not supported for views that slice rows (height), as it requires recomputation"))
    end
    abs_start = msav.cols.start + interval.start - 1
    abs_interval = abs_start:(abs_start + length(interval) - 1)
    return @view root(msav).base_count[:, abs_interval]
end
get_base_count(msav::MSAView) = get_base_count(msav, 1:length(msav))

"""
    msadepth(msa::AbstractMSA, pos::Int)
    msadepth(msa::AbstractMSA, interval::UnitRange{Int})
    msadepth(msa::AbstractMSA)

Calculate sequence depth (coverage) at positions. Depth is the sum of base probabilities, capped at 1.0.

# Arguments
- `msa::AbstractMSA`: The MSA.
- `pos::Int`: Single position.
- `interval::UnitRange{Int}`: Range of positions.
- If no position or interval is provided, calculate depth for all positions.

# Returns
- `Float64` for a single position.
- `Vector{Float64}` for multiple positions.

See also [`msadet`](@ref), [`get_base_count`](@ref).
"""
function msadepth(msa::AbstractMSA, pos::Int)::Float64
    min(1.0, sum(get_base_count(msa, pos)))
end
function msadepth(msa::AbstractMSA, interval::UnitRange{Int})::Vector{Float64}
    [msadepth(msa, pos) for pos in interval]
end
function msadepth(msa::AbstractMSA)::Vector{Float64}
    return msadepth(msa, 1:length(msa))
end

"""
    msadet(msa::AbstractMSA, pos::Int)
    msadet(msa::AbstractMSA, interval::UnitRange{Int})
    msadet(msa::AbstractMSA)

Calculate sequence determinacy (entropy inverse) at positions. Determinacy is the maximum base frequency normalized by total coverage.

# Arguments
- `msa::AbstractMSA`: The MSA.
- `pos::Int`: Single position.
- `interval::UnitRange{Int}`: Range of positions.
- If no position or interval is provided, calculate determinacy for all positions.

# Returns
- `Float64` for a single position (0.0 to 1.0).
- `Vector{Float64}` for multiple positions.

See also [`msadepth`](@ref), [`get_base_count`](@ref).
"""
function msadet(msa::AbstractMSA, pos::Int)::Float64
    v = get_base_count(msa, pos)
    s = min(1.0, sum(v))
    s == 0.0 ? 0.0 : maximum(v) / s
end
function msadet(msa::AbstractMSA, interval::UnitRange{Int})::Vector{Float64}
    return [msadet(msa, pos) for pos in interval]
end
function msadet(msa::AbstractMSA)::Vector{Float64}
    return msadet(msa, 1:length(msa))
end

"""
    consensus_major(msa::AbstractMSA, pos::Int)
    consensus_major(msa::AbstractMSA, interval::UnitRange{Int}=1:width(msa))

Generate a majority-rule consensus sequence using simple majority rule, ignoring gap characters.

# Arguments
- `msa::AbstractMSA`: The MSA.
- `pos::Int`: Single position.
- `interval::UnitRange{Int}=1:width(msa)`: Range of positions.

# Returns
- `Char` for a single position (most common base).
- `GappedOligo` for multiple positions.

See also [`consensus_degen`](@ref), [`get_base_count`](@ref).
"""
function consensus_major(msa::AbstractMSA, pos::Int)
    p = get_base_count(msa, pos)
    if sum(p) == 0
        return '-'
    end
    return NON_DEGEN_BASES[argmax(p)]
end
function consensus_major(msa::AbstractMSA, interval::UnitRange{Int}=1:width(msa))
    seq = join(consensus_major(msa, j) for j in interval)
    desc = "Major consensus for $(nseqs(msa)) seq MSA"
    return GappedOligo(seq, desc)
end

"""
    consensus_degen(msa::AbstractMSA, pos::Int; slack::Real=0.0)
    consensus_degen(msa::AbstractMSA, interval::UnitRange{Int}=1:width(msa); slack::Real=0.0)

Generate a degenerate consensus sequence allowing ambiguity. Bases with frequency > `slack` are included in the degeneracy.

# Arguments
- `msa::AbstractMSA`: The MSA.
- `pos::Int`: Single position.
- `interval::UnitRange{Int}=1:width(msa)`: Range of positions.
- `slack::Real=0.0`: Minimum frequency threshold for inclusion.

# Returns
- `Char` for a single position (IUPAC ambiguity code).
- `GappedOligo` for multiple positions.

See also [`consensus_major`](@ref), [`get_base_count`](@ref).
"""
function consensus_degen(msa::AbstractMSA, pos::Int; slack::Real=0.0)::Char
    0 ≤ slack < 1 || throw(ArgumentError("slack must be in [0,1)"))
    p = get_base_count(msa, pos)
    if sum(p) == 0
        return '-'
    end
    active = findall(>(slack), p)
    isempty(active) && return '-'
    bs = NON_DEGEN_BASES[active]
    return IUPAC_V2B[bs]
end
function consensus_degen(msa::AbstractMSA, interval::UnitRange{Int}=1:width(msa); slack::Real=0.0)::GappedOligo
    seq = join(consensus_degen(msa, j; slack=slack) for j in interval)
    desc = "Degenerate consensus for $(nseqs(msa)) seq MSA"
    return GappedOligo(seq, desc)
end

"""
    dry_msa(msa::AbstractMSA; gap_content::Real=1.0)

Remove columns and rows with excessive gap content. Columns with no non-gap characters are always removed. Rows with gap proportion > `gap_content` are removed.

# Arguments
- `msa::AbstractMSA`: The MSA.
- `gap_content::Real=1.0`: Maximum allowed gap proportion (default: 1.0, keep all).

# Returns
- A new [`MSA`](@ref) with filtered sequences and columns.
"""
function dry_msa(msa::AbstractMSA; gap_content::Real=1.0)
    0 ≤ gap_content ≤ 1 || throw(ArgumentError("gap_content must be in [0,1]"))
    non_gap_cols = [j for j in 1:length(msa) if any(>(0.0), get_base_count(msa, j))]
    if isempty(non_gap_cols)
        return msa[:, 1:0]
    end
    num_cols = length(non_gap_cols)
    kept_rows = Int[]
    for i in 1:nseqs(msa)
        gap_count = sum((1 for j in non_gap_cols if getsequence(msa, i, j) == '-'), init=0)
        prop = num_cols > 0 ? gap_count / num_cols : 0.0
        if prop < gap_content
            push!(kept_rows, i)
        end
    end
    if isempty(kept_rows)
        return MSA(GappedOligo[]; bootstrap=bval(msa))
    end
    new_seqs = Vector{GappedOligo}(undef, length(kept_rows))
    for (k, row) in enumerate(kept_rows)
        sub_str = join(getsequence(msa, row, j) for j in non_gap_cols)
        new_seqs[k] = GappedOligo(sub_str, description(getsequence(msa, row)))
    end
    return MSA(new_seqs; bootstrap=bval(msa))
end

"""
    nucleotide_diversity(msa::AbstractMSA; ignore_gaps::Bool=true, max_pairs::Int=10000)

Calculate average pairwise nucleotide diversity. Uses probabilistic distance for degenerate bases. For large MSAs (>200 sequences), samples random pairs.

# Arguments
- `msa::AbstractMSA`: The MSA.
- `ignore_gaps::Bool=true`: Whether to skip gap-gap comparisons.
- `max_pairs::Int=10000`: Maximum pairs to sample for large MSAs.

# Returns
- `Float64`: Average pairwise distance.
"""
function nucleotide_diversity(msa::AbstractMSA; ignore_gaps::Bool=true, max_pairs::Int=10000)::Float64
    L = length(msa)
    L == 0 && return 0.0
    n = nseqs(msa)
    n < 2 && return 0.0

    total_possible = n * (n - 1) ÷ 2
    compute_all = n <= 200 || max_pairs >= total_possible
    npairs = compute_all ? total_possible : min(max_pairs, total_possible)

    pairs = Vector{Tuple{Int,Int}}(undef, npairs)
    @inbounds if compute_all
        idx = 1
        for i in 1:n-1, j in i+1:n
            pairs[idx] = (i, j)
            idx += 1
        end
    else
        for idx in 1:npairs
            i, j = rand(1:n), rand(1:n)
            while i == j
                i, j = rand(1:n), rand(1:n)
            end
            pairs[idx] = minmax(i, j)
        end
    end

    total_diffs = Threads.Atomic{Float64}(0.0)
    Threads.@threads for (i, j) in pairs
        d = _pairwise_distance(msa, i, j; ignore_gaps)
        Threads.atomic_add!(total_diffs, d)
    end

    return total_diffs[] / npairs
end

"""
    _pairwise_distance(msa::AbstractMSA, i::Int, j::Int; ignore_gaps::Bool=true)

Calculate pairwise distance between two sequences using probabilistic matching for degenerate bases.

# Arguments
- `msa::AbstractMSA`: The MSA.
- `i::Int`, `j::Int`: Sequence indices.
- `ignore_gaps::Bool=true`: Whether to skip gap positions.

# Returns
- `Float64`: Normalized distance (0.0 to 1.0).
"""
function _pairwise_distance(msa::AbstractMSA, i::Int, j::Int; ignore_gaps::Bool=true)::Float64
    seq_i = getsequence(msa, i)
    seq_j = getsequence(msa, j)
    total_sites = 0
    diff_sum = 0.0
    for k in 1:length(msa)
        c_i = seq_i[k]
        c_j = seq_j[k]
        if c_i == '-' || c_j == '-' 
            ignore_gaps && continue
            # Treat gap as mismatch
            total_sites += 1
            diff_sum += 1.0
            continue
        end
        total_sites += 1
        if c_i == c_j
            continue  # Exact match
        end
        # Probabilistic mismatch for IUPAC/degenerate
        probs_i = get(IUPAC_PROBS, c_i, (.0,.0,.0,.0))
        probs_j = get(IUPAC_PROBS, c_j, (.0,.0,.0,.0))
        match_prob = sum(probs_i .* probs_j)
        diff_sum += 1.0 - match_prob
    end
    total_sites == 0 ? 0.0 : diff_sum / total_sites
end

include("show_msa.jl")
end # module

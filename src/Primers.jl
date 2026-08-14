module Primers

export AbstractPrimer
export Primer
export construct_primers, best_pairs, export_evrogen
export miniblast, MiniBlastHit
export add_illumina_adapters, reannotated
export setAdapters!, getAdapters

using ..Oligos
using ..Alignments

using ProgressMeter
using Statistics

const GLOBAL_ADAPTERS = Ref{Union{Nothing,Pair{Oligo,Oligo}}}(nothing)

"""
    setAdapters!() -> Nothing
    setAdapters!(adapters::Pair{<:Oligo, <:Oligo}) -> Nothing
    setAdapters!(adapters::Pair{<:AbstractString, <:AbstractString}) -> Nothing

Set the global adapter sequences to be automatically appended to the 5' ends of primers 
during `construct_primers`. 

# Arguments
- `()`: Resets the global adapters to `nothing` (no adapters will be added).
- `adapters::Pair{<:Oligo, <:Oligo}`: A custom pair of valid `Oligo` sequences.
- `adapters::Pair{<:AbstractString, <:AbstractString}`: A custom pair of strings (will be converted to `Oligo`).

# Details
When global adapters are set, `construct_primers` will automatically concatenate them to 
candidate primers, recalculate ΔG at the primer's mean Tm, and discard candidates where 
the adapter worsens ΔG by more than `max_dg_drop`.
"""
function setAdapters!()::Nothing
    GLOBAL_ADAPTERS[] = nothing
    return nothing
end

function setAdapters!(adapters::Pair{<:Oligo,<:Oligo})::Nothing
    GLOBAL_ADAPTERS[] = adapters
    return nothing
end

function setAdapters!(adapters::Pair{<:AbstractString,<:AbstractString})::Nothing
    setAdapters!(Oligo(adapters.first, "Custom Fwd Adapter") => Oligo(
        adapters.second,
        "Custom Rev Adapter",
    ))
end

getAdapters()::Union{Nothing,Pair{Oligo,Oligo}} = GLOBAL_ADAPTERS[]

"""
    AbstractPrimer{T<:Union{Oligo,DegenOligo}}

Represent the abstract supertype for PCR primers.

See also [`Primer`](@ref).
"""
abstract type AbstractPrimer{T<:Union{Oligo,DegenOligo}} end

"""
    Primer{T} <: AbstractPrimer{T}

Represent a concrete PCR primer. Store the consensus sequence, its position in the MSA,
and thermodynamic properties (Tm, ΔG, GC content).

See also [`AbstractPrimer`](@ref), [`construct_primers`](@ref).
"""
struct Primer{T} <: AbstractPrimer{T}
    msa::AbstractMSA
    pos::UnitRange{Int}
    is_forward::Bool
    consensus::T
    tail_length::Int
    tm::@NamedTuple{
        mean::Float64,
        conf::Tuple{Float64,Float64},
        min::Float64,
        max::Float64,
    }
    dg::Float64
    gc::Float64
    slack::Float64
    adapter::Union{Nothing,Oligo}
end

"""
    Primer{T}(msa::AbstractMSA, pos::UnitRange{Int}, is_forward::Bool, consensus, tail_length::Int, tm, dg::Float64, gc::Float64, slack::Float64) where T

Constructs a `Primer` with type parameter `{T}` without specifying an adapter.
Automatically converts the `consensus` to type `T` and sets the `adapter` field to `nothing`.
"""
function Primer{T}(
    msa::AbstractMSA,
    pos::UnitRange{Int},
    is_forward::Bool,
    consensus,
    tail_length::Int,
    tm,
    dg::Float64,
    gc::Float64,
    slack::Float64,
)::Primer{T} where {T<:Union{Oligo,DegenOligo}}
    return Primer{T}(
        msa,
        pos,
        is_forward,
        convert(T, consensus),
        tail_length,
        tm,
        dg,
        gc,
        slack,
        nothing,
    )
end

"""
    Primer(msa::AbstractMSA, pos::UnitRange{Int}, is_forward::Bool, consensus::T, tail_length::Int, tm, dg::Float64, gc::Float64, slack::Float64) where T

Constructs a `Primer` without specifying the type parameter `{T}` or an adapter.
Infers the type `T` directly from the `consensus` sequence and sets the `adapter` field to `nothing`.
"""
function Primer(
    msa::AbstractMSA,
    pos::UnitRange{Int},
    is_forward::Bool,
    consensus::T,
    tail_length::Int,
    tm,
    dg::Float64,
    gc::Float64,
    slack::Float64,
)::Primer{T} where {T<:Union{Oligo,DegenOligo}}
    return Primer{T}(
        msa,
        pos,
        is_forward,
        consensus,
        tail_length,
        tm,
        dg,
        gc,
        slack,
        nothing,
    )
end

"""
    Primer(msa::AbstractMSA, interval::UnitRange{Int}; kwargs...)

Construct a `Primer` object for a given interval in the MSA, calculating its thermodynamic properties.

If global adapters are set via [`setAdapters!`](@ref), they are automatically appended to the 5' end. 
The ΔG of the full sequence (adapter + primer) is calculated at the primer's mean Tm. 
If the adapter worsens ΔG by more than `max_dg_drop`, a warning is issued.

# Arguments
- `msa::AbstractMSA`: The multiple sequence alignment.
- `interval::UnitRange{Int}`: The position range of the primer in the MSA.
- `is_forward::Bool=true`: Design a forward (`true`) or reverse (`false`) primer.
- `tail_length::Int=3`: Length of the 3' tail region.
- `max_samples::Int=1000`: Number of samples for Monte Carlo estimation of Tm and ΔG.
- `tm_conf_int=0.8`: Confidence interval for Tm.
- `tm_conds=:pcr`: Thermodynamic conditions for Tm calculation.
- `dg_temp=37.0`: Temperature for ΔG calculation (used for the primer without adapter).
- `slack=0.0`: Minimum frequency threshold for including a base in the degenerate consensus.
- `max_dg_drop::Real=1.0`: Threshold for warning if the adapter worsens ΔG significantly.
- `descr`: Description string for the primer.

See also [`construct_primers`](@ref), [`consensus_degen`](@ref), [`setAdapters!`](@ref).
"""
function Primer(
    msa::AbstractMSA,
    interval::UnitRange{Int};
    is_forward::Bool=true,
    tail_length::Int=3,
    max_samples::Int=1000,
    tm_conf_int=0.8,
    tm_conds=:pcr,
    dg_temp=37.0,
    slack=0.0,
    max_dg_drop::Real=1.0,
    adapter_pair = GLOBAL_ADAPTERS[],
    descr="Primer for $(nseqs(msa)) seq MSA at positions $interval",
)::Primer{DegenOligo}
    _cons = consensus_degen(msa, interval; slack=slack)
    gapped_cons = is_forward ? _cons : _ext_revcomp(_cons)

    if hasgaps(gapped_cons)
        throw(ArgumentError(
            "Cannot construct primer for interval $interval: " *
            "consensus contains gaps (excessive gaps in MSA region).",
        ))
    end

    underlying_oligo = DegenOligo(String(gapped_cons), string(descr))

    Tm = _ext_tm(
        underlying_oligo;
        max_samples=max_samples,
        conf_int=tm_conf_int,
        conditions=tm_conds,
    )
    delta_G = _ext_dg(underlying_oligo; max_samples=max_samples, temp=dg_temp)
    GC = _ext_gc_content(underlying_oligo)

    final_dg = delta_G
    adapter_to_add = nothing

    if !isnothing(adapter_pair)
        adapter_to_add = is_forward ? adapter_pair.first : adapter_pair.second
        full_oligo = adapter_to_add * underlying_oligo

        final_dg = _ext_dg(full_oligo; max_samples=max_samples, temp=Tm.mean)

        if (final_dg - delta_G) < -max_dg_drop
            @warn "Primer at $interval with adapter has a significant ΔG drop: $(round(final_dg - delta_G, digits = 2)) kcal/mol"
        end
    end

    return Primer{DegenOligo}(
        msa,
        interval,
        is_forward,
        underlying_oligo,
        tail_length,
        Tm,
        final_dg,
        GC,
        slack,
        adapter_to_add,
    )
end

"""
    reannotated(primer::AbstractPrimer, annotation::AbstractString) -> Primer
    reannotated(pair::Pair{<:AbstractPrimer, <:AbstractPrimer}, annotation::AbstractString) -> Pair{Primer}

Create a new primer (or primer pair) with the updated description (annotation).
Since Julia structs are immutable, a new object is returned rather than mutating the existing one in-place.

# Arguments
- `primer` / `pair`: A single primer or a primer pair.
- `annotation::AbstractString`: The new description string.

# Returns
- A new `Primer` or `Pair{Primer, Primer}` with the updated description.
"""
function reannotated(primer::Primer{T}, annotation::AbstractString)::Primer{T} where {T<:Union{Oligo,DegenOligo}}
    new_consensus = T(String(primer.consensus), String(annotation))

    return Primer{T}(
        primer.msa,
        primer.pos,
        primer.is_forward,
        new_consensus,
        primer.tail_length,
        primer.tm,
        primer.dg,
        primer.gc,
        primer.slack,
        primer.adapter,
    )
end

function reannotated(
    pair::Pair{<:AbstractPrimer,<:AbstractPrimer},
    annotation::AbstractString,
)
    return reannotated(pair.first, annotation) => reannotated(pair.second, annotation)
end

# These are overloaded in ext/SeqFoldExt.jl to load SeqFold.jl library dynamically
_ext_revcomp(args...; kwargs...) = error(
    "Primer construction requires SeqFold library to be loaded.\n" *
    "In order to get this functionality, please `]add SeqFold` to your project\n" *
    "and load it with `using SeqFold` before constructing primers.",
)
_ext_tm(args...; kwargs...) = error(
    "Primer construction requires SeqFold library to be loaded.\n" *
    "In order to get this functionality, please `]add SeqFold` to your project\n" *
    "and load it with `using SeqFold` before constructing primers.",
)
_ext_dg(args...; kwargs...) = error(
    "Primer construction requires SeqFold library to be loaded.\n" *
    "In order to get this functionality, please `]add SeqFold` to your project\n" *
    "and load it with `using SeqFold` before constructing primers.",
)
_ext_gc_content(args...; kwargs...) = error(
    "Primer construction requires SeqFold library to be loaded.\n" *
    "In order to get this functionality, please `]add SeqFold` to your project\n" *
    "and load it with `using SeqFold` before constructing primers.",
)

Base.String(primer::AbstractPrimer)::String = isnothing(primer.adapter) ?
                                              String(primer.consensus) :
                                              String(primer.adapter) * String(primer.consensus)
Base.length(primer::AbstractPrimer)::Int = length(primer.consensus) +
                                           (isnothing(primer.adapter) ? 0 : length(primer.adapter))
Base.isempty(primer::AbstractPrimer)::Bool = isempty(primer.consensus)
Base.iterate(primer::AbstractPrimer, state...) = iterate(String(primer), state...)
Base.getindex(primer::AbstractPrimer, i::Int)::Char = getindex(String(primer), i)
Base.getindex(primer::AbstractPrimer, r::UnitRange{Int}) = getindex(String(primer), r)

Base.convert(::Type{DegenOligo}, primer::AbstractPrimer) = convert(DegenOligo, primer.consensus)

Oligos.n_unique_oligos(primer::AbstractPrimer)::BigInt = n_unique_oligos(primer.consensus)
Oligos.n_deg_pos(primer::AbstractPrimer)::Int = n_deg_pos(primer.consensus)
Oligos.description(primer::AbstractPrimer)::String = description(primer.consensus)
Oligos.hasgaps(::AbstractPrimer)::Bool = false
Oligos.nondegens(primer::AbstractPrimer) = nondegens(primer.consensus)
Oligos.oligo_range(primer::AbstractPrimer)::UnitRange{Int} = primer.pos

"""
    _has_nonspecific_match(primer_seq::AbstractString, msa::AbstractMSA, target_interval; min_identity=0.8) -> Bool

Check if a degenerate primer sequence has high-probability matches outside the target interval in the MSA.
Evaluates both forward and reverse complement orientations.
"""
function _has_nonspecific_match(
    primer_seq::AbstractString,
    msa::AbstractMSA,
    target_interval::UnitRange{Int};
    min_identity::Float64=0.8,
)::Bool
    plen = length(primer_seq)
    L = width(msa)
    plen > L && return false

    p_oligo = Oligos.DegenOligo(primer_seq)
    rc_oligo = _ext_revcomp(p_oligo)
    rc_seq = String(rc_oligo)

    threshold = plen * min_identity

    # Local function to slide a sequence over the MSA base count matrix
    function _check_seq(seq::AbstractString)::Bool
        for startpos in 1:(L-plen+1)
            # Skip if the window overlaps with the target interval
            window_end = startpos + plen - 1
            if startpos <= target_interval.stop && window_end >= target_interval.start
                continue
            end

            match_score = 0.0
            valid = true

            for j in 1:plen
                col_idx = startpos + j - 1
                col_probs = get_base_count(msa, col_idx)

                # Skip windows with significant gaps (depth < 0.5)
                if sum(col_probs) < 0.5
                    valid = false
                    break
                end

                p_char = seq[j]
                p_probs = get(Oligos.IUPAC_PROBS, p_char, (0.0, 0.0, 0.0, 0.0))

                # Accumulate probability of match
                match_score += sum(p_probs .* col_probs)

                # Early exit if max possible score is below threshold
                if match_score + (plen - j) < threshold
                    valid = false
                    break
                end
            end

            if valid && match_score >= threshold
                return true
            end
        end
        return false
    end

    # Check both forward primer binding (reverse complement match)
    # and reverse primer binding (forward match)
    return _check_seq(primer_seq) || _check_seq(rc_seq)
end

"""
    MiniBlastHit

Represents a single hit returned by the `miniblast` function.

# Fields
- `pos::UnitRange{Int}`: The start and end positions of the match in the MSA (1-based).
- `strand::Symbol`: The strand of the match (`:forward` or `:reverse`).
- `identity::Float64`: The average match probability (0.0 to 1.0).
"""
struct MiniBlastHit
    pos::UnitRange{Int}
    strand::Symbol
    identity::Float64
end

"""
    miniblast(target_msa::AbstractMSA, query_oligo::AbstractDegen, threshold=0.75) -> Vector{MiniBlastHit}

Search for high-probability matches of `query_oligo` within `target_msa` using a probabilistic sliding window.
Evaluates both forward and reverse complement orientations of the query.

# Arguments
- `target_msa::AbstractMSA`: The multiple sequence alignment to search within.
- `query_oligo::AbstractDegen`: The degenerate oligo to search for.
- `threshold::Real=0.75`: Minimum average match probability (identity) required to report a hit.

# Returns
- `Vector{MiniBlastHit}`: A list of matches sorted by identity (descending).
"""
function miniblast(
    target_msa::AbstractMSA,
    query_oligo::AbstractDegen,
    threshold::Real=0.75,
)::Vector{MiniBlastHit}
    q_ungapped_str = filter(!=('-'), String(query_oligo))
    plen = length(q_ungapped_str)
    L = width(target_msa)
    (iszero(plen) || plen > L) && return MiniBlastHit[]

    q_oligo_ungapped = DegenOligo(q_ungapped_str)
    rc_str = String(_ext_revcomp(q_oligo_ungapped))

    abs_match_threshold = plen * threshold
    hits = MiniBlastHit[]

    function _search_strand(seq::AbstractString, strand::Symbol)
        for startpos in 1:(L-plen+1)
            match_score = 0.0
            valid = true

            for j in 1:plen
                col_idx = startpos + j - 1
                col_probs = get_base_count(target_msa, col_idx)

                if sum(col_probs) < 0.5
                    valid = false
                    break
                end

                p_char = seq[j]
                p_probs = get(Oligos.IUPAC_PROBS, p_char, (0.0, 0.0, 0.0, 0.0))

                match_score += sum(p_probs .* col_probs)

                if match_score + (plen - j) < abs_match_threshold
                    valid = false
                    break
                end
            end

            if valid && match_score >= abs_match_threshold
                identity = match_score / plen
                push!(hits, MiniBlastHit(startpos:(startpos+plen-1), strand, identity))
            end
        end
    end

    _search_strand(q_ungapped_str, :forward)
    _search_strand(rc_str, :reverse)

    sort!(hits, by=h -> h.identity, rev=true)
    return hits
end

"""
    construct_primers(msa::AbstractMSA; kwargs...) -> Vector{Primer{DegenOligo}}

Construct a list of candidate primers from an MSA based on thermodynamic, conservation, and specificity filters.

# Arguments
- `msa::AbstractMSA`: The multiple sequence alignment.
- `is_forward::Bool=true`: Design forward (`true`) or reverse (`false`) primers.
- `length_range::UnitRange{Int}=17:23`: Allowed primer lengths.
- `tail_length::Int=3`: Length of the 3' tail region.
- `head_degen_pos::Int=5`: Maximum allowed degenerate positions in the 5' head region.
- `tail_degen_pos::Int=0`: Maximum allowed degenerate positions in the 3' tail region.
- `slack::Real=0.05`: Minimum frequency threshold for including a base in the degenerate consensus.
- `gc_range::UnitRange{Int}=40:60`: Allowed GC content percentage range.
- `tm_range::UnitRange{Int}=55:60`: Allowed melting temperature (Tm) range.
- `min_delta_g::Real=-5.0`: Minimum allowed free energy (ΔG) at `dg_temp`.
- `min_msadepth::Float64=0.75`: Minimum sequence depth (coverage) required across the primer region.
- `max_oligo_variants::Int=100`: Maximum number of unique sequences the degenerate primer can represent.
- `max_samples::Int=5000`: Number of samples for Monte Carlo estimation of Tm and ΔG.
- `tm_conf_int::Real=0.2`: Confidence interval for Tm.
- `tm_conds=:pcr`: Thermodynamic conditions for Tm calculation.
- `dg_temp::Real=mean(tm_range)`: Temperature for ΔG calculation.
- `offtarget_reject_threshold::Real=0.75`: Maximum allowed average match probability for off-target binding within the MSA. If a candidate primer matches another region in the MSA with an average probability greater than or equal to this threshold, it is discarded. Checks both forward and reverse complement orientations.
- `nested_pair::Union{Nothing, Tuple{Pair{<:AbstractPrimer, <:AbstractPrimer}, Int}}=nothing`: An optional tuple specifying a flanking primer pair and an offset for nested PCR design.
  - If `nothing` or `offset == 0`, constructs primers across the entire MSA.
  - If `offset < 0`, constructs primers strictly inside the flanking pair's amplicon boundaries, shrunk by the absolute value of the offset.
  - If `offset > 0`, constructs forward primers upstream of the flanking amplicon (with the given offset) and reverse primers downstream.

# Returns
- `Vector{Primer{DegenOligo}}`: A list of valid candidate primers.

See also [`best_pairs`](@ref), [`Primer`](@ref), [`consensus_degen`](@ref).
"""
function construct_primers(
    msa::AbstractMSA;
    is_forward::Bool=true,
    length_range::UnitRange{Int}=17:23,
    tail_length::Int=3,
    tail_degen_pos::Int=0,
    head_degen_pos::Int=5,
    slack::Real=0.05,
    gc_range::UnitRange{Int}=40:60,
    tm_range::UnitRange{Int}=55:60,
    min_delta_g::Real=-5.0,
    min_msadepth::Float64=0.75,
    max_oligo_variants::Int=100,
    max_samples::Int=5000,
    tm_conf_int::Real=0.2,
    tm_conds=:pcr,
    dg_temp::Real=mean(tm_range),
    offtarget_reject_threshold::Real=0.75,
    adapter_pair = GLOBAL_ADAPTERS[],
    max_dg_drop::Real=1.0,
    nested_pair::Union{
        Nothing,
        Tuple{Pair{<:AbstractPrimer,<:AbstractPrimer},Int},
    }=nothing,
)::Vector{Primer{DegenOligo}}
    0 ≤ slack < 1 || throw(ArgumentError("slack must be in [0,1)"))
    0 ≤ min_msadepth ≤ 1 || throw(ArgumentError("min_msadepth must be in [0,1]"))
    1 ≤ max_oligo_variants || throw(ArgumentError("max_oligo_variants must be at least 1"))
    2 ≤ minimum(length_range) ||
        throw(ArgumentError("lower bound of length_range must be ≥ 2nt"))
    0 ≤ tail_length ≤ minimum(length_range) ||
        throw(ArgumentError("tail_length must be [0, length_range.start]"))
    0 ≤ gc_range.start ≤ gc_range.stop ≤ 100 ||
        throw(ArgumentError("gc_range must be in [0, 100]"))
    0 ≤ tm_range.start ≤ tm_range.stop ≤ 100 ||
        throw(ArgumentError("tm_range must be in [0, 100]"))
    0 ≤ offtarget_reject_threshold ≤ 1 ||
        throw(ArgumentError("offtarget_reject_threshold must be in [0,1]"))

    primers = Primer{DegenOligo}[]
    L = length(msa)
    base_count = get_base_count(msa)

    search_start::Int = 1
    search_end::Int = L

    if !isnothing(nested_pair)
        flanking_pair, offset = nested_pair
        amp_start = flanking_pair.first.pos.start
        amp_stop = flanking_pair.second.pos.stop

        1 <= amp_start <= amp_stop <= L ||
            throw(ArgumentError(
                "Flanking pair positions ($amp_start:$amp_stop) are out of MSA bounds (1:$L)",
            ))

        if offset < 0
            search_start = amp_start - offset
            search_end = amp_stop + offset

            if search_start > search_end
                throw(ArgumentError(
                    "Nested offset ($offset) is too large, leaving no valid region between $search_start and $search_end",
                ))
            end
            if search_start < 1 || search_end > L
                throw(ArgumentError(
                    "Nested search bounds ($search_start:$search_end) exceed MSA boundaries (1:$L)",
                ))
            end
            if search_end - search_start + 1 < minimum(length_range)
                throw(ArgumentError(
                    "Nested region length ($(search_end - search_start + 1)) is smaller than minimum primer length ($(minimum(length_range)))",
                ))
            end
        elseif offset > 0
            if is_forward
                search_start = 1
                search_end = amp_start - offset
                if search_end < 1
                    throw(ArgumentError(
                        "Not enough space upstream of the flanking amplicon (amp_start=$amp_start) to construct outer forward primers with offset $offset",
                    ))
                end
                if search_end - search_start + 1 < minimum(length_range)
                    throw(ArgumentError(
                        "Not enough space upstream to construct outer forward primers. Available: $(search_end - search_start + 1)bp, Min primer length: $(minimum(length_range))",
                    ))
                end
            else
                search_start = amp_stop + offset
                search_end = L
                if search_start > L
                    throw(ArgumentError(
                        "Not enough space downstream of the flanking amplicon (amp_stop=$amp_stop) to construct outer reverse primers with offset $offset",
                    ))
                end
                if search_end - search_start + 1 < minimum(length_range)
                    throw(ArgumentError(
                        "Not enough space downstream to construct outer reverse primers. Available: $(search_end - search_start + 1)bp, Min primer length: $(minimum(length_range))",
                    ))
                end
            end
        end
    end

    prog = Progress(
        length(length_range);
        desc=rpad(is_forward ? "Constructing F…" : "Constructing R…", 20),
        color=:white,
        barlen=10,
    )
    l = ReentrantLock()

    Threads.@threads for len in length_range
        len > L && continue
        tail_len = min(tail_length, len)
        head_len = len - tail_len
        head_len < 0 && continue

        actual_search_start = max(1, search_start)
        actual_search_end = min(L - len + 1, search_end - len + 1)

        if actual_search_start > actual_search_end
            continue
        end

        for startpos in actual_search_start:actual_search_end
            interval::UnitRange{Int} = startpos:(startpos+len-1)
            depths = msadepth(msa, interval)
            any(<(min_msadepth), depths) && continue

            if is_forward
                head_interval = startpos:(startpos+head_len-1)
                tail_interval = (startpos+head_len):(startpos+len-1)
            else
                head_interval = (startpos+tail_len):(startpos+len-1)
                tail_interval = startpos:(startpos+tail_len-1)
            end

            if head_len > 0
                head_freqs = @view base_count[:, head_interval]
                head_deg = sum(count(>(slack), col) > 1 for col in eachcol(head_freqs))
                head_deg > head_degen_pos && continue
            end

            if tail_len > 0
                tail_freqs = @view base_count[:, tail_interval]
                tail_deg = sum(count(>(slack), col) > 1 for col in eachcol(tail_freqs))
                tail_deg > tail_degen_pos && continue
            end

            _cons = consensus_degen(msa, interval; slack=slack)
            hasgaps(_cons) && continue

            cons_str = String(parent(_cons))
            if _has_nonspecific_match(
                cons_str,
                msa,
                interval;
                min_identity=offtarget_reject_threshold,
            )
                continue
            end

            gapped_cons = is_forward ? _cons : _ext_revcomp(_cons)

            cons = DegenOligo(gapped_cons)
            n_unique_oligos(cons) > max_oligo_variants && continue

            cons = n_unique_oligos(cons) == 1 ? Oligo(cons) : cons

            gc = _ext_gc_content(cons)
            !(gc_range.start / 100 <= gc <= gc_range.stop / 100) && continue

            dg_val = _ext_dg(cons; max_samples=max_samples, temp=dg_temp)
            dg_val < min_delta_g && continue

            Tm = _ext_tm(
                cons;
                max_samples=max_samples,
                conf_int=tm_conf_int,
                conditions=tm_conds,
            )
            (tm_range.stop < first(Tm.conf) || last(Tm.conf) < tm_range.start) && continue

            final_dg = dg_val
            adapter_to_add = nothing

            if !isnothing(adapter_pair)
                adapter_to_add = is_forward ?
                                 adapter_pair.first :
                                 adapter_pair.second
                full_oligo = adapter_to_add * cons

                final_dg = _ext_dg(full_oligo; max_samples=max_samples, temp=Tm.mean)

                if (final_dg - dg_val) < -max_dg_drop
                    continue
                end
            end

            primer = Primer{DegenOligo}(
                msa,
                interval,
                is_forward,
                DegenOligo(cons),
                tail_len,
                Tm,
                final_dg,
                gc,
                slack,
                adapter_to_add,
            )

            lock(l)
            try
                push!(primers, primer)
                next!(prog)
            finally
                unlock(l)
            end
        end
    end

    return primers
end

"""
    best_pairs(forwards::Vector{<:Primer}, reverses::Vector{<:Primer};
        amplicon_len::UnitRange{Int}=0:9999,
        max_tm_diff::Real=4.0,
        nested_pair=nothing
    ) -> Vector{Pair{Primer{DegenOligo}}}

Find the best matching pairs of forward and reverse primers.

# Arguments
- `forwards::Vector{<:Primer}`: A list of forward primers.
- `reverses::Vector{<:Primer}`: A list of reverse primers.
- `amplicon_len::UnitRange{Int}=0:9999`: Allowed range for the total amplicon length.
- `max_tm_diff::Real=4.0`: Maximum allowed difference in mean Tm between forward and reverse primers.
- `nested_pair::Union{Nothing, Tuple{Pair{<:AbstractPrimer, <:AbstractPrimer}, Int}}=nothing`: An optional tuple specifying a flanking primer pair and an offset for nested PCR design.
  - If `nothing` or `offset == 0`, performs normal pairing.
  - If `offset < 0`, only pairs entirely inside the flanking pair's amplicon (minus the offset margin) are considered.
  - If `offset > 0`, only pairs with the forward primer upstream and reverse primer downstream of the flanking amplicon (plus the offset margin) are considered.

# Returns
- `Vector{Pair{Primer{DegenOligo}, Primer{DegenOligo}}}`: A sorted list of valid primer pairs, ordered by the smallest difference in mean Tm.

See also [`construct_primers`](@ref), [`Primer`](@ref).
"""
function best_pairs(
    forwards::Vector{<:Primer},
    reverses::Vector{<:Primer};
    amplicon_len::UnitRange{Int}=0:9999,
    max_tm_diff::Real=4.0,
    nested_pair::Union{
        Nothing,
        Tuple{Pair{<:AbstractPrimer,<:AbstractPrimer},Int},
    }=nothing,
)::Vector{Pair{Primer{DegenOligo},Primer{DegenOligo}}}
    pairs = Pair{Primer{DegenOligo},Primer{DegenOligo}}[]
    (isempty(forwards) || isempty(reverses)) && return pairs

    all(p -> p.is_forward, forwards) ||
        throw(ArgumentError("All forwards must be forward primers"))
    all(p -> !p.is_forward, reverses) ||
        throw(ArgumentError("All reverses must be reverse primers"))

    anymsa = root(first(forwards).msa)
    all(root(p.msa) == anymsa for p in forwards) ||
        throw(ArgumentError("All primers must refer to the same MSA"))
    all(root(p.msa) == anymsa for p in reverses) ||
        throw(ArgumentError("All primers must refer to the same MSA"))

    L = length(anymsa)
    use_nested = false
    nested_amp_start = 0
    nested_amp_stop = 0
    nested_offset = 0

    if !isnothing(nested_pair)
        flanking_pair, offset = nested_pair

        if root(flanking_pair.first.msa) != anymsa ||
           root(flanking_pair.second.msa) != anymsa
            throw(ArgumentError("Flanking pair must refer to the same MSA as the primers"))
        end

        nested_amp_start = flanking_pair.first.pos.start
        nested_amp_stop = flanking_pair.second.pos.stop
        nested_offset = offset

        1 <= nested_amp_start <= nested_amp_stop <= L ||
            throw(ArgumentError(
                "Flanking pair positions ($nested_amp_start:$nested_amp_stop) are out of MSA bounds (1:$L)",
            ))

        if nested_offset < 0
            valid_start = nested_amp_start - nested_offset
            valid_end = nested_amp_stop + nested_offset

            valid_start > valid_end &&
                throw(ArgumentError(
                    "Nested offset ($nested_offset) is too large, leaving no valid region between $valid_start and $valid_end",
                ))
            (valid_start < 1 || valid_end > L) &&
                throw(ArgumentError(
                    "Nested search bounds ($valid_start:$valid_end) exceed MSA boundaries (1:$L)",
                ))

            if minimum(amplicon_len) > valid_end - valid_start + 1
                throw(ArgumentError(
                    "Specified minimum amplicon length ($(minimum(amplicon_len))) is larger than the available nested region ($(valid_end - valid_start + 1)bp)",
                ))
            end
        elseif nested_offset > 0
            min_theoretical_amplicon = (nested_amp_stop + nested_offset) -
                                       (nested_amp_start - nested_offset) +
                                       1
            if maximum(amplicon_len) < min_theoretical_amplicon
                throw(ArgumentError(
                    "Specified amplicon_len max ($(maximum(amplicon_len))) is smaller than the minimum possible outer amplicon length ($min_theoretical_amplicon)",
                ))
            end
        end

        use_nested = true
    end

    @showprogress desc=rpad("Paring…", 20) enabled=(length(forwards)>1000) barlen=10 for f in forwards
        for r in reverses
            if f.pos.stop >= r.pos.start
                # overlapping primers
                continue
            end
            if max_tm_diff < abs(f.tm.mean - r.tm.mean)
                # too big melting T difference
                continue
            end

            amplicon = r.pos.stop - f.pos.start + 1
            if !(amplicon in amplicon_len)
                continue
            end

            if use_nested
                if nested_offset < 0
                    if f.pos.start < nested_amp_start - nested_offset ||
                       r.pos.stop > nested_amp_stop + nested_offset
                        continue
                    end
                elseif nested_offset > 0
                    if f.pos.stop >= nested_amp_start - nested_offset ||
                       r.pos.start <= nested_amp_stop + nested_offset
                        continue
                    end
                end
            end

            push!(pairs, f => r)
        end
    end

    sort!(pairs; by=p -> abs(p.first.tm.mean - p.second.tm.mean))
    return pairs
end

"""
    add_illumina_adapters(pair, fwd_adapter, rev_adapter; max_dg_drop) -> Union{Nothing, Pair{Primer}}
    add_illumina_adapters(pairs, fwd_adapter, rev_adapter; max_dg_drop::Real=1.0) -> Vector{Pair{Primer}}

Concatenates Illumina adapters to the 5' ends of primers and filters out candidates where 
the adapter significantly worsens the thermodynamic stability (ΔG) of either primer.

# Arguments
- `pair` / `pairs`: A single primer pair or a vector of primer pairs.
- `fwd_adapter::AbstractString`: The 5' adapter sequence for the forward primer.
- `rev_adapter::AbstractString`: The 5' adapter sequence for the reverse primer.
- `max_dg_drop::Real=1.0`: The maximum allowed decrease in ΔG (in kcal/mol, more negative) 
  caused by the addition of the adapter. If `ΔG_full - ΔG_primer < -max_dg_drop`, 
  the pair is discarded.

# Returns
- For a single `Pair`: Returns a `Pair{Primer{DegenOligo}, Primer{DegenOligo}}` if valid, otherwise `nothing`.
- For an `AbstractVector` of `Pair`: Returns a `Vector{Pair{Primer{DegenOligo}, Primer{DegenOligo}}}`.

# Throws
- `ArgumentError`: If the adapter sequences contain invalid or degenerate characters.

# Details
Thermodynamic stability (ΔG) of the full sequence is calculated at the primer's mean Tm. 
The vectorized method optimizes performance by caching ΔG calculations for primers that 
appear in multiple pairs.
"""
function add_illumina_adapters(
    pair::Pair{<:AbstractPrimer,<:AbstractPrimer},
    fwd_adapter::AbstractString,
    rev_adapter::AbstractString;
    max_dg_drop::Real=1.0,
)::Union{Nothing,Pair{Primer{DegenOligo},Primer{DegenOligo}}}
    res = add_illumina_adapters([pair], fwd_adapter, rev_adapter; max_dg_drop=max_dg_drop)
    return isempty(res) ? nothing : first(res)
end

function add_illumina_adapters(
    pairs::AbstractVector{<:Pair{<:AbstractPrimer,<:AbstractPrimer}},
    fwd_adapter::AbstractString,
    rev_adapter::AbstractString;
    max_dg_drop::Real=1.0,
)::Vector{Pair{Primer{DegenOligo},Primer{DegenOligo}}}

    fwd_ad = Oligo(fwd_adapter, "Illumina_Fwd_Adapter")
    rev_ad = Oligo(rev_adapter, "Illumina_Rev_Adapter")

    cache = IdDict{AbstractPrimer,Union{Nothing,Primer{DegenOligo}}}()

    function _process_primer(p::AbstractPrimer, adapter::Oligo)::Union{Nothing,Primer{DegenOligo}}
        if haskey(cache, p)
            return cache[p]
        end

        temp = p.tm.mean
        full_oligo = adapter * p.consensus
        dg_full = _ext_dg(full_oligo; max_samples=1000, temp=temp)

        if (dg_full - p.dg) < -max_dg_drop
            cache[p] = nothing
            return nothing
        end

        new_p = Primer{DegenOligo}(
            p.msa,
            p.pos,
            p.is_forward,
            p.consensus,
            p.tail_length,
            p.tm,
            dg_full,
            p.gc,
            p.slack,
            adapter,
        )

        cache[p] = new_p
        return new_p
    end

    full_pairs = Pair{Primer{DegenOligo},Primer{DegenOligo}}[]

    for pp in pairs
        f_new = _process_primer(pp.first, fwd_ad)
        isnothing(f_new) && continue

        r_new = _process_primer(pp.second, rev_ad)
        isnothing(r_new) && continue

        push!(full_pairs, f_new => r_new)
    end

    return full_pairs
end

function Base.convert(::Type{Primer{T1}}, p::Primer{T2})::Primer{T1} where {T1<:Union{Oligo,DegenOligo},T2<:Union{Oligo,DegenOligo}}
    T2 === T1 && return p
    return Primer{T1}(
        p.msa,
        p.pos,
        p.is_forward,
        convert(T1, p.consensus),
        p.tail_length,
        p.tm,
        p.dg,
        p.gc,
        p.slack,
        p.adapter,
    )
end

function Base.convert(
    ::Type{Pair{Primer{T3},Primer{T3}}},
    p::Pair{Primer{T1},Primer{T2}},
)::Pair{Primer{T3},Primer{T3}} where {T1<:Union{Oligo,DegenOligo},T2<:Union{Oligo,DegenOligo},T3<:Union{Oligo,DegenOligo}}
    T1 === T3 && T2 === T3 && return p
    return Pair(convert(Primer{T3}, p.first), convert(Primer{T3}, p.second))
end

function Base.convert(
    ::Type{Pair{T1,T2}},
    p::Pair{<:AbstractPrimer,<:AbstractPrimer},
)::Pair{T1,T2} where {T1<:AbstractPrimer,T2<:AbstractPrimer}
    return convert(T1, p.first) => convert(T2, p.second)
end

include("show_primers.jl")
include("export_primers.jl")

end # module
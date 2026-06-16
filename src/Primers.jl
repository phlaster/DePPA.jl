module Primers

export AbstractPrimer
export Primer
export construct_primers, best_pairs

using ..Oligos
using ..Alignments

using ProgressMeter
using Statistics

"""
    AbstractPrimer{T<:Union{Oligo,DegenOligo}}

Represent the abstract supertype for PCR primers.

See also [`Primer`](@ref).
"""
abstract type AbstractPrimer{T<:Union{Oligo,DegenOligo}} end

"""
    Primer{T} <: AbstractPrimer{T}

Represent a concrete PCR primer. Store the consensus sequence, its position in the MSA,
and thermodynamic properties (Tm, dG, GC content).

See also [`AbstractPrimer`](@ref), [`construct_primers`](@ref).
"""
struct Primer{T} <: AbstractPrimer{T}
    msa::AbstractMSA
    pos::UnitRange{Int}
    is_forward::Bool
    consensus::T
    tail_length::Int
    tm::@NamedTuple{mean::Float64, conf::Tuple{Float64, Float64}, min::Float64, max::Float64}
    dg::Float64
    gc::Float64
    slack::Float64
end

"""
    Primer(msa::AbstractMSA, interval::UnitRange{Int}; kwargs...)

Construct a `Primer` object for a given interval in the MSA, calculating its thermodynamic properties.

# Arguments
- `msa::AbstractMSA`: The multiple sequence alignment.
- `interval::UnitRange{Int}`: The position range of the primer in the MSA.
- `is_forward::Bool=true`: Design a forward (`true`) or reverse (`false`) primer.
- `tail_length::Int=3`: Length of the 3' tail region.
- `max_samples::Int=1000`: Number of samples for Monte Carlo estimation of Tm and dG.
- `tm_conf_int=0.8`: Confidence interval for Tm.
- `tm_conds=:pcr`: Thermodynamic conditions for Tm calculation.
- `dg_temp=37.0`: Temperature for dG calculation.
- `slack=0.0`: Minimum frequency threshold for including a base in the degenerate consensus.
- `descr`: Description string for the primer.

See also [`construct_primers`](@ref), [`consensus_degen`](@ref).
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
    descr="Primer for $(nseqs(msa)) seq MSA at positions $interval"
)
    _cons = consensus_degen(msa, interval; slack=slack)
    gapped_cons = is_forward ? _cons : _ext_revcomp(_cons)
    
    underlying_oligo = DegenOligo(String(gapped_cons), string(descr))
    
    Tm = _ext_tm(underlying_oligo; max_samples=max_samples, conf_int=tm_conf_int, conditions=tm_conds)
    dG = _ext_dg(underlying_oligo; max_samples=max_samples, temp=dg_temp)
    GC = _ext_gc_content(underlying_oligo)
    Primer(msa, interval, is_forward, underlying_oligo, tail_length, Tm, dG, GC, slack)
end

# These are overloaded in ext/SeqFoldExt.jl to load SeqFold.jl library dynamically
_ext_revcomp(args...; kwargs...) = error(
    "Primer construction requires SeqFold library to be loaded.\n" *
    "In order to get this functionality, please `]add SeqFold` to your project\n" *
    "and load it with `using SeqFold` before constructing primers.")
_ext_tm(args...; kwargs...) = error(
    "Primer construction requires SeqFold library to be loaded.\n" *
    "In order to get this functionality, please `]add SeqFold` to your project\n" *
    "and load it with `using SeqFold` before constructing primers.")
_ext_dg(args...; kwargs...) = error(
    "Primer construction requires SeqFold library to be loaded.\n" *
    "In order to get this functionality, please `]add SeqFold` to your project\n" *
    "and load it with `using SeqFold` before constructing primers.")
_ext_gc_content(args...; kwargs...) = error(
    "Primer construction requires SeqFold library to be loaded.\n" *
    "In order to get this functionality, please `]add SeqFold` to your project\n" *
    "and load it with `using SeqFold` before constructing primers."
)

Base.String(primer::AbstractPrimer) = String(primer.consensus)
Base.length(primer::AbstractPrimer) = length(primer.consensus)
Base.isempty(primer::AbstractPrimer) = isempty(primer.consensus)
Base.iterate(primer::AbstractPrimer, state...) = iterate(primer.consensus, state...)
Base.getindex(primer::AbstractPrimer, i::Int) = getindex(primer.consensus, i)
Base.getindex(primer::AbstractPrimer, r::UnitRange{Int}) = getindex(primer.consensus, r)

Base.convert(::Type{DegenOligo}, primer::AbstractPrimer) = primer.consensus

Oligos.n_unique_oligos(primer::AbstractPrimer) = n_unique_oligos(primer.consensus)
Oligos.n_deg_pos(primer::AbstractPrimer) = n_deg_pos(primer.consensus)
Oligos.description(primer::AbstractPrimer) = description(primer.consensus)
Oligos.hasgaps(::AbstractPrimer) = false
Oligos.nondegens(primer::AbstractPrimer) = nondegens(primer.consensus)
Oligos.oligo_range(primer::AbstractPrimer) = primer.pos

"""
    construct_primers(msa::AbstractMSA; kwargs...) -> Vector{Primer{DegenOligo}}

Construct a list of candidate primers from an MSA based on thermodynamic and conservation filters.

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
- `min_delta_g::Real=-5.0`: Minimum allowed free energy (dG) at `dg_temp`.
- `min_msadepth::Float64=0.75`: Minimum sequence depth (coverage) required across the primer region.
- `max_oligo_variants::Int=100`: Maximum number of unique sequences the degenerate primer can represent.
- `max_samples::Int=5000`: Number of samples for Monte Carlo estimation of Tm and dG.
- `tm_conf_int::Real=0.2`: Confidence interval for Tm.
- `tm_conds=:pcr`: Thermodynamic conditions for Tm calculation.
- `dg_temp::Real=mean(tm_range)`: Temperature for dG calculation.

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
    dg_temp::Real=mean(tm_range)
)::Vector{Primer{DegenOligo}}
    0 ≤ slack < 1 || throw(ArgumentError("slack must be in [0,1)"))
    0 ≤ min_msadepth ≤ 1 || throw(ArgumentError("min_msadepth must be in [0,1]"))
    1 ≤ max_oligo_variants || throw(ArgumentError("max_oligo_variants must be at least 1"))
    2 ≤ minimum(length_range) || throw(ArgumentError("lower bound of length_range must be ≥ 2nt"))
    0 ≤ tail_length ≤ minimum(length_range) || throw(ArgumentError("tail_length must be [0, length_range.start]"))
    0 ≤ gc_range.start ≤ gc_range.stop ≤ 100 || throw(ArgumentError("gc_range must be in [0, 100]"))
    0 ≤ tm_range.start ≤ tm_range.stop ≤ 100 || throw(ArgumentError("tm_range must be in [0, 100]"))
    
    primers = Primer{DegenOligo}[]
    L = length(msa)
    base_count = get_base_count(msa)
    prog = Progress(length(length_range);
        desc=is_forward ? "Constructing F..." : "Constructing R...",
        color=:white,
        barlen=10
    )
    l = ReentrantLock()
    
    Threads.@threads for len in length_range
        len > L && continue
        tail_len = min(tail_length, len)
        head_len = len - tail_len
        head_len < 0 && continue

        for startpos in 1:(L - len + 1)
            interval::UnitRange{Int} = startpos:(startpos + len - 1)
            depths = msadepth(msa, interval)
            any(<(min_msadepth), depths) && continue
            
            if is_forward
                head_interval = startpos:(startpos + head_len - 1)
                tail_interval = (startpos + head_len):(startpos + len - 1)
            else
                head_interval = (startpos + tail_len):(startpos + len - 1)
                tail_interval = startpos:(startpos + tail_len - 1)
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
            gapped_cons = is_forward ? _cons : _ext_revcomp(_cons)
            hasgaps(gapped_cons) && continue
            
            cons = DegenOligo(gapped_cons)
            n_unique_oligos(cons) > max_oligo_variants && continue

            cons = n_unique_oligos(cons) == 1 ? Oligo(cons) : cons
            
            gc = _ext_gc_content(cons)
            !(gc_range.start / 100 <= gc <= gc_range.stop / 100) && continue
            
            dg_val = _ext_dg(cons; max_samples=max_samples, temp=dg_temp)
            dg_val < min_delta_g && continue
            
            Tm = _ext_tm(cons; max_samples=max_samples, conf_int=tm_conf_int, conditions=tm_conds)
            (tm_range.stop < first(Tm.conf) || last(Tm.conf) < tm_range.start) && continue
            
            primer = Primer{DegenOligo}(msa, interval, is_forward, cons, tail_len, Tm, dg_val, gc, slack)

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
    best_pairs(forwards::Vector{<:Primer}, reverses::Vector{<:Primer}; amplicon_len::UnitRange{Int}=0:9999, max_tm_diff::Real=4.0) -> Vector{Pair{Primer{DegenOligo}}}

Find the best matching pairs of forward and reverse primers.

# Arguments
- `forwards::Vector{<:Primer}`: A list of forward primers.
- `reverses::Vector{<:Primer}`: A list of reverse primers.
- `amplicon_len::UnitRange{Int}=0:9999`: Allowed range for the total amplicon length.
- `max_tm_diff::Real=4.0`: Maximum allowed difference in mean Tm between forward and reverse primers.

# Returns
- `Vector{Pair{Primer{DegenOligo}}}`: A sorted list of valid primer pairs, ordered by the smallest difference in mean Tm.

See also [`construct_primers`](@ref), [`Primer`](@ref).
"""
function best_pairs(
    forwards::Vector{<:Primer},
    reverses::Vector{<:Primer};
    amplicon_len::UnitRange{Int}=0:9999,
    max_tm_diff::Real=4.0
)::Vector{Pair{Primer{DegenOligo}}}
    pairs = Pair{Primer{DegenOligo}}[]
    (isempty(forwards) || isempty(reverses)) && return pairs

    all(p -> p.is_forward, forwards)  || throw(ArgumentError("All forwards must be forward primers"))
    all(p -> !p.is_forward, reverses) || throw(ArgumentError("All reverses must be reverse primers"))
    
    anymsa = root(first(forwards).msa)
    all(root(p.msa) == anymsa for p in forwards) || throw(ArgumentError("All primers must refer to the same MSA"))
    all(root(p.msa) == anymsa for p in reverses) || throw(ArgumentError("All primers must refer to the same MSA"))
    
    @showprogress desc="Matching primer pairs..." enabled=(length(forwards)>1000) barlen=10 for f in forwards
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
            if amplicon in amplicon_len
                push!(pairs, f => r)
            end
        end
    end
    
    sort!(pairs; by = p -> abs(p.first.tm.mean - p.second.tm.mean))
    return pairs
end

function Base.convert(::Type{Primer{T1}}, p::Primer{T2}) where {T1, T2}
    T2 === T1 && return p
    return Primer(
        p.msa, p.pos, p.is_forward,
        convert(T1, p.consensus),
        p.tail_length, p.tm, p.dg, p.gc, p.slack
    )
end

function Base.convert(::Type{Pair{Primer{T3}, Primer{T3}}}, p::Pair{Primer{T1}, Primer{T2}}) where {T1, T2, T3}
    T1 === T3 && T2 === T3 && return p
    return Pair(convert(Primer{T3}, p.first), convert(Primer{T3}, p.second))
end

Base.convert(::Type{Pair{Primer}}, p::Pair{<:AbstractPrimer, <:AbstractPrimer}) = p

include("show_primers.jl")

end # module
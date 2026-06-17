module Oligos
include("utils.jl")

using Statistics

export AbstractOligo, AbstractDegen, AbstractGapped
export Oligo, DegenOligo, GappedOligo, OligoView
export oligo_range, description
export NonDegenIterator, nondegens
export hasgaps, getgaps, n_deg_pos, n_unique_oligos
export sampleChar, sampleView, sampleNondeg, sample_max_gc, sample_min_gc
export unfolded_proportion

"""
    AbstractOligo <: AbstractString

Represent the abstract supertype for all oligomer types in DePPA.
"""
abstract type AbstractOligo <: AbstractString end

"""
    AbstractDegen <: AbstractOligo

Represent the abstract supertype for oligomers that may contain degenerate (IUPAC) bases.
"""
abstract type AbstractDegen <: AbstractOligo end

"""
    AbstractGapped <: AbstractDegen

Represent the abstract supertype for oligomers that may contain gap characters (`-`).
"""
abstract type AbstractGapped <: AbstractDegen end

"""
    description(oligo::AbstractOligo) -> String

Return the description string associated with the oligomer.
"""
description(::AbstractString) = ""

"""
    Oligo(seq::AbstractString, descr::Union{AbstractString,Integer}="")

Represent a non-degenerate nucleic acid sequence (only `A`, `C`, `G`, `T`).

Sequences are automatically converted to uppercase.
"""
struct Oligo <: AbstractOligo
    seq::String
    description::String

    function Oligo(seq::AbstractString, descr::Union{AbstractString,Integer})
        if !isempty(seq) && !isvalid(Oligo, seq)
            invalid_chars = join(setdiff(Set(seq), NON_DEGEN_BASES), ", ")
            error("Oligo contains unallowed characters: $invalid_chars")
        end
        
        seq = uppercase(seq)
        descr_str = string(descr)
        
        new(seq, descr_str)
    end
end

"""
    DegenOligo(seq::AbstractString, descr::Union{AbstractString,Integer}="")

Represent a degenerate nucleic acid sequence, allowing IUPAC ambiguity codes.
"""
struct DegenOligo <: AbstractDegen
    seq::String
    n_deg_pos::Int
    n_unique_oligos::BigInt
    description::String

    function DegenOligo(seq::AbstractString, descr::Union{AbstractString,Integer})
        if !isempty(seq) && !isvalid(DegenOligo, seq)
            invalid_chars = join(setdiff(Set(seq), ALL_BASES), ", ")
            error("DegenOligo contains unallowed characters: $invalid_chars")
        end
        
        seq = uppercase(seq)
        descr_str = string(descr)

        n_deg = count(char -> char in DEGEN_BASES, seq)
        n_unique = reduce(*, (IUPAC_COUNTS[char] for char in seq), init=BigInt(1))
        
        new(seq, n_deg, n_unique, descr_str)
    end
end

"""
    GappedOligo(seq::AbstractString, descr::Union{AbstractString,Integer}="")

Represent a gapped nucleic acid sequence, allowing gap characters (`-`).
"""
struct GappedOligo <: AbstractGapped
    parent::DegenOligo
    gaps::Vector{Pair{Int, Int}}
    total_length::Int

    function GappedOligo(seq::AbstractString, descr::Union{AbstractString,Integer})
        if !isempty(seq) && !isvalid(GappedOligo, seq)
            invalid_chars = join(setdiff(Set(seq), BASES_W_GAPS), ", ")
            error("DegenOligo contains unallowed characters: $invalid_chars")
        end
        
        seq = uppercase(seq)
        total_len = length(seq)
        descr_str = string(descr)
        parent_seq = filter(!=('-'), seq)
        parent_oligo = DegenOligo(parent_seq, descr_str)
        
        gaps = Pair{Int,Int}[]
        parent_pos = 0
        i = 1
        
        @inbounds while i <= total_len
            if seq[i] != '-'
                parent_pos += 1
                i += 1
            else
                gap_start = parent_pos + 1
                len = 0
                while i <= total_len && seq[i] == '-'
                    len += 1
                    i += 1
                end
                push!(gaps, gap_start => len)
            end
        end

        return new(parent_oligo, gaps, total_len)
    end
end

"""
    OligoView{T<:AbstractOligo} <: AbstractOligo

Represent a lightweight view into an [`AbstractOligo`](@ref) as a contiguous subsequence.
"""
struct OligoView{T<:AbstractOligo} <: AbstractOligo
    parent::T
    range::UnitRange{Int}

    function OligoView(oligo::T, interval::UnitRange{Int}) where T<:Union{Oligo, DegenOligo, GappedOligo}
        @boundscheck checkbounds(oligo, interval)
        new{T}(oligo, interval)
    end
end

"""
    NonDegenIterator{T<:AbstractOligo}

Iterate over all unique non-degenerate sequences represented by a degenerate oligomer.
"""
struct NonDegenIterator{T<:AbstractOligo}
    oligo::T
    n_variants::Integer
    NonDegenIterator(o::T) where T<:AbstractOligo = new{T}(o, n_unique_oligos(o))
end


##########################
#      Base methods      #
##########################
#### NonDegenIterator ####
function Base.iterate(iter::NonDegenIterator)
    oligo = parent(iter)
    length(oligo) == 0 && return nothing
    
    options = [IUPAC_B2V[c] for c in String(oligo)]
    lens = length.(options)
    indices = ones(Int, length(oligo))
    buffer = Vector{Char}(undef, length(oligo))
    
    for j in 1:length(oligo)
        buffer[j] = options[j][indices[j]]
    end
    
    state = (indices, options, lens, buffer, length(oligo))
    descr = string("Non-degen sample from: ", description(oligo))
    return Oligo(String(buffer), descr), state
end
function Base.iterate(iter::NonDegenIterator, state)
    indices, options, lens, buffer, n = state
    n == 0 && return nothing

    pos = n
    @inbounds while pos > 0
        indices[pos] += 1
        indices[pos] <= lens[pos] && break
        indices[pos] = 1
        pos -= 1
    end
    pos == 0 && return nothing

    @inbounds @simd for j in 1:n
        buffer[j] = options[j][indices[j]]
    end
    d = description(parent(iter))
    descr = isempty(d) ? "Non-degen sample" : "Non-degen sample from: $d"
    return Oligo(String(buffer), descr), (indices, options, lens, buffer, n)
end
Base.length(iter::NonDegenIterator) = iter.n_variants
Base.eltype(::Type{<:NonDegenIterator}) = Oligo
Base.parent(ndi::NonDegenIterator) = ndi.oligo




######### Oligos ##########
Base.:(==)(o::AbstractOligo, s::SubString{<:Base.AnnotatedString}) = String(o) == String(s)
Base.:(==)(s::SubString{<:Base.AnnotatedString}, o::AbstractOligo) = o == s
Base.:(==)(o::AbstractOligo, s::Base.AnnotatedString) = String(o) == String(s)
Base.:(==)(s::Base.AnnotatedString, o::AbstractOligo) = o == s
Base.:(==)(a::GappedOligo, b::GappedOligo) = (parent(a) == parent(b)) && (getgaps(a) == getgaps(b)) && (length(a) == length(b))

Base.getindex(oligo::AbstractOligo, r::UnitRange{Int}) = OligoView(parent(oligo), r)
Base.getindex(oligo::AbstractOligo, i::Int) = getindex(String(oligo), i)
function Base.getindex(ov::OligoView, r::UnitRange{Int})
    @boundscheck checkbounds(ov, r)
    ov_range_start::Int = oligo_range(ov).start
    actual_start = ov_range_start + r.start - 1
    actual_stop = ov_range_start + r.stop - 1
    actual_range = actual_start:actual_stop
    return OligoView(parent(ov), actual_range)
end
function Base.getindex(oligo::GappedOligo, r::UnitRange{Int})
    @boundscheck begin
        if !(1 <= first(r) <= last(r) <= length(oligo))
            throw(BoundsError(oligo, r))
        end
    end
    return OligoView(oligo, r)
end
function Base.getindex(ov::OligoView, i::Int)
    adjusted_i = oligo_range(ov).start + i - 1
    @boundscheck checkbounds(ov, i)
    return getindex(parent(ov), adjusted_i)
end

Base.uppercase(oligo::AbstractOligo) = String(oligo)
Base.String(oligo::AbstractOligo) = oligo.seq
Base.String(ov::OligoView) = String(parent(ov))[oligo_range(ov)]
function Base.String(go::GappedOligo)
    parent_len = length(parent(go))
    buffer = Vector{Char}(undef, length(go))
    ungapped_pos = 1 
    gap_idx = 1
    buffer_idx = 1
    gaps = getgaps(go)
    L = length(go)

    @inbounds while buffer_idx <= L && ungapped_pos <= parent_len
        if gap_idx <= length(gaps) && ungapped_pos == gaps[gap_idx].first
            _, len = gaps[gap_idx]
            for i in buffer_idx:(buffer_idx + len - 1)
                buffer[i] = '-'
            end
            buffer_idx += len
            gap_idx += 1
            continue
        end
        buffer[buffer_idx] = parent(go).seq[ungapped_pos]
        ungapped_pos += 1
        buffer_idx += 1
    end
    
    @inbounds while buffer_idx <= L && gap_idx <= length(gaps)
        _, len = gaps[gap_idx]
        for i in buffer_idx:(buffer_idx + len - 1)
            buffer[i] = '-'
        end
        buffer_idx += len
        gap_idx += 1
    end
    
    if ungapped_pos <= parent_len || buffer_idx <= L
        throw(ArgumentError("Mismatch in sequence construction: ungapped_pos=$ungapped_pos (expected $(parent_len + 1)), buffer_idx=$buffer_idx (expected $(L + 1))"))
    end
    return String(buffer)
end

Base.iterate(go::GappedOligo) = length(go) == 0 ? nothing : iterate(go, (1, 1, 1, 0))
function Base.iterate(go::GappedOligo, state::NTuple{4, Int})
    pos, ungapped_pos, gap_idx, remaining = state
    if pos > length(go)
        return nothing
    end
    gaps = getgaps(go)
    if remaining > 0
        char = '-'
        new_remaining = remaining - 1
        new_ungapped = ungapped_pos
        new_gap_idx = new_remaining == 0 ? gap_idx + 1 : gap_idx
    elseif gap_idx <= length(gaps) && ungapped_pos == gaps[gap_idx].first
        _, len = gaps[gap_idx]
        char = '-'
        new_remaining = len - 1
        new_ungapped = ungapped_pos
        new_gap_idx = new_remaining == 0 ? gap_idx + 1 : gap_idx
    else
        if ungapped_pos > length(parent(go))
            throw(BoundsError(parent(go), ungapped_pos))
        end
        char = parent(go)[ungapped_pos]
        new_ungapped = ungapped_pos + 1
        new_remaining = 0
        new_gap_idx = gap_idx
    end
    return char, (pos + 1, new_ungapped, new_gap_idx, new_remaining)
end
Base.iterate(oligo::AbstractOligo) = length(oligo) == 0 ? nothing : (oligo[1], 2)
Base.iterate(oligo::AbstractOligo, i::Int) = i>length(oligo) ? nothing : (oligo[i],i+1)

Base.parent(oligo::AbstractOligo) = oligo
Base.parent(oligo::OligoView) = oligo.parent
Base.parent(oligo::GappedOligo) = oligo.parent

Base.ncodeunits(oligo::AbstractOligo) = length(oligo)

Base.codeunit(::AbstractOligo) = UInt8
Base.codeunit(oligo::AbstractOligo, i::Integer) = codeunit(String(oligo), i)
Base.codeunit(ov::OligoView, i::Integer) = codeunit(parent(ov), oligo_range(ov).start + i - 1)
function Base.codeunit(go::GappedOligo, i::Integer)
    @boundscheck 1 <= i <= length(go) || throw(BoundsError(go, i))
    parent_len = length(parent(go))
    ungapped_pos = 1
    gapped_pos = 0
    
    for (start, len) in go.gaps
        num_parent_chars_before_gap = (start - 1) - (ungapped_pos - 1)
        
        if gapped_pos < i <= gapped_pos + num_parent_chars_before_gap
            offset = i - gapped_pos
            return UInt8(parent(go)[ungapped_pos + offset - 1])
        end
        
        gapped_pos += num_parent_chars_before_gap
        ungapped_pos += num_parent_chars_before_gap
        
        if gapped_pos < i <= gapped_pos + len
            return UInt8('-')
        end
        gapped_pos += len
    end
    
    offset = i - gapped_pos
    if ungapped_pos + offset - 1 > parent_len
        throw(BoundsError(go, i))
    end
    return UInt8(parent(go)[ungapped_pos + offset - 1])
end

Base.isvalid(::AbstractOligo) = true
Base.isvalid(oligo::AbstractOligo, i::Int) = 1 <= i <= length(oligo)
Base.isvalid(::Type{Oligo}, s::AbstractString) = all(c -> uppercase(c) in NON_DEGEN_BASES, s)
Base.isvalid(::Type{DegenOligo}, s::AbstractString) = all(c -> uppercase(c) in ALL_BASES, s)
Base.isvalid(::Type{GappedOligo}, s::AbstractString) = all(c -> uppercase(c) in BASES_W_GAPS, s)

Base.length(oligo::AbstractOligo) = length(String(oligo))
Base.length(go::GappedOligo) = go.total_length
Base.length(ov::OligoView) = length(oligo_range(ov))
Base.isempty(oligo::AbstractOligo) = length(oligo) == 0
Base.lastindex(oligo::AbstractOligo) = length(oligo)

function Base.convert(::Type{T}, o::AbstractOligo) where T<:AbstractOligo
    if T === typeof(o)
        return o
    elseif T <: AbstractGapped
        return GappedOligo(String(o), description(o))
    elseif T <: AbstractDegen
        return DegenOligo(String(o), description(o))
    else T <: AbstractOligo
        return Oligo(String(o), description(o))
    end
end

_is_degenerate_type(::Type{T}) where {T<:AbstractOligo} = false
_is_degenerate_type(::Type{T}) where {T<:AbstractDegen} = true
_is_degenerate_type(::Type{OligoView{T}}) where {T} = _is_degenerate_type(T)
_has_gaps_type(::Type{T}) where {T<:AbstractOligo} = T <: Union{GappedOligo, OligoView{<:GappedOligo}}

function Base.promote_rule(::Type{T1}, ::Type{T2}) where {T1<:AbstractOligo, T2<:AbstractOligo}
    is_deg_type = _is_degenerate_type(T1) || _is_degenerate_type(T2)
    is_gapped_type = _has_gaps_type(T1) || _has_gaps_type(T2)
    return is_gapped_type ? GappedOligo : is_deg_type ? DegenOligo : Oligo
end
Base.promote_rule(::Type{<:Union{String, SubString}}, ::Type{T}) where {T<:AbstractOligo} = T
Base.promote_rule(::Type{T}, ::Type{<:Union{String, SubString}}) where {T<:AbstractOligo} = T


const EMPTY_OLIG = Oligo("", "")
Oligo() = EMPTY_OLIG
Oligo(oligo::Oligo) = oligo
Oligo(seq::AbstractString) = Oligo(seq, description(seq))

const EMPTY_DEGENERATE = DegenOligo("", "")
DegenOligo() = EMPTY_DEGENERATE
DegenOligo(oligo::DegenOligo) = oligo
DegenOligo(seq::AbstractString) = DegenOligo(seq, description(seq))

const EMPTY_GAPPED = GappedOligo(DegenOligo(), "")
GappedOligo() = EMPTY_GAPPED
GappedOligo(oligo::GappedOligo) = oligo
GappedOligo(seq::AbstractString) = GappedOligo(seq, description(seq))
DegenOligo(go::GappedOligo) = DegenOligo(String(go), description(go))



##########################
#      Type methods      #
##########################
"""
    oligo_range(oligo::AbstractOligo) -> UnitRange{Int}

Return the index range of the oligomer or view.
"""
oligo_range(oligo::AbstractOligo) = 1:length(oligo)
oligo_range(ov::OligoView) = ov.range

description(oligo::AbstractOligo) = parent(oligo).description
description(ov::OligoView{<:GappedOligo}) = parent(parent(ov)).description

"""
    hasgaps(oligo::AbstractOligo) -> Bool

Return `true` if the sequence contains gap characters (`-`), and `false` otherwise.

See also [`getgaps`](@ref), [`GappedOligo`](@ref).
"""
hasgaps(::AbstractOligo) = false
hasgaps(ov::OligoView{Oligo}) = false
hasgaps(ov::OligoView{DegenOligo}) = false
hasgaps(go::GappedOligo) = !isempty(go.gaps)
hasgaps(ov::OligoView) = any(c == '-' for c in ov)

const EMPTY_GAPS_VECTOR = Vector{Pair{Int, Int}}()

"""
    getgaps(oligo::AbstractOligo) -> Vector{Pair{Int, Int}}

Return a vector of gap locations and lengths for a gapped oligomer.

See also [`hasgaps`](@ref), [`GappedOligo`](@ref).
"""
getgaps(::AbstractOligo) = EMPTY_GAPS_VECTOR
getgaps(go::GappedOligo) = go.gaps
getgaps(ov::OligoView{GappedOligo}) = ov |> GappedOligo |> getgaps

"""
    n_unique_oligos(oligo::AbstractOligo) -> BigInt

Return the total number of unique non-degenerate sequences represented by the oligomer.

See also [`n_deg_pos`](@ref), [`nondegens`](@ref).
"""
n_unique_oligos(::AbstractOligo) = BigInt(1)
n_unique_oligos(d::DegenOligo) = d.n_unique_oligos
n_unique_oligos(ov::OligoView) = reduce(*, (IUPAC_COUNTS[base] for base in ov), init=BigInt(1))
n_unique_oligos(go::GappedOligo) = n_unique_oligos(parent(go))

"""
    n_deg_pos(oligo::AbstractOligo) -> Int

Return the number of degenerate (ambiguity) positions in the sequence.

See also [`n_unique_oligos`](@ref), [`nondegens`](@ref).
"""
n_deg_pos(::AbstractOligo) = 0
n_deg_pos(d::DegenOligo) = d.n_deg_pos
n_deg_pos(ov::OligoView) = count(char -> char in DEGEN_BASES, ov)
n_deg_pos(go::GappedOligo) = n_deg_pos(parent(go))

"""
    nondegens(oligo::AbstractOligo)

Return an iterator over all unique non-degenerate sequences represented by the oligomer.

For non-degenerate sequences, return a tuple containing the sequence itself.

See also [`NonDegenIterator`](@ref), [`sampleNondeg`](@ref).
"""
nondegens(oligo::Oligo) = isempty(oligo) ? Tuple{}() : (oligo,)
nondegens(go::GappedOligo) = hasgaps(go) ?
    error("Cannot iterate over sequence with gaps") :
    nondegens(DegenOligo(go))

nondegens(deg::DegenOligo) = n_deg_pos(deg) == 0 ?
    nondegens(Oligo(deg)) :
    NonDegenIterator(deg)
nondegens(ov::OligoView{T}) where T = nondegens(T(ov))

"""
    sampleChar(oligo::AbstractOligo) -> Char

Return a randomly sampled character from the oligomer sequence.

See also [`sampleView`](@ref), [`sampleNondeg`](@ref).
"""
function sampleChar(oligo::AbstractOligo)
    n = length(oligo)
    n == 0 && throw(ArgumentError("Cannot sample character from empty oligomer"))
    return oligo[rand(1:n)]
end

"""
    sampleView(oligo::AbstractOligo, len::Int) -> OligoView

Return a random contiguous subsequence (view) of the specified length.

See also [`OligoView`](@ref), [`sampleChar`](@ref).
"""
function sampleView(oligo::AbstractOligo, len::Int)
    n = length(oligo)
    len <= 0 && throw(ArgumentError("Length must be positive, got $len"))
    len > n && throw(ArgumentError("Requested view length $len exceeds oligomer length $n"))
    start = rand(1:(n - len + 1))
    return oligo[start:start+len-1]
end


_base_oligo_type(::Type{T}) where {T<:Oligo} = Oligo
_base_oligo_type(::Type{T}) where {T<:DegenOligo} = DegenOligo
_base_oligo_type(::Type{GappedOligo}) = GappedOligo
_base_oligo_type(::Type{OligoView{U}}) where {U} = _base_oligo_type(U)
_base_oligo_type(::Type{T}) where {T<:AbstractOligo} = T  # fallback

"""
    sampleNondeg(oligo::AbstractOligo) -> AbstractOligo

Return a randomly sampled non-degenerate sequence from the possible variants of the oligomer.

See also [`nondegens`](@ref), [`sample_max_gc`](@ref), [`sample_min_gc`](@ref).
"""
sampleNondeg(o::Oligo) = o
sampleNondeg(o::OligoView{Oligo}) = Oligo(o)
function sampleNondeg(d::T) where T <: AbstractOligo
    OutType = _base_oligo_type(T)
    isempty(d) && return OutType()
    buffer = Vector{Char}(undef, length(d))
    @inbounds for (i, c) in enumerate(d)
        options = get(IUPAC_B2V, c, ('-',))
        buffer[i] = rand(options)
    end
    d = description(d)
    descr = isempty(d) ? "Non-degen sample" : "Non-degen sample of $d"
    return OutType(String(buffer), descr)
end

"""
    sample_max_gc(oligo::AbstractOligo) -> AbstractOligo

Return a non-degenerate sequence sampled from the oligomer, maximizing GC content.

See also [`sample_min_gc`](@ref), [`sampleNondeg`](@ref).
"""
function sample_max_gc(d::T) where T <: AbstractOligo
    (isempty(d) || T == Oligo) && return d
    
    buffer = Vector{Char}(undef, length(d))
    @inbounds for (i, c) in enumerate(d)
        options = get(MAX_GC_OPTIONS, c, ('-',))
        buffer[i] = rand(options)
    end
    
    base_descr = description(d)
    descr = isempty(base_descr) ? "Max GC content sample" : "Max GC content sample of $base_descr"
    
    return T(String(buffer), descr)
end

"""
    sample_min_gc(oligo::AbstractOligo) -> AbstractOligo

Return a non-degenerate sequence sampled from the oligomer, minimizing GC content.

See also [`sample_max_gc`](@ref), [`sampleNondeg`](@ref).
"""
function sample_min_gc(d::T) where T <: AbstractOligo
    (isempty(d) || T == Oligo) && return d
    
    buffer = Vector{Char}(undef, length(d))
    @inbounds for (i, c) in enumerate(d)
        options = get(MIN_GC_OPTIONS, c, ('-',))
        buffer[i] = rand(options)
    end
    
    base_descr = description(d)
    descr = isempty(base_descr) ? "Min GC content sample" : "Min GC content sample of $base_descr"
    
    return T(String(buffer), descr)
end

function unfolded_proportion(args...; kwargs...)
    error(
        "`unfolded_proportion` function requires SeqFold library to be loaded.\n" *
        "In order to get this functionality, please `]add SeqFold` to your project\n" *
        "and load it with `using SeqFold`."
    )
end


include("show_oligos.jl")
end # module

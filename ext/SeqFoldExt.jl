module SeqFoldExt

using Statistics

using SeqFold

using DEPPA.Oligos
using DEPPA.Primers

"""
    SeqFold.revcomp(oligo::T) where T <: AbstractOligo -> T

Return the reverse complement of an [`AbstractOligo`](@ref).

Delegates to `SeqFold.revcomp` for the underlying string representation, using
[`Oligos.DNA_COMP_TABLE_DEG`](@ref) for fast, degenerate-aware complementation.

See also [`SeqFold.complement`](@ref), [`SeqFold.tm`](@ref).
"""
SeqFold.revcomp(oligo::T) where T <: AbstractOligo = T(SeqFold.revcomp(String(oligo); table=Oligos.DNA_COMP_TABLE_DEG), "Reverse complement of " * description(oligo))

"""
    SeqFold.complement(oligo::T) where T <: AbstractOligo -> T

Return the complement of an [`AbstractOligo`](@ref).

Delegates to `SeqFold.complement` for the underlying string representation, using
[`Oligos.DNA_COMP_TABLE_DEG`](@ref) for fast, degenerate-aware complementation.

See also [`SeqFold.revcomp`](@ref).
"""
SeqFold.complement(oligo::T) where T <: AbstractOligo = T(SeqFold.complement(String(oligo); table=Oligos.DNA_COMP_TABLE_DEG), "Complement of " * description(oligo))

"""
    SeqFold.gc_content(oligo::AbstractOligo) -> Float64

Calculate the GC content of an [`AbstractOligo`](@ref) as a fraction between 0.0 and 1.0.

For gapped oligomers, delegate to the parent sequence. For degenerate bases, use the
expected GC content from [`Oligos.IUPAC_GC_CONTENT`](@ref).

See also [`SeqFold.gc_cache`](@ref).
"""
function SeqFold.gc_content(oligo::AbstractOligo)::Float64
    hasgaps(oligo) && return SeqFold.gc_content(parent(oligo))
    isempty(oligo) && return NaN
    total_gc = sum(get(Oligos.IUPAC_GC_CONTENT, c, 0.5) for c in oligo, init=0.0)
    return total_gc / length(oligo)
end

"""
    logsumexp(x::Vector{T}) -> T

Compute the logarithm of the sum of exponentials of the elements in `x` in a numerically stable way.
"""
function logsumexp(x::Vector{T}) where T
    m = maximum(x)
    isinf(m) && return m
    s = sum(exp(val - m) for val in x, init=zero(T))
    return m + log(s)
end

"""
    SeqFold.dg(oligo::AbstractOligo; temp::Real=37.0, max_samples::Int=1000, mode::Symbol=:average) -> Float64

Compute the minimum free energy (ΔG, kcal/mol⁻¹) of an [`AbstractOligo`](@ref) at a specified temperature.

For degenerate oligomers, calculate the free energy by sampling or enumerating all unique
non-degenerate variants.

# Arguments
- `oligo::AbstractOligo`: The oligomer to fold.
- `temp::Real=37.0`: The temperature (°C) at which to perform the folding.
- `max_samples::Int=1000`: Maximum number of variants to sample if the total number exceeds this limit.
- `mode::Symbol=:average`: Calculation mode. `:average` computes the Boltzmann-weighted average ΔG,
  while `:worstcase` computes the minimum (most negative) ΔG among the sampled variants.

# Returns
- `Float64`: The calculated free energy (ΔG) in kcal/mol⁻¹, rounded to two decimal places.

See also [`SeqFold.fold`](@ref), [`SeqFold.tm`](@ref), [`Oligos.unfolded_proportion`](@ref).
"""
function SeqFold.dg(oligo::AbstractOligo; temp::Real=37.0, max_samples::Int=1000, mode::Symbol=:average)::Float64
    hasgaps(oligo) && error("Folding not supported for gapped sequences")
    isempty(oligo) && return NaN
    total_variants = n_unique_oligos(oligo)
    total_variants == 1 && return SeqFold.dg(String(oligo); temp=temp)

    T_K::Float64 = temp + 273.15
    R::Float64 = 1.9872e-3
    RT::Float64 = R * T_K
    N::Int = clamp(total_variants, 1, max_samples)
    ΔGs::Float64 = Inf64
    if mode === :average
        log_terms = Vector{Float16}(undef, N)
        if total_variants <= max_samples
            for (i, o) in enumerate(nondegens(oligo))
                ΔG = SeqFold.dg(String(o); temp=temp)
                log_terms[i] = -ΔG / RT
            end
        else
            for i in 1:N
                o = sampleNondeg(oligo)
                ΔG = SeqFold.dg(String(o); temp=temp)
                log_terms[i] = -ΔG / RT
            end
        end
        log_Q = logsumexp(log_terms) - log(N)
        ΔGs = -RT * log_Q
    elseif mode === :worstcase
        if total_variants <= max_samples
            for o in nondegens(oligo)
                ΔGs = min(SeqFold.dg(String(o); temp=temp), ΔGs)
            end
        else
            for _ in 1:max_samples
                o = sampleNondeg(oligo)
                ΔGs = min(SeqFold.dg(String(o); temp=temp), ΔGs)
            end
        end
    else
        error("Invalid mode: $mode. Supported modes are :average and :worstcase.")
    end
    return round(ΔGs, digits=2)
end

"""
    SeqFold.dg_cache(oligo::AbstractOligo; temp::Real=37.0) -> Matrix{Float64}

Compute a matrix of free energy values for all possible subsequences of an [`AbstractOligo`](@ref).

Delegates to `SeqFold.dg_cache` for the underlying non-degenerate, ungapped string representation.
Gapped sequences are not supported.

See also [`SeqFold.dg`](@ref), [`SeqFold.tm_cache`](@ref).
"""
function SeqFold.dg_cache(oligo::AbstractOligo; temp::Real=37.0)::Matrix{Float64}
    hasgaps(oligo) && error("Free energy cache not supported for gapped sequences")
    return SeqFold.dg_cache(String(Oligo(oligo)); temp=temp)
end

"""
    SeqFold.tm(oligo1::AbstractOligo, oligo2::AbstractOligo; conditions=:pcr, conf_int::Real=0.8, max_samples::Int=1000, kwargs...) -> NamedTuple

Calculate the melting temperature (Tm, °C) for the duplex formation of two [`AbstractOligo`](@ref) sequences.

For degenerate oligomers, calculate the Tm by sampling or enumerating all unique non-degenerate
variant pairs. Returns a `NamedTuple` containing the mean Tm, confidence interval, and absolute min/max bounds.

# Arguments
- `oligo1::AbstractOligo`: The first DNA sequence.
- `oligo2::AbstractOligo`: The second DNA sequence, must be the same length as `oligo1`.
- `conditions=:pcr`: Buffer conditions specification (see `SeqFold.tm`).
- `conf_int::Real=0.8`: Confidence interval for the Tm distribution (0.0, 1.0].
- `max_samples::Int=1000`: Maximum number of variant pairs to sample if the total combinations exceed this limit.
- `kwargs...`: Additional parameters to override preset conditions.

# Returns
- `NamedTuple`: `(mean=Float64, conf=Tuple{Float64, Float64}, min=Float64, max=Float64)`.

See also [`SeqFold.tm`](@ref), [`SeqFold.dg`](@ref).
"""
function SeqFold.tm(
    oligo1::AbstractOligo,
    oligo2::AbstractOligo; 
    
    conditions=:pcr,
    conf_int::Real=0.8,
    max_samples::Int=1000,
    kwargs...
)
    (hasgaps(oligo1) || hasgaps(oligo2)) && error("Melting temperature calculation not supported for gapped sequences")
    !(0 < conf_int <= 1) && throw(ArgumentError("conf_int must be in range (0, 1]"))

    i = n_unique_oligos(oligo1)
    j = n_unique_oligos(oligo2)

    if i == 1 && j == 1
        mean_Tm = SeqFold.tm(String(oligo1), String(oligo2); conditions=conditions, kwargs...)
        return (
            mean = mean_Tm,
            conf = (mean_Tm, mean_Tm),
            min = mean_Tm,
            max = mean_Tm
        )
    end

    T = zeros(Float64, min(i * j, max_samples))
    Tm_min = -Inf64
    Tm_max = Inf64
    if i * j > max_samples
        for k in 1:max_samples
            o1 = sampleNondeg(oligo1)
            o2 = sampleNondeg(oligo2)
            @inbounds T[k] = SeqFold.tm(String(o1), String(o2); conditions=conditions, kwargs...)
        end
        Tm_min_sample, Tm_max_sample = extrema(T)
        Tm_min = min(SeqFold.tm(String(sample_min_gc(oligo1)), String(sample_min_gc(oligo2)); conditions=conditions, kwargs...), Tm_min_sample)
        Tm_max = max(SeqFold.tm(String(sample_max_gc(oligo1)), String(sample_max_gc(oligo2)); conditions=conditions, kwargs...), Tm_max_sample)
    else
        counter = 1
        for o1 in nondegens(oligo1)
            for o2 in nondegens(oligo2)
                @inbounds T[counter] = SeqFold.tm(String(o1), String(o2); conditions=conditions, kwargs...)
                counter += 1
            end
        end
        Tm_min, Tm_max = extrema(T)
    end

    mean_Tm = mean(T)
    alpha = (1 - conf_int) / 2
    low = quantile(T, alpha)
    high = quantile(T, 1 - alpha)

    return (
        mean = round(mean_Tm, digits=1),
        conf = (round(low, digits=1), round(high, digits=1)),
        min = Tm_min,
        max = Tm_max
    )
end

"""
    SeqFold.tm(oligo::AbstractOligo; conditions=:pcr, conf_int::Real=0.9, max_samples::Int=1000, kwargs...) -> NamedTuple

Calculate the melting temperature (Tm, °C) for the self-complementary duplex formation of an [`AbstractOligo`](@ref).

Matches the sequence with its exact complement. See [`SeqFold.tm`](@ref) for details on arguments and return values.

See also [`SeqFold.tm`](@ref), [`SeqFold.complement`](@ref).
"""
SeqFold.tm(
    oligo::AbstractOligo;
    
    conditions=:pcr,
    conf_int::Real=0.9,
    max_samples::Int=1000,
    kwargs...
) = SeqFold.tm(oligo, SeqFold.complement(oligo); conditions=conditions, conf_int=conf_int, max_samples=max_samples, kwargs...)

"""
    SeqFold.tm_cache(oligo1::AbstractOligo, oligo2::AbstractOligo; conditions=:pcr, kwargs...) -> Matrix{Float64}

Compute a matrix of melting temperatures for all possible subsequences of an [`AbstractOligo`](@ref) pair.

Delegates to `SeqFold.tm_cache` for the underlying non-degenerate, ungapped string representation.
Gapped sequences are not supported.

See also [`SeqFold.tm_cache`](@ref), [`SeqFold.tm`](@ref).
"""
function SeqFold.tm_cache(oligo1::AbstractOligo, oligo2::AbstractOligo; conditions=:pcr, kwargs...)::Matrix{Float64}
    (hasgaps(oligo1) || hasgaps(oligo2)) && error("Melting temperature cache not supported for gapped sequences")
    SeqFold.tm_cache(String(Oligo(oligo1)), String(Oligo(oligo2)); conditions=conditions, kwargs...)
end

"""
    SeqFold.tm_cache(oligo::AbstractOligo; conditions=:pcr, kwargs...) -> Matrix{Float64}

Compute a matrix of melting temperatures for all possible subsequences of an [`AbstractOligo`](@ref) matched with its complement.

See also [`SeqFold.tm_cache`](@ref), [`SeqFold.complement`](@ref).
"""
SeqFold.tm_cache(
    oligo::AbstractOligo;
    
    conditions=:pcr,
    kwargs...
)::Matrix{Float64} = SeqFold.tm_cache(oligo, SeqFold.complement(oligo); conditions=conditions, kwargs...)

"""
    SeqFold.dot_bracket(oligo::AbstractOligo, structs::Vector{SeqFold.Structure}) -> String

Generate the dot-bracket notation representation of a predicted nucleic acid secondary structure for an [`AbstractOligo`](@ref).

Delegates to `SeqFold.dot_bracket` for the underlying string representation.

See also [`SeqFold.dot_bracket`](@ref), [`SeqFold.fold`](@ref).
"""
SeqFold.dot_bracket(
    oligo::AbstractOligo, structs::Vector{SeqFold.Structure}
) = SeqFold.dot_bracket(String(Oligo(oligo)), structs)

"""
    SeqFold.gc_cache(oligo::AbstractOligo) -> Matrix{Float64}

Compute a matrix of GC scores for all possible subsequences of an [`AbstractOligo`](@ref).

Delegates to `SeqFold.gc_cache` for the underlying non-degenerate string representation.

See also [`SeqFold.gc_cache`](@ref), [`SeqFold.gc_content`](@ref).
"""
SeqFold.gc_cache(oligo::AbstractOligo)::Matrix{Float64} = SeqFold.gc_cache(String(Oligo(oligo)))


"""
    SeqFold.tm(primer::AbstractPrimer) -> NamedTuple

Return the precomputed melting temperature (Tm) statistics stored in an [`AbstractPrimer`](@ref).

See also [`SeqFold.dg`](@ref), [`SeqFold.gc_content`](@ref).
"""
SeqFold.tm(primer::AbstractPrimer) = primer.tm

"""
    SeqFold.dg(primer::AbstractPrimer) -> Float64

Return the precomputed minimum free energy (ΔG) stored in an [`AbstractPrimer`](@ref).

See also [`SeqFold.tm`](@ref), [`SeqFold.gc_content`](@ref).
"""
SeqFold.dg(primer::AbstractPrimer) = primer.dg

"""
    SeqFold.gc_content(primer::AbstractPrimer) -> Float64

Return the precomputed GC content stored in an [`AbstractPrimer`](@ref).

See also [`SeqFold.tm`](@ref), [`SeqFold.dg`](@ref).
"""
SeqFold.gc_content(primer::AbstractPrimer) = primer.gc

"""
    Oligos.unfolded_proportion(oligo::AbstractOligo; temp::Real, max_samples::Int) -> Float64

Calculate the expected proportion of the oligomer population that remains unfolded (single-stranded) at a given temperature.

For degenerate oligomers, calculate the unfolded fraction by sampling or enumerating all unique
non-degenerate variants and averaging their individual unfolded probabilities based on their ΔG.

# Arguments
- `oligo::AbstractOligo`: The oligomer to evaluate.
- `temp::Real`: The temperature (°C) at which to evaluate the folding equilibrium.
- `max_samples::Int`: Maximum number of variants to sample if the total number exceeds this limit.

# Returns
- `Float64`: The unfolded proportion, clamped between 0.0 and 1.0.

See also [`SeqFold.dg`](@ref), [`Oligos.sampleNondeg`](@ref).
"""
function Oligos.unfolded_proportion(oligo; temp, max_samples)
    hasgaps(oligo) && error("Folding not supported for gapped sequences")
    isempty(oligo) && return NaN
    total_variants = n_unique_oligos(oligo)
    T_K::Float64 = temp + 273.15
    R = 1.9872e-3
    RT = R * T_K
    if total_variants == 1
        ΔG = SeqFold.dg(String(oligo); temp=temp)
        K_f = exp(-ΔG / RT)
        return clamp(inv(1 + K_f), 0.0, 1.0)
    end
    N = clamp(total_variants, 1, max_samples)
    unfolded_fractions = Vector{Float64}(undef, N)
    if total_variants <= max_samples
        idx = 1
        for o in nondegens(oligo)
            ΔG = SeqFold.dg(String(o); temp=temp)
            K_f = exp(-ΔG / RT)
            unfolded_fractions[idx] = inv(1 + K_f)
            idx += 1
        end
        avg_unfolded = mean(unfolded_fractions)
    else
        for i in 1:N
            o = sampleNondeg(oligo)
            ΔG = SeqFold.dg(String(o); temp=temp)
            K_f = exp(-ΔG / RT)
            unfolded_fractions[i] = 1 / (1 + K_f)
        end
        avg_unfolded = mean(unfolded_fractions)
    end
    return clamp(avg_unfolded, 0.0, 1.0)
end

"""
    Primers._ext_revcomp(o::AbstractOligo)

Compute the reverse complement of an [`AbstractOligo`](@ref) using `SeqFold`.

This is an internal extension hook used by the `Primers` module to dynamically load
thermodynamic calculations without creating a hard dependency on `SeqFold`.
"""
function Primers._ext_revcomp(o::AbstractOligo)
    SeqFold.revcomp(o)
end

"""
    Primers._ext_tm(oligo::AbstractOligo; max_samples, conf_int, conditions)

Calculate the melting temperature (Tm) of an [`AbstractOligo`](@ref) using `SeqFold`.

This is an internal extension hook used by the `Primers` module to dynamically load
thermodynamic calculations without creating a hard dependency on `SeqFold`.
"""
function Primers._ext_tm(oligo::AbstractOligo; max_samples, conf_int, conditions)
    SeqFold.tm(oligo; max_samples=max_samples, conf_int=conf_int, conditions=conditions)
end

"""
    Primers._ext_dg(oligo::AbstractOligo; max_samples, temp)

Compute the minimum free energy (ΔG) of an [`AbstractOligo`](@ref) using `SeqFold`.

This is an internal extension hook used by the `Primers` module to dynamically load
thermodynamic calculations without creating a hard dependency on `SeqFold`.
"""
function Primers._ext_dg(oligo::AbstractOligo; max_samples, temp)
    SeqFold.dg(oligo; max_samples=max_samples, temp=temp)
end

"""
    Primers._ext_gc_content(oligo::AbstractOligo)

Calculate the GC content of an [`AbstractOligo`](@ref) using `SeqFold`.

This is an internal extension hook used by the `Primers` module to dynamically load
thermodynamic calculations without creating a hard dependency on `SeqFold`.
"""
function Primers._ext_gc_content(oligo::AbstractOligo)
    SeqFold.gc_content(oligo)
end

end  # module

"""
Basic KPM kinetic calculator for reactions.

Kinetic calculator that uses KPM to predict activation
energies for reactions and uses these predictions within
the Arrhenius equation to calculate rate constants.

Basic calculator uses a flat `RT/h` Arrhenius prefactor
for all reactions. This is usually not accurate enough
to enable sensible kinetic simulations.

Implemented conditions:
* Temperature (`T`, unit: K)

Requires:
* Reaction energies (`rd.dH`, unit: eV)
* Trained KPM model (`kpm.model_path`, loaded when `KPMRun` is instantiated)

Has support for dispatching with/without a maximum rate constant
`k_max` and scaling by time unit `t_unit` (assuming rates are
provided in units of /s).
"""
mutable struct KPMBasicCalculator{kmType, uType, tType} <: KineticaCore.AbstractKineticCalculator
    Ea::Vector{uType}
    kpm::KPMRun
    k_max::kmType
    t_unit::String
    t_mult::tType
    uncertainty::Bool
end

"""
    calculator = KPMBasicCalculator(kpm[, uncertainty, k_max, t_unit])

Outer constructor method for basic KPM calculator.
"""
function KPMBasicCalculator(kpm::KPMRun;
        uncertainty::Bool=false,
        k_max::Union{Nothing, uType}=nothing, 
        t_unit::String="s") where {uType <: AbstractFloat}

    EaType = uncertainty ? Measurement : Float64
    t_mult = tconvert(t_unit, "s")
    return KPMBasicCalculator(EaType[], kpm, k_max, t_unit, t_mult, uncertainty)
end

function KineticaCore.setup_network!(sd::SpeciesData, rd::RxData, calc::KPMBasicCalculator)
    calc.Ea = calc.kpm(sd, rd)
    if eltype(calc.Ea) <: Measurement && !calc.uncertainty
        calc.Ea = Measurements.value.(calc.Ea)
    end
end

function Base.splice!(calc::KPMBasicCalculator, rids::Vector{Int})
    splice!(calc.Ea, rids)
end

"""
    rates = calculator(; T)

Calculate rates with basic KPM calculator.

Requires temperature (`T`) as a keyword argument.

Automatically dispatches to a method with correct formula for
`k_max`-aware calculation if this is defined in the underlying
`KPMBasicCalculator`.
"""
function (calc::KPMBasicCalculator{uType, uType, tType})(; T::Number) where {uType, tType}
    A = (Constants.R * T) / Constants.h
    k_r = A * exp.(-calc.Ea / (Constants.R * T)) * calc.t_mult
    return 1.0 ./ ((1.0 / calc.k_max) .+ (1.0 ./ k_r))
end

function (calc::KPMBasicCalculator{Nothing, uType, tType})(; T::Number) where {uType, tType}
    A = (Constants.R * T) / Constants.h
    k = A * exp.(-calc.Ea / (Constants.R * T)) * calc.t_mult
    return k
end

function KineticaCore.has_conditions(::KPMBasicCalculator, symbols::Vector{Symbol})
    return all([sym in [:T] for sym in symbols])
end

KineticaCore.allows_continuous(::KPMBasicCalculator) = true


"""
Collision theory-based KPM kinetic calculator for reactions.

Kinetic calculator that uses KPM to predict activation
energies for reactions and uses these predictions within
the Arrhenius equation to calculate rate constants.

Collision theory-based calculator approximates Arrhenius
prefactors on a per-reaction basis with a hard sphere
approximation of collision frequency. This is most accurate
for small, spherical reactants, and becomes less realistic
the further from this ideal the reactants get.

Additionally, allows for using one of a selection of steric
factors for correcting discrepancies between collision theory
prefactors and real Arrhenius prefactors.

Implemented conditions:
* Temperature (`T`, unit: K)

Requires:
* Reaction energies (`rd.dH`, unit: eV)
* Trained KPM model (`kpm.model_path`, loaded when `KPMRun` is instantiated)

Has support for dispatching with/without a maximum rate constant
`k_max` and scaling by time unit `t_unit` (assuming rates are
provided in units of /s).
"""
mutable struct KPMCollisionCalculator{kmType, EaType, uType, tType} <: KineticaCore.AbstractKineticCalculator
    Ea::Vector{EaType}
    μ::Vector{uType}
    σ::Vector{uType}
    ρ::Vector{uType}
    kpm::KPMRun
    k_max::kmType
    t_unit::String
    t_mult::tType
    inert_species::Union{Nothing, Vector{String}}
    steric_factor::Symbol
    steric_factor_params::Union{Nothing, uType, Vector{uType}}
    uncertainty::Bool
end

"""
    calculator = KPMCollisionCalculator(kpm[, inert_species, steric_factor, uncertainty, k_max, t_unit])

Outer constructor method for collision theory-based KPM calculator.

Collision theory requires all reactions to result from a
collision between two or more species, so unimolecular
reactions require a collision partner. If `inert_species`
has been set, the CRN will be modified to use this species as
the collision partner. Otherwise, an average collision partner
will be estimated from all species in the CRN and it will
remain unmodified.
"""
function KPMCollisionCalculator(kpm::KPMRun;
        inert_species::Union{Nothing, Vector{String}} = nothing,
        steric_factor::Symbol=:none,
        steric_factor_params::Union{Nothing, Float64, Vector{Float64}}=nothing,
        uncertainty::Bool=false,
        k_max::Union{Nothing, kmType}=nothing, 
        t_unit::String="s") where {kmType <: AbstractFloat}

    EaType = uncertainty ? Measurement : Float64
    uType = Float64 # May need a way for specifying this in future
    t_mult = tconvert(t_unit, "s")

    return KPMCollisionCalculator(EaType[], uType[], uType[], uType[],
        kpm, k_max, t_unit, t_mult, inert_species, steric_factor, 
        steric_factor_params, uncertainty)
end

function KineticaCore.setup_network!(sd::SpeciesData, rd::RxData, calc::KPMCollisionCalculator)
    if !isnothing(calc.inert_species)
        @info "Inserting inert species into unimolecular reactions."
        KineticaCore.insert_inert!(rd, sd, calc.inert_species)
    end

    Ea = calc.kpm(sd, rd)
    if eltype(Ea) <: Measurement && !calc.uncertainty
        calc.Ea = Measurements.value.(Ea)
    else
        calc.Ea = Ea
    end

    get_species_stats!(sd)
    calc.μ, calc.σ = calc_collision_params(rd, sd)
    calc.ρ = calc_steric_factors(rd, sd, calc.steric_factor; params=calc.steric_factor_params)
end

function Base.splice!(calc::KPMCollisionCalculator, rids::Vector{Int})
    splice!(calc.Ea, rids)
    splice!(calc.μ, rids)
    splice!(calc.σ, rids)
    splice!(calc.ρ, rids)
end

"""
    rates = calculator(; T)

Calculate rates with collision theory-based KPM calculator.

Requires temperature (`T`) as a keyword argument.

Automatically dispatches to a method with correct formula for
`k_max`-aware calculation if this is defined in the underlying
`KPMCollisionCalculator`.
"""
function (calc::KPMCollisionCalculator{uType, uType, tType})(; T::Number) where {uType, tType}
    k_b_dm = Constants.k_b * 1000
    
    k_r = calc.σ .* calc.ρ .* sqrt.((8*k_b_dm*T)./(pi*calc.μ)) .* exp.(-calc.Ea / (Constants.R * T)) * Constants.N_A^2 * calc.t_mult
    return 1.0 ./ ((1.0 / calc.k_max) .+ (1.0 ./ k_r))
end

function (calc::KPMCollisionCalculator{Nothing, uType, tType})(; T::Number) where {uType, tType}
    k_b_dm = Constants.k_b * 1000
    
    k = calc.σ .* calc.ρ .* sqrt.((8*k_b_dm*T)./(pi*calc.μ)) .* exp.(-calc.Ea / (Constants.R * T)) * Constants.N_A^2 * calc.t_mult
    return k
end

function KineticaCore.has_conditions(::KPMCollisionCalculator, symbols::Vector{Symbol})
    return all([sym in [:T] for sym in symbols])
end

KineticaCore.allows_continuous(::KPMCollisionCalculator) = true


"""
Collision theory-based KPM kinetic calculator for reactions, with translational entropy.

Kinetic calculator that uses KPM to predict activation
energies for reactions and uses these predictions within
the Arrhenius equation to calculate rate constants.

Collision theory-based calculator approximates Arrhenius
prefactors on a per-reaction basis with a hard sphere
approximation of collision frequency. This is most accurate
for small, spherical reactants, and becomes less realistic
the further from this ideal the reactants get.

Additionally, calculates translational entropy change for
each reaction to correct for lack of entropic contribution
in collision theory prefactors.

Implemented conditions:
* Temperature (`T`, unit: K)

Requires:
* Reaction energies (`rd.dH`, unit: eV)
* Trained KPM model (`kpm.model_path`, loaded when `KPMRun` is instantiated)

Has support for dispatching with/without a maximum rate constant
`k_max` and scaling by time unit `t_unit` (assuming rates are
provided in units of /s).
"""
mutable struct KPMCollisionEntropyCalculator{kmType, EaType, uType, tType} <: KineticaCore.AbstractKineticCalculator
    Ea::Vector{EaType}
    μ::Vector{uType}
    σ::Vector{uType}
    ΔS_t::Vector{uType}
    ΔN_m::Vector{Int}
    kpm::KPMRun
    k_max::kmType
    t_unit::String
    t_mult::tType
    inert_species::Union{Nothing, Vector{String}}
    uncertainty::Bool
end

"""
    calculator = KPMCollisionEntropyCalculator(kpm[, inert_species, steric_factor, uncertainty, k_max, t_unit])

Outer constructor method for collision theory-based KPM calculator with translational entropy.

Collision theory requires all reactions to result from a
collision between two or more species, so unimolecular
reactions require a collision partner. If `inert_species`
has been set, the CRN will be modified to use this species as
the collision partner. Otherwise, an average collision partner
will be estimated from all species in the CRN and it will
remain unmodified.
"""
function KPMCollisionEntropyCalculator(kpm::KPMRun;
        inert_species::Union{Nothing, Vector{String}} = nothing,
        uncertainty::Bool=false,
        k_max::Union{Nothing, kmType}=nothing, 
        t_unit::String="s") where {kmType <: AbstractFloat}

    EaType = uncertainty ? Measurement : Float64
    uType = Float64 # May need a way for specifying this in future
    t_mult = tconvert(t_unit, "s")

    return KPMCollisionEntropyCalculator(EaType[], uType[], uType[], uType[], Int[],
        kpm, k_max, t_unit, t_mult, inert_species, uncertainty)
end

function KineticaCore.setup_network!(sd::SpeciesData, rd::RxData, calc::KPMCollisionEntropyCalculator)
    if !isnothing(calc.inert_species)
        @info "Inserting inert species into unimolecular reactions."
        KineticaCore.insert_inert!(rd, sd, calc.inert_species)
    end

    Ea = calc.kpm(sd, rd)
    if eltype(Ea) <: Measurement && !calc.uncertainty
        calc.Ea = Measurements.value.(Ea)
    else
        calc.Ea = Ea
    end

    get_species_stats!(sd)
    calc.μ, calc.σ = calc_collision_params(rd, sd)
    calc.ΔS_t, calc.ΔN_m = calc_entropy_change(rd, sd)
end

function Base.splice!(calc::KPMCollisionEntropyCalculator, rids::Vector{Int})
    splice!(calc.Ea, rids)
    splice!(calc.μ, rids)
    splice!(calc.σ, rids)
    splice!(calc.ΔS_t, rids)
    splice!(calc.ΔN_m, rids)
end

"""
    rates = calculator(; T)

Calculate rates with collision theory-based KPM calculator, with translational entropy.

Requires temperature (`T`) as a keyword argument.

Automatically dispatches to a method with correct formula for
`k_max`-aware calculation if this is defined in the underlying
`KPMCollisionCalculator`.
"""
function (calc::KPMCollisionEntropyCalculator{uType, uType, tType})(; T::Number) where {uType, tType}
    k_b_dm = Constants.k_b * 1000

    entropic_temperature = 3*Constants.R*log(ℯ, T)/2
    ΔS_t = calc.ΔS_t .+ (calc.ΔN_m * entropic_temperature)
    
    k_r = calc.σ .* sqrt.((8*k_b_dm*T)./(pi*calc.μ)) .* exp.(-calc.Ea / (Constants.R * T)) .* exp.(ΔS_t / Constants.R) * Constants.N_A^2 * calc.t_mult
    return 1.0 ./ ((1.0 / calc.k_max) .+ (1.0 ./ k_r))
end

function (calc::KPMCollisionEntropyCalculator{Nothing, uType, tType})(; T::Number) where {uType, tType}
    k_b_dm = Constants.k_b * 1000
    
    entropic_temperature = 3*Constants.R*log(ℯ, T)/2
    ΔS_t = calc.ΔS_t .+ (calc.ΔN_m * entropic_temperature)
    
    k = calc.σ .* sqrt.((8*k_b_dm*T)./(pi*calc.μ)) .* exp.(-calc.Ea / (Constants.R * T)) .* exp.(ΔS_t / Constants.R) * Constants.N_A^2 * calc.t_mult
    return k
end

function KineticaCore.has_conditions(::KPMCollisionEntropyCalculator, symbols::Vector{Symbol})
    return all([sym in [:T] for sym in symbols])
end

KineticaCore.allows_continuous(::KPMCollisionEntropyCalculator) = false
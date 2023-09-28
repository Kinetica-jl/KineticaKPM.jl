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
struct KPMBasicCalculator{kmType, uType, tType} <: KineticaCore.AbstractKineticCalculator
    Ea::Vector{uType}
    k_max::kmType
    t_unit::String
    t_mult::tType
end

"""
    calculator = KPMBasicCalculator(rd, sd, kpm[, k_max, t_unit])

Outer constructor method for basic KPM calculator.
"""
function KPMBasicCalculator(rd::RxData, sd::SpeciesData, kpm::KPMRun;
        uncertainty::Bool=false,
        k_max::Union{Nothing, uType}=nothing, 
        t_unit::String="s") where {uType <: AbstractFloat}

    Ea = kpm(rd, sd)
    if !uncertainty
        Ea = Measurements.value.(Ea)
    end

    t_mult = tconvert(t_unit, "s")
    return KPMBasicCalculator(Ea, k_max, t_unit, t_mult)
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
struct KPMCollisionCalculator{kmType, EaType, uType, tType} <: KineticaCore.AbstractKineticCalculator
    Ea::Vector{EaType}
    μ::Vector{uType}
    σ::Vector{uType}
    ρ::Vector{uType}
    k_max::kmType
    t_unit::String
    t_mult::tType
    steric_factor::Symbol
end

"""
    calculator = KPMCollisionCalculator(rd, sd, kpm[, steric_type, k_max, t_unit])

Outer constructor method for collision theory-based KPM calculator.

Collision theory requires all reactions to result from a
collision between two or more species, so unimolecular
reactions require a collision partner. If `pars.inert_species`
has been set, the CRN will be modified to use this species as
the collision partner. Otherwise, an average collision partner
will be estimated from all species in the CRN and it will
remain unmodified.
"""
function KPMCollisionCalculator(rd::RxData, sd::SpeciesData,
        pars::ODESimulationParams, kpm::KPMRun;
        steric_factor::Symbol=:none,
        uncertainty::Bool=false,
        k_max::Union{Nothing, uType}=nothing, 
        t_unit::String="s") where {uType <: AbstractFloat}

    if !isnothing(pars.inert_species)
        @info "Inserting inert species into unimolecular reactions."
        KineticaCore.insert_inert!(rd, sd, pars.inert_species)
    end

    Ea = kpm(rd, sd)
    if !uncertainty
        Ea = Measurements.value.(Ea)
    end

    get_species_stats!(sd)
    μ, σ = calc_collision_params(rd, sd)
    ρ = calc_steric_factors(rd, sd, steric_factor)
    t_mult = tconvert(t_unit, "s")

    return KPMCollisionCalculator(Ea, μ, σ, ρ, k_max, t_unit, t_mult, steric_factor)
end

"""
    rates = calculator(; T)

Calculate rates with basic KPM calculator.

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
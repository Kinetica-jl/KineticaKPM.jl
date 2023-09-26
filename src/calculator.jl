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
        k_max::Union{Nothing, uType}=nothing, 
        t_unit::String="s") where {uType <: AbstractFloat}
    Ea = kpm(rd, sd)
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
    R = 8.314462618 # Gas constant (J/K/mol)
    h = 6.626070e-34 # Planck constant (Js)

    A = (R * T) / h
    k_r = A * exp.(-calc.Ea / (R * T)) * calc.t_mult
    return 1.0 ./ ((1.0 / calc.k_max) .+ (1.0 ./ k_r))
end

function (calc::KPMBasicCalculator{Nothing, uType, tType})(; T::Number) where {uType, tType}
    R = 8.314462618 # Gas constant (J/K/mol)
    h = 6.626070e-34 # Planck constant (Js)

    A = (R * T) / h
    k = A * exp.(-calc.Ea / (R * T)) * calc.t_mult
    return k
end

function has_conditions(::KPMBasicCalculator, symbols::Vector{Symbol})
    return all([sym in [:T] for sym in symbols])
end

allows_continuous(::KPMBasicCalculator) = true
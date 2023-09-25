"""
"""
struct KPMBasicCalculator{kmType, uType, tType} <: KineticaCore.AbstractKineticCalculator
    Ea::Vector{uType}
    k_max::kmType
    t_unit::String
    t_mult::tType
end

"""
"""
function KPMBasicCalculator(rd::RxData, sd::SpeciesData, kpm::KPMRun;
        k_max::Union{Nothing, uType}=nothing, 
        t_unit::String="s") where {uType <: AbstractFloat}
    Ea = kpm(rd, sd)
    t_mult = tconvert(t_unit, "s")
    return KPMBasicCalculator(Ea, k_max, t_unit, t_mult)
end

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
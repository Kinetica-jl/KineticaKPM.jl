"""
    μ, σ = calc_collision_params(rd, sd[, easy_units])

Calculates the collision theory parameters for given reactions.

Collision theory requires calculation of a reduced mass μ and a
collision cross-section σ for each pair of reactants in a reaction.

Calculates these values, using an 'average collision partner' for
unimolecular reactions. Possible to get around this by implementing
a non-reactive collision partner with unit concentration to all 
unimolecular reactions, effectively making them bimolecular for the
purposes of these calculations.

Must be preceded by a call to `get_frag_stats!()` to populate species
weights and radii.
"""
function calc_collision_params(rd::RxData, sd::SpeciesData; easy_units::Bool=true)
    m_avg = mean(values(sd.cache[:weights]))
    r_avg = mean(values(sd.cache[:radii]))
    μ_AB = zeros(Float64, rd.nr)
    d_AB = zeros(Float64, rd.nr)
    for i in 1:rd.nr
        if length(rd.reacs[i]) == 1
            m_A = sd.cache[:weights][rd.id_reacs[i][1]]
            r_A = sd.cache[:radii][rd.id_reacs[i][1]]
            if rd.stoic_reacs[i][1] == 1
                m_B = m_avg
                r_B = r_avg
            else
                m_B = m_A
                r_B = r_A
            end
            μ_AB[i] = (m_A * m_B)/(m_A + m_B)
            d_AB[i] = r_A + r_B
        else
            m_reacs = [sd.cache[:weights][rd.id_reacs[i][j]] for j in 1:length(rd.reacs[i])]
            μ_AB[i] = prod(m_reacs)/sum(m_reacs)
            d_AB[i] = sum([sd.cache[:radii][rd.id_reacs[i][j]] for j in 1:length(rd.reacs[i])])
        end
    end

    if easy_units
        μ_kg = μ_AB * Constants.amu_to_kg
        d_dm = d_AB * Constants.Ang_to_dm
        σ_dm = pi * (d_dm.^2)
        return μ_kg, σ_dm
    else
        σ_AB = pi * (d_AB.^2)
        return μ_AB, σ_AB
    end
end


"""
    ρ = calc_steric_factors(rd, sd, steric_factor)

Calculates steric factors for all reactions in `rd`, using the requested steric factor function.

Valid steric factors are
* `:basic` - Calculates steric factors as `1/(α_A * α_B)`, where `α_i = n_i² + 5r_i(n_i - 1)`.
* `:exp` - Calculates steric factors as `1/(α_A * α_B)^β`, where `α_i = n_i² + 5r_i(n_i - 1)`. Requires passing a value of β via `params`
* `:logistic` - Calculates steric factors with an adjustable bivariate logistic distribution. Requires passing a value of β via `params`
* `:dlogistic` - Calculates steric factors with 2 adjustable bivariate logistic distributions. Requires passing a vector of [β_assoc, β_dissoc] via `params`
* `:none` - Makes all steric factors equal to 1, removing them from the rate equation.
"""
function calc_steric_factors(rd::RxData, sd::SpeciesData, steric_factor::Symbol; params=nothing)
    if !(steric_factor in [:basic, :exp, :logistic, :dlogistic, :none])
        throw(ArgumentError("Unknown steric factor type: $steric_factor"))
    end

    if steric_factor in [:exp, :logistic]
        if isnothing(params)
            error("Steric factor $steric_factor requires a value for beta, which must be passed through steric_factor_params.")
        end
        ρ = calc_steric_factors(rd, sd, Val(steric_factor), params)
    elseif steric_factor == :dlogistic
        if isnothing(params)
            error("Steric factor $steric_factor requires values for associative and dissociative beta, which must be passed through steric_factor_params.")
        end
        ρ = calc_steric_factors(rd, sd, Val(steric_factor), params[1], params[2])
    else
        ρ = calc_steric_factors(rd, sd, Val(steric_factor))
    end
    return ρ
end

"""
    ρ = calc_steric_factors(rd, sd, Val(:none))

Makes all steric factors equal to 1, removing them from the rate equation.
"""
function calc_steric_factors(rd::RxData, sd::SpeciesData, ::Val{:none})
    ρ = ones(Float64, rd.nr)
    return ρ
end

"""
    ρ = calc_steric_factors(rd, uniq_frags, Val(:basic))

Calculates steric factors as `1/(α_A * α_B)`, where `α_i = n_i + 5r_i(n_i - 1)`.
"""
function calc_steric_factors(rd::RxData, sd::SpeciesData, ::Val{:basic})
    n_avg = mean([sd.xyz[i]["N_atoms"] for i in 1:sd.n])
    r_avg = mean(values(sd.cache[:radii]))
    ρ = zeros(Float64, rd.nr)

    for i in 1:rd.nr
        n_A = sd.xyz[rd.id_reacs[i][1]]["N_atoms"]
        r_A = sd.cache[:radii][rd.id_reacs[i][1]]
        if length(rd.reacs[i]) == 1
            if rd.stoic_reacs[i][1] == 1
                n_B = n_avg
                r_B = r_avg
            else
                n_B = n_A
                r_B = r_A
            end
        else
            n_B = sd.xyz[rd.id_reacs[i][2]]["N_atoms"]
            r_B = sd.cache[:radii][rd.id_reacs[i][2]]
        end
        α_A = n_A^2 + (n_A - 1) * 5 * r_A
        α_B = n_B^2 + (n_B - 1) * 5 * r_B
        α = (α_A * α_B) / 100

        ρ[i] = 1 / (100 * α)
    end

    return ρ
end

"""
    ρ = calc_steric_factors(rd, uniq_frags, Val(:exp), β)

Calculates steric factors as `1/(α_A * α_B)^β`, where `α_i = n_i² + 5r_i(n_i - 1)`.
"""
function calc_steric_factors(rd::RxData, sd::SpeciesData, ::Val{:exp}, β::Float64)
    n_avg = mean([sd.xyz[i]["N_atoms"] for i in 1:sd.n])
    r_avg = mean(values(sd.cache[:radii]))
    ρ = zeros(Float64, rd.nr)

    for i in 1:rd.nr
        n_A = sd.xyz[rd.id_reacs[i][1]]["N_atoms"]
        r_A = sd.cache[:radii][rd.id_reacs[i][1]]
        if length(rd.reacs[i]) == 1
            if rd.stoic_reacs[i][1] == 1
                n_B = n_avg
                r_B = r_avg
            else
                n_B = n_A
                r_B = r_A
            end
        else
            n_B = sd.xyz[rd.id_reacs[i][2]]["N_atoms"]
            r_B = sd.cache[:radii][rd.id_reacs[i][2]]
        end
        α_A = n_A^2 + (n_A - 1) * 5 * r_A
        α_B = n_B^2 + (n_B - 1) * 5 * r_B
        α = (α_A * α_B)^β

        ρ[i] = 1 / α
    end

    return ρ
end

"""
    ρ = calc_steric_factors(rd, uniq_frags, Val(:logistic), β)

Calculates steric factors on an adjustable bivariate logistic distribution.

The steric factor ρ for a reaction is calculated as

    ρ = D / ((1+exp(βα_A))(1+exp(βα_B)))

where 

    D = (1+exp(β))^2, 
    α_i = n_i + 5r_i(n_i-1), 
    
and A and B are the bimolecular reactants.
"""
function calc_steric_factors(rd::RxData, sd::SpeciesData, ::Val{:logistic}, β::Float64)
    n_avg = mean([sd.xyz[i]["N_atoms"] for i in 1:sd.n])
    r_avg = mean(values(sd.cache[:radii]))
    ρ = zeros(Float64, rd.nr)
    D = (1+exp(β))^2

    for i in 1:rd.nr
        n_A = sd.xyz[rd.id_reacs[i][1]]["N_atoms"]
        r_A = sd.cache[:radii][rd.id_reacs[i][1]]
        if length(rd.reacs[i]) == 1
            if rd.stoic_reacs[i][1] == 1
                n_B = n_avg
                r_B = r_avg
            else
                n_B = n_A
                r_B = r_A
            end
        else
            n_B = sd.xyz[rd.id_reacs[i][2]]["N_atoms"]
            r_B = sd.cache[:radii][rd.id_reacs[i][2]]
        end
        α_A = n_A + (n_A - 1) * 5 * r_A
        α_B = n_B + (n_B - 1) * 5 * r_B

        ρ[i] = D / ((1+exp(β*α_A))*(1+exp(β*α_B)))
    end

    return ρ
end

"""
    ρ = calc_steric_factors(rd, uniq_frags, Val(:dlogistic), β_assoc, β_dissoc)

Calculates steric factors on 2 adjustable bivariate logistic distributions.

Creates 2 logistic distributions: one for associative reactions, and the other
for dissociative reactions. The gradient of each distribution can be tuned with
its respective β parameter.

The steric factor ρ for an associative reaction is calculated as

    ρ = D_assoc / ((1+exp(β_assoc*α_A))(1+exp(β_assoc*α_B)))

while the steric factor for a dissociative reaction is calculated as
    
    ρ = D_dissoc / ((1+exp(β_dissoc*α_A))(1+exp(β_dissoc*α_B)))

where 

    D_assoc = (1+exp(β_assoc))^2, 
    D_dissoc = (1+exp(β_dissoc))^2,
    α_i = n_i + 5r_i(n_i-1), 
    
and A and B are the bimolecular reactants. Reactions which are neither
associative nor dissociative use the associative steric factor by default.
"""
function calc_steric_factors(rd::RxData, sd::SpeciesData, ::Val{:dlogistic}, β_assoc::Float64, β_dissoc::Float64)
    n_avg = mean([sd.xyz[i]["N_atoms"] for i in 1:sd.n])
    r_avg = mean(values(sd.cache[:radii]))
    ρ = zeros(Float64, rd.nr)
    D_assoc = (1+exp(β_assoc))^2
    D_dissoc = (1+exp(β_dissoc))^2

    for i in 1:rd.nr
        ΔN_m = sum(rd.stoic_prods[i]) - sum(rd.stoic_reacs[i])

        n_A = sd.xyz[rd.id_reacs[i][1]]["N_atoms"]
        r_A = sd.cache[:radii][rd.id_reacs[i][1]]
        if length(rd.reacs[i]) == 1
            if rd.stoic_reacs[i][1] == 1
                n_B = n_avg
                r_B = r_avg
            else
                n_B = n_A
                r_B = r_A
            end
        else
            n_B = sd.xyz[rd.id_reacs[i][2]]["N_atoms"]
            r_B = sd.cache[:radii][rd.id_reacs[i][2]]
        end
        α_A = n_A + (n_A - 1) * 5 * r_A
        α_B = n_B + (n_B - 1) * 5 * r_B

        if ΔN_m > 0 # n_reactants < n_products therefore dissociative
            ρ[i] = D_dissoc / ((1+exp(β_dissoc*α_A))*(1+exp(β_dissoc*α_B)))
        else
            ρ[i] = D_assoc / ((1+exp(β_assoc*α_A))*(1+exp(β_assoc*α_B)))
        end
    end

    return ρ
end


"""
"""
function calc_entropy_change(rd::RxData, sd::SpeciesData)
    ΔN_m = [sum(rd.stoic_prods[i]) - sum(rd.stoic_reacs[i]) for i in 1:rd.nr]

    ΔS_t = zeros(rd.nr)
    for i in 1:rd.nr
        mass_term = prod([sd.cache[:weights][spec] for spec in rd.id_reacs[i]]) / prod([sd.cache[:weights][spec] for spec in rd.id_prods[i]])
        ΔS_t[i] = ΔN_m[i] * (3*Constants.R*log(ℯ, 2*pi*Constants.k_b/(Constants.h^2))/2 +
                  Constants.R*log(ℯ, 1000.0) + 5*Constants.R/2) + 3*Constants.R*mass_term/2
    end

    return ΔS_t, ΔN_m
end
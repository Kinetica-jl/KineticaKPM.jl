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
* `:basic` - Calculates steric factors as `1/(α_A * α_B)`, where `α_i = n_i + 5r_i(n_i - 1)`.
* `:none` - Makes all steric factors equal to 1, removing them from the rate equation.
"""
function calc_steric_factors(rd::RxData, sd::SpeciesData, steric_factor::Symbol)
    if !(steric_factor in [:basic, :none])
        throw(ArgumentError("Unknown steric factor type: $steric_factor"))
    end
        ρ = calc_steric_factors(rd, sd, Val(steric_factor))
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
        α_A = n_A + (n_A - 1) * 5 * r_A
        α_B = n_B + (n_B - 1) * 5 * r_B
        α = (α_A * α_B) / 100

        ρ[i] = 1 / (100 * α)
    end

    return ρ
end
"""
"""
struct KPMArgs
    model::String
    reactant::String
    product::String
    enthalpy::String
    outfile::String
    direction::String
    uncertainty::String
    fix_radicals::String
    suppress_rdlogs::String
    verbose::String
end

"""
"""
mutable struct KPMRun
    model_path::String
    predictor::PyObject

    function KPMRun(model_path)
        args = KPMArgs(
            model_path,
            "reacs.xyz",
            "prods.xyz",
            "dH.txt",
            "preds.txt",
            "both",
            "True",
            "True",
            "True",
            "True"
        )
        _predict = pyimport("KPM.predict")
        pysys.stdout = typeof(current_logger().stream) == IOStream ? current_logger().stream : pysys.__stdout__
        predictor = _predict.ModelPredictor(args)
        flush_log()
        pysys.stdout = pysys.__stdout__

        return new(model_path, predictor)
    end
end

"""
"""
function (self::KPMRun)(rd::RxData, sd::SpeciesData; rdir::Union{String, Nothing}=nothing)
    calc_dir = isnothing(rdir) ? mktempdir(; prefix="kinetica_kpm_") : rdir
    @info "Predicting Ea of all reactions in $(calc_dir)"
    
    @info "Optimising all reactant/product systems."
    flush_log()
    reac_molsys = []
    prod_molsys = []
    for rcount in 1:rd.nr
        reac_xyzs = [sd.xyz[idx] for (i, idx) in enumerate(rd.id_reacs[rcount]) for _ in 1:rd.stoic_reacs[rcount][i]]
        push!(reac_molsys, system_from_mols(reac_xyzs))
    
        prod_xyzs = [sd.xyz[idx] for (i, idx) in enumerate(rd.id_prods[rcount]) for _ in 1:rd.stoic_prods[rcount][i]]
        push!(prod_molsys, system_from_mols(prod_xyzs))
    end
    @info "Reactant/product systems optimised."

    write_frames(joinpath(calc_dir, "reacs.xyz"), reac_molsys)
    write_frames(joinpath(calc_dir, "prods.xyz"), prod_molsys)

    dH_conv = rd.dH .* Constants.eV_to_kcal_per_mol
    writedlm(joinpath(calc_dir, "dH.txt"), vcat(dH_conv...))

    @info "Predicting Ea..."
    outfile = open(joinpath(calc_dir, "kpm.out"), "w")
    pysys.stdout = PyTextIO(outfile)
    curr_dir = pwd()
    cd(calc_dir)
    diffs = self.predictor.process_xyzs()
    self.predictor.predict(diffs)
    cd(curr_dir)
    pysys.stdout = pysys.__stdout__
    close(outfile)
    @info "Prediction complete.\n"

    Ea = read_predictions(joinpath(calc_dir, "preds.txt"))

    return Ea    
end

"""
"""
function (self::KPMRun)(rd::RxData, sd::SpeciesData, rcount::Int; rdir::Union{String, Nothing}=nothing)
    calc_dir = isnothing(rdir) ? mktempdir(; prefix="kinetica_kpm_") : rdir
    @info "Predicting Ea of reaction $rcount in $(calc_dir)"

    reac_xyzs = [sd.xyz[idx] for (i, idx) in enumerate(rd.id_reacs[rcount]) for _ in 1:rd.stoic_reacs[rcount][i]]
    system_from_mols(reac_xyzs, joinpath(calc_dir, "reacs.xyz"))
    @info "Optimised reactant system."

    prod_xyzs = [sd.xyz[idx] for (i, idx) in enumerate(rd.id_prods[rcount]) for _ in 1:rd.stoic_prods[rcount][i]]
    system_from_mols(prod_xyzs, joinpath(calc_dir, "prods.xyz"))
    @info "Optimised product system."

    dH = rd.dH[rcount] * Constants.eV_to_kcal_per_mol
    open(joinpath(calc_dir, "dH.txt"), "w") do f
        write(f, dH)
    end

    @info "Predicting Ea..."
    outfile = open(joinpath(calc_dir, "kpm.out"), "w")
    pysys.stdout = PyTextIO(outfile)
    curr_dir = pwd()
    cd(calc_dir)
    diffs = self.predictor.process_xyzs()
    self.predictor.predict(diffs)
    cd(curr_dir)
    pysys.stdout = pysys.__stdout__
    close(outfile)
    @info "Prediction complete.\n"
    
    Ea = read_predictions(joinpath(calc_dir, "preds.txt"))

    return Ea
end


"""
    Ea = read_predictions(predfile)

Reads in KPM predictions from a preds.txt
"""
function read_predictions(predfile::String)
    # Read in activation energies as a vector of strings.
    Eact_strs = vec(readdlm(predfile, '\n', String, '\n';
        comments=true, comment_char='#', skipblanks=true))
    
    # Parse with/without uncertainties attached.
    if occursin('±', Eact_strs[1])
        # Remove the uncertainty component from the strings if present.
        splits = split.(Eact_strs, " ")
        Eact_strs = getindex.(splits, 1)
        uncert_strs = getindex.(splits, 3)
        Eacts = parse.(Float64, Eact_strs)
        Eact_uncerts = parse.(Float64, uncert_strs)
        Ea = measurement.(Eacts, Eact_uncerts)
    else
        Ea = parse.(Float64, Eact_strs)
    end

    return Ea
end
"""
Data container for KPM initialisation arguments.

Not intended to be user-accessible, instead called internally
when constructing KPMRun.
"""
struct KPMArgs
    model::String
    reactant::String
    product::String
    enthalpy::String
    outfile::Union{String, Nothing}
    direction::String
    uncertainty::String
    fix_radicals::String
    suppress_rdlogs::String
    verbose::String
end

"""
    kpm = KPMRun(model_path)

KPM runner for prediction of reaction activation energies.

Handles instantiation of the underlying Scikit-Learn model
in Python, and can be called to use this model for predictions.
"""
mutable struct KPMRun
    model_path::String
    predictor::Py

    function KPMRun(model_path)
        args = KPMArgs(
            model_path,
            "reacs.xyz",
            "prods.xyz",
            "dH.txt",
            nothing,
            "forward",
            "True",
            "True",
            "True",
            "True"
        )
        _predict = pyimport("KPM.predict")
        pysys.stdout = typeof(current_logger().stream) == IOStream ? pytextio(current_logger().stream) : pysys.__stdout__
        predictor = _predict.ModelPredictor(args)
        flush_log()
        pysys.stdout = pysys.__stdout__

        return new(model_path, predictor)
    end
end

"""
    Ea = kpm(rd[, rcount, rdir])

Predicts activation energies using KPM.

Assembles joines SMILES of arectants and products, creates
reaction difference fingerprints and performs predictions.

If `rcount` is not provided, predicts Ea for all reactions
currently in `rd`. If it is provided, only predicts Ea for
the specified reaction.

If `rdir` is provided, saves all KPM inputs/outputs to this
directory (which must exist before calling). Otherwise, does
all IO within a temporary directory that is deleted once
calculations are finished.
"""
function (self::KPMRun)(sd::SpeciesData, rd::RxData; rdir::Union{String, Nothing}=nothing)
    calc_dir = isnothing(rdir) ? mktempdir(; prefix="kinetica_kpm_") : rdir
    @info " - Predicting Ea of all reactions in $(calc_dir)"
    flush_log()

    outfile = open(joinpath(calc_dir, "kpm.out"), "w")
    pysys.stdout = pytextio(outfile)

    rsmi = [join(sort(reduce(vcat, [[sd.toStr[rids[i]] for _ in 1:rstoics[i]] for i in axes(rids, 1)])), ".") for (rids, rstoics) in zip(rd.id_reacs, rd.stoic_reacs)]
    self.predictor.rsmi = pylist(rsmi)
    rmol = [rdChem.MolFromSmiles(smi) for smi in rsmi]
    psmi = [join(sort(reduce(vcat, [[sd.toStr[pids[i]] for _ in 1:pstoics[i]] for i in axes(pids, 1)])), ".") for (pids, pstoics) in zip(rd.id_prods, rd.stoic_prods)]
    self.predictor.psmi = pylist(psmi)
    pmol = [rdChem.MolFromSmiles(smi) for smi in psmi]
    self.predictor.dH_arr = pylist(rd.dH .* Constants.eV_to_kcal_per_mol)
    self.predictor.num_reacs = rd.nr

    diffs = kpm_utils.descriptors.calc_diffs(rd.nr, self.predictor.descriptor_type, rmol, pmol, 
        self.predictor.dH_arr, self.predictor.morgan_radius, self.predictor.morgan_num_bits)
    Ea, uncerts = self.predictor.predict(diffs)

    pysys.stdout = pysys.__stdout__
    close(outfile)
    @info " - Prediction complete.\n"
    flush_log()

    Ea = pyconvert(Vector, Ea) * Constants.kcal_to_J
    uncerts = pyconvert(Vector, uncerts) * Constants.kcal_to_J
    Ea_final = measurement.(Ea, uncerts)

    return Ea_final
end

function (self::KPMRun)(sd::SpeciesData, rd::RxData, rcount::Int; rdir::Union{String, Nothing}=nothing)
    calc_dir = isnothing(rdir) ? mktempdir(; prefix="kinetica_kpm_") : rdir
    @info "Predicting Ea of reaction $rcount in $(calc_dir)"
    flush_log()

    outfile = open(joinpath(calc_dir, "kpm.out"), "w")
    pysys.stdout = outfile

    rsmi = [join(sort(reduce(vcat, [[sd.toStr[rd.id_reacs[rcount][i]] for _ in 1:rd.stoic_reacs[rcount][i]] for i in axes(rd.id_reacs[rcount], 1)])), ".")]
    self.predictor.rsmi = pylist(rsmi)
    rmol = [rdChem.MolFromSmiles(smi) for smi in rsmi]
    psmi = [join(sort(reduce(vcat, [[sd.toStr[rd.id_prods[rcount][i]] for _ in 1:rd.stoic_prods[rcount][i]] for i in axes(rd.id_prods[rcount], 1)])), ".")]
    self.predictor.psmi = pylist(psmi)
    pmol = [rdChem.MolFromSmiles(smi) for smi in psmi]
    self.predictor.dH_arr = pylist([rd.dH[rcount] .* Constants.eV_to_kcal_per_mol])
    self.predictor.num_reacs = 1

    diff = kpm_utils.descriptors.calc_diffs(1, self.predictor.descriptor_type, rmol, pmol, 
        self.predictor.dH_arr, self.predictor.morgan_radius, self.predictor.morgan_num_bits)
    Ea, uncerts = self.predictor.predict(diff)

    pysys.stdout = pysys.__stdout__
    close(outfile)
    @info "Prediction complete.\n"
    flush_log()
    
    Ea = pyconvert(Vector, Ea) * Constants.kcal_to_J
    uncerts = pyconvert(Vector, uncerts) * Constants.kcal_to_J
    Ea_final = measurement(Ea[1], uncerts[1])

    return Ea_final
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
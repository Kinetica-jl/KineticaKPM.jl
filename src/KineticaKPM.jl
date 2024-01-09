module KineticaKPM

using Kinetica
Constants = Kinetica.Constants

using Logging
using PyCall
using ExtXYZ
using Measurements
using DelimitedFiles
using Statistics

const version = VersionNumber(0, 3, 0)
const rdChem = PyNULL()
const kpm_utils = PyNULL()
function __init__()
    copy!(rdChem, pyimport("rdkit.Chem"))
    copy!(kpm_utils, pyimport("KPM.utils"))
end

include("runner.jl")
export KPMRun
include("collision_utils.jl")
include("calculator.jl")
export KPMBasicCalculator, KPMCollisionCalculator, KPMCollisionEntropyCalculator

end

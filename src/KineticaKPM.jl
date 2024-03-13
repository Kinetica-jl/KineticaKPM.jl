"""
KineticaKPM.jl

UK Ministry of Defence © Crown Owned Copyright 2024/AWE​
"""
module KineticaKPM

using Kinetica
Constants = Kinetica.Constants

using Logging
using PythonCall
using ExtXYZ
using Measurements
using DelimitedFiles
using Statistics

const version = VersionNumber(0, 3, 1)
const rdChem = PythonCall.pynew()
const kpm_utils = PythonCall.pynew()
const pysys = PythonCall.pynew()
function __init__()
    PythonCall.pycopy!(rdChem, pyimport("rdkit.Chem"))
    PythonCall.pycopy!(kpm_utils, pyimport("KPM.utils"))
    PythonCall.pycopy!(pysys, pyimport("sys"))
end

include("runner.jl")
export KPMRun
include("collision_utils.jl")
include("calculator.jl")
export KPMBasicCalculator, KPMCollisionCalculator, KPMCollisionEntropyCalculator

end

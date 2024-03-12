using Conda

function flush_log(logger=Base.global_logger())
    flush(logger.stream)
end

PROJECT_DIR = dirname(dirname(@__FILE__))

# Check git is installed.
gitcheck = run(Cmd(`which git`, ignorestatus=true))
if gitcheck.exitcode == 1
    error("Git not found on the system. Install Git before running this script.")
end

# Get dependencies
depsdir = joinpath(PROJECT_DIR, "submodules")
kpmdir = joinpath(depsdir, "KineticPredictorModel")
if readdir(kpmdir) != String[]
    @info string("Submodules already initialised, checking for updates...")
    gitcmd = run(Cmd(`git submodule update`, dir=PROJECT_DIR))
    if gitcmd.exitcode == 1
        error("Error updating Git submodules.")
    end
    @info string("Submodules up to date.")
else
    @info string("Obtaining and initialising submodules...")
    gitcmd = run(Cmd(`git submodule init`, dir=PROJECT_DIR))
    if gitcmd.exitcode == 1
        error("Error initialising Git submodules.")
    end
    gitcmd = run(Cmd(`git submodule update --init`, dir=PROJECT_DIR))
    if gitcmd.exitcode == 1
        error("Error initialising Git submodules.")
    end
    @info string("All submodules initialised.")
end
flush_log()

# Ensure the Conda environment has the neccessary dependencies.
ignore_python_deps = lowercase(get(ENV, "KINETICA_BUILD_IGNORE_CONDA", "TRUE"))
if ignore_python_deps == "true"
    @info "KINETICA_BUILD_IGNORE_CONDA is TRUE, ignoring Python dependencies."
else
    @info "Setting up Python dependencies..."
    flush_log()
    kpmdeps = readlines(joinpath(kpmdir, "requirements.txt"))
    Conda.add(kpmdeps, channel="conda-forge")
    Conda.pip_interop(true)
    Conda.pip("install", ["--no-deps", "-e", kpmdir])
    @info "Python setup complete."
end

@info string("KineticaKPM.jl setup complete.")
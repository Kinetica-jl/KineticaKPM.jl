# KineticaKPM

[![CC BY-NC-SA 4.0][cc-by-nc-sa-shield]][cc-by-nc-sa]

A modular kinetic calculator addon for [Kinetica.jl](https://github.com/Kinetica-jl/Kinetica.jl), enabling fast prediction of temperature-dependent reaction rate constants.

KineticaKPM.jl uses the [KineticPredictorModel (KPM)](https://github.com/joegilkes/KineticPredictorModel) Python package to predict activation energies from the chemical structures of reactants and products. It calculates rate constants using a variety of Arrhenius-like equations, mainly using the collision theory approximation to calculate Arrhenius prefactors.

## Documentation

For information on installation, usage and development of KineticaKPM.jl, see the main documentation of Kinetica.jl (in progress).

## License

This work is licensed under a [Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License][cc-by-nc-sa]

[cc-by-nc-sa]: http://creativecommons.org/licenses/by-nc-sa/4.0/
[cc-by-nc-sa-shield]: https://img.shields.io/badge/License-CC%20BY--NC--SA%204.0-lightgrey.svg
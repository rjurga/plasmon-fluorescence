# Plasmonic modification of the decay rates of a dipole near a metal sphere

## Introduction

The environment of a quantum optical emitter can strongly affect its emission properties[[1](https://doi.org/10.1103/PhysRev.69.674)]. In particular, plasmonic environments can result in important enhancement or quenching[[2](https://doi.org/10.1103/PhysRevLett.96.113002)]. In the case of a dipole near a metal sphere, exact expressions for the decay rates are available from the electrodynamical theory[[3](https://doi.org/10.1016/0039-6028(88)90776-5)]. This program is a numerical implementation of these expressions.

## Features

The total, radiative and nonradiative decay rates as well as the quantum efficiency are computed and normalized to their values in free space. The results can be saved in a text file and plotted in an interactive window.

The expressions of the normalized total and radiative decay rates are sums over angular modes. They can be plotted as a function of the maximum angular mode to verify convergence.

## Convergence

A high number (> 100) of angular modes are required for convergence of the normalized total decay rate when the dipole is close (~ 1 nm) to the surface of the metal sphere. However the numerical evaluation of the expressions may fail if the angular modes are too high. To allow reaching higher angular modes, numerical values of the Spherical Bessel functions are converted to NumPy's clongdouble type. Please read [NumPy's documentation](https://docs.scipy.org/doc/numpy/user/basics.types.html#extended-precision) for further information about compatibility with your hardware.

## Metals

Currently only the Drude model for the description of the metal is implemented. Experimental permittivity data from scientific publications are subject to licenses that may prevent distribution with this software. However additional permittivities can be easily implemented by the user in the data_processing.permittivity function.

## References

[[1](https://doi.org/10.1103/PhysRev.69.674)] E. M. Purcell, Spontaneous Emission Probabilities at Radio Frequencies, Phys. Rev. **69**, 674 – Published 1 June 1946.

[[2](https://doi.org/10.1103/PhysRevLett.96.113002)] Pascal Anger, Palash Bharadwaj, and Lukas Novotny, Enhancement and Quenching of Single-Molecule Fluorescence, Phys. Rev. Lett. **96**, 113002 – Published 21 March 2006.

[[3](https://doi.org/10.1016/0039-6028(88)90776-5)] Young Sik Kim, P.T. Leung, Thomas F. George, Classical decay rates for molecules in the presence of a spherical surface: A complete treatment, Surface Science, Volume 195, Issue 1, 1988, Pages 1-14, ISSN 0039-6028.

# Installation and usage

## Installation

Download or clone the repository through the GitHub webpage or with your git client.

## Usage

The file `parameters.py` contains all the parameters that the user can change and their descriptions. Then run the `data_processing.py` script with the following command:
```python data_processing.py```

## Requirements

The software requirements below show the versions used for development and testing. It may work on older versions without guarantee.

- Python 3.6.1
- NumPy 1.12.1
- SciPy 0.19.0
- Matplotlib 2.0.2

# Contributing

If you wish to improve this software, you can:

- Report bugs through GitHub.
- Request features through GitHub.
- Fork the repository to implement changes and create a pull request.

Although the goals of the software are narrow and well-defined, features and pull requests that fit in this scope will be considered.

# License

This software is provided under the MIT License. Please read the [full license](LICENSE) for more information.
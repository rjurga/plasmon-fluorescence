# Decay rates of a dipole near a metal sphere

## Introduction

The environment of a quantum optical emitter can strongly affect its emission properties<sup>[[1](https://doi.org/10.1103/PhysRev.69.674)]</sup>. In particular, plasmonic environments can result in an important enhancement or quenching<sup>[[2](https://doi.org/10.1103/PhysRevLett.96.113002)]</sup>. In the case of a dipole near a metal sphere, exact expressions for the decay rates are available from the electrodynamical theory<sup>[[3](https://doi.org/10.1016/0039-6028(88)90776-5)]</sup>. This can be extended to take into account the nonlocal response of the metal<sup>[[4](https://doi.org/10.1103/PhysRevLett.31.1434), [5](https://dx.doi.org/10.1103/PhysRevB.11.2871), [6](https://doi.org/10.1021/jp204261u), [7](https://doi.org/10.1021/nn406153k), [8](https://doi.org/10.1088/0953-8984/27/18/183204)]</sup>. This software is a numerical implementation of these expressions.

## Features

Computation of:
- Total decay rate.
- Radiative decay rate.
- Nonradiative decay rate.
- Quantum efficiency.

These quantities are normalized to their values in free space.

Metal dielectric functions can be given by:
- The Drude model.
- Experimental data imported from text files.

Models of the optical response of the metal:
- Local.
- Nonlocal:
    * Hydrodynamic Theory<sup>[[9](https://doi.org/10.1002/cphc.201200992)]</sup>.
    * Generalized Nonlocal Optical Response<sup>[[10](https://doi.org/10.1038/ncomms4809)]</sup>.

The results can be:
- Saved in a text file.
- Plotted in an interactive window:
    * 1D plots as a function of the emission frequency with a fixed dipole-sphere distance.
    * 1D plots as a function of the dipole-sphere distance with a fixed emission frequency.
    * 2D map plots as a function of the emission frequency and the dipole-sphere distance.
    * 1D plots as a function of the maximum angular mode to verify convergence.

## Metals

Metals are described through their dielectric functions. The following ones are implemented:
- The Drude model.
- Evaporated gold as measured by Olmon *et al*.\*<sup>[[11](https://doi.org/10.1103/PhysRevB.86.235147)]</sup>
- Template-stripped gold as measured by Olmon *et al*.\*<sup>[[11](https://doi.org/10.1103/PhysRevB.86.235147)]</sup>
- Single-crystal gold as measured by Olmon *et al*.\*<sup>[[11](https://doi.org/10.1103/PhysRevB.86.235147)]</sup>
- Silver as measured by Yang *et al*.\*<sup>[[12](https://doi.org/10.1103/PhysRevB.91.235137)]</sup>

Additional dielectric functions can be easily implemented by the user in the `data_processing.permittivity` function by imitating the provided implementations.

\* The files containing the data of the dielectric functions are not provided. The user is required to obtain them from the linked references and place them in the `Metals` directory.

## Convergence

The expressions of the normalized total and radiative decay rates are sums over angular modes. A high number (> 100) of angular modes are required for convergence of the normalized total decay rate when the dipole is close (~ 1 nm) to the surface of the metal sphere. However the numerical evaluation of the expressions may fail if the angular modes orders are too high. To allow reaching higher angular modes orders, numerical values of the Spherical Bessel functions are converted to NumPy's clongdouble type. Please read [NumPy's documentation](https://docs.scipy.org/doc/numpy/user/basics.types.html#extended-precision) for further information about compatibility with your hardware.

On the test hardware, this software can reach up to modes of order 111. This is enough to ensure convergence in most situations except at very short dipole-sphere distances. Comparisons with finite element method calculations lead to a relative difference of no more than 3% in the worst tested situation (dipole-sphere distance of 1 nm). A much better agreement is easily reached in more favorable situations.

The software can plot the decay rates as a function of the maximum angular mode to let you verify the convergence.

## Validation

The software is validated through several tests in the `tests.py` file. The tests include:
- Comparison of Mie coefficients and special functions with numerical values from Wolfram Alpha.
- Comparison of normalized total and radiative decay rates with results of finite element method computations. With a maximum angular mode number of 111, the relative difference in the performed tests is bounded from above by 3%. This worst case is achieved when the dipole is as close as 1 nm from the surface of the metal sphere.

## References

[[1](https://doi.org/10.1103/PhysRev.69.674)] E. M. Purcell, Spontaneous Emission Probabilities at Radio Frequencies, Phys. Rev. 69, 674 – Published 1 June 1946.

[[2](https://doi.org/10.1103/PhysRevLett.96.113002)] Pascal Anger, Palash Bharadwaj, and Lukas Novotny, Enhancement and Quenching of Single-Molecule Fluorescence, Phys. Rev. Lett. 96, 113002 – Published 21 March 2006.

[[3](https://doi.org/10.1016/0039-6028(88)90776-5)] Young Sik Kim, P.T. Leung, Thomas F. George, Classical decay rates for molecules in the presence of a spherical surface: A complete treatment, Surface Science, Volume 195, Issue 1, 1988, Pages 1-14, ISSN 0039-6028.

[[4](https://doi.org/10.1103/PhysRevLett.31.1434)] R. Ruppin, Optical Properties of a Plasma Sphere, Phys. Rev. Lett. 31, 1434 – Published 10 December 1973.

[[5](https://dx.doi.org/10.1103/PhysRevB.11.2871)] R. Ruppin, Optical properties of small metal spheres, Phys. Rev. B 11, 2871 – Published 15 April 1975.

[[6](https://doi.org/10.1021/jp204261u)] Christin David and F. Javier García de Abajo, Spatial Nonlocality in the Optical Response of Metal Nanoparticles, The Journal of Physical Chemistry C 2011 115 (40), 19470-19475.

[[7](https://doi.org/10.1021/nn406153k)] Thomas Christensen, Wei Yan, Søren Raza, Antti-Pekka Jauho, N. Asger Mortensen, and Martijn Wubs, Nonlocal Response of Metallic Nanospheres Probed by Light, Electrons, and Atoms, ACS Nano 2014 8 (2), 1745-1758.

[[8](https://doi.org/10.1088/0953-8984/27/18/183204)] Søren Raza, Sergey I Bozhevolnyi, Martijn Wubs, and N. Asger Mortensen, Nonlocal optical response in metallic nanostructures, 2015 J. Phys.: Condens. Matter 27 183204.

[[9](https://doi.org/10.1002/cphc.201200992)] Ciracì, C., Pendry, J. B. and Smith, D. R., Hydrodynamic Model for Plasmonics: A Macroscopic Approach to a Microscopic Problem. ChemPhysChem, 14: 1109–1116 (2013).

[[10](https://doi.org/10.1038/ncomms4809)] N. A. Mortensen, S. Raza, M. Wubs, T. Søndergaard and S. I. Bozhevolnyi, A generalized non-local optical response theory for plasmonic nanostructures, Nature Communications 5, 3809 (2014).

[[11](https://doi.org/10.1103/PhysRevB.86.235147)] Robert L. Olmon, Brian Slovick, Timothy W. Johnson, David Shelton, Sang-Hyun Oh, Glenn D. Boreman, and Markus B. Raschke, Optical dielectric function of gold, Phys. Rev. B 86, 235147 – Published 28 December 2012.

[[12](https://doi.org/10.1103/PhysRevB.91.235137)] Honghua U. Yang, Jeffrey D'Archangel, Michael L. Sundheimer, Eric Tucker, Glenn D. Boreman, and Markus B. Raschke, Optical dielectric function of silver, Phys. Rev. B 91, 235137 – Published 22 June 2015.

# Installation and usage

## Installation

Download or clone the repository through the GitHub webpage or with your git client.

## Usage

The file `parameters.py` contains all the parameters that the user can change and their descriptions. Then run the `data_processing.py` script with the following command:
```python data_processing.py```.

## Requirements

The software requirements below show the versions used for development and testing. It may work on older versions without guarantee.

- Python 2.7.13 and 3.6.1
- NumPy 1.13.1
- SciPy 0.19.1
- Matplotlib 2.0.2

# Contributing

If you wish to improve this software, you can:

- Report bugs through GitHub.
- Request features through GitHub.
- Fork the repository to implement changes and create a pull request.

Although the goals of this software are narrow and well-defined, features and pull requests that fit in this scope will be considered.

# License

This software is provided under the MIT License. Please read the [full license](LICENSE) for more information.
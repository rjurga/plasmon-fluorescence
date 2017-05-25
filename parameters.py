# Whether to save the results in results.txt.
save = True

# Whether to plot the results in an interactive window.
show = True

# Whether to plot the decay rates as a function of the max angular mode order.
# The plot is made for the smallest given distance and emission parameter.
show_convergence = True

# Maximum angular mode order used in the computations.
n_max = 111


# Permittivity of the embedding medium.
eps_medium = 1.0

# Material of the metal sphere.
# Only one of the options below should be uncommented.
metal = 'Drude'
# metal = 'Olmon evaporated gold'
# metal = 'Olmon template-stripped gold'
# metal = 'Olmon single-crystal gold'
# metal = 'Yang silver'

# Plasma frequency in eV.
# Only used if the metal is described with Drude.
hbar_omega_p = 8.1

# Damping rate in eV.
# Only used if the metal is described with Drude.
hbar_gamma = 0.047


# Radius of the metal sphere.
radius = 30

# Unit in which the radius is given.
# Valid values: 'm', 'cm', 'mm', 'um', 'nm', 'A'.
radius_unit = 'nm'

# Orientation of the dipole with respect to the sphere.
# Valid values: 'radial', 'tangential', 'averaged'.
orientation = 'averaged'

# Quantum efficiency of the dipole in vacuum.
q_0 = 1.0


# Smallest distance between the dipole and the surface of the sphere.
distance_min = 1.0

# Largest distance between the dipole and the surface of the sphere.
distance_max = 10.0

# Number of points between the smallest and largest distances.
# If it is 1, only the smallest distance is considered.
distance_n = 20

# Unit in which the distances are given.
# Valid values: 'm', 'cm', 'mm', 'um', 'nm', 'A'.
distance_unit = 'nm'


# Smallest value of the emission parameter of the dipole.
emission_min = 1.0

# Largest value of the emission parameter of the dipole.
emission_max = 4.0

# Number of points between the smallest and largest emission parameter values.
# If it is 1, only the smallest emission parameter value is considered.
emission_n = 20

# Type of the given emission parameter.
# Valid values: 'omega', 'hbar omega (J)', 'hbar omega (eV)', 'frequency (Hz)',
#               'wavelength (m)', 'wavelength (nm)'.
emission_label = 'hbar omega (eV)'

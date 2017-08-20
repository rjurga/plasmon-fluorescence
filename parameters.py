# GENERAL PARAMETERS

params = {
    # Whether to save the results in results.txt.
    'save results': True,

    # Whether to plot the results in an interactive window.
    'show results': True,

    # Whether to plot the decay rates as a function of the max angular mode order.
    # The plot is made for the smallest given distance and emission parameter.
    'show convergence': True,

    # Maximum angular mode order used in the computations.
    'n_max': 111,
}


# MATERIALS PARAMETERS

materials = {
    # Permittivity of the embedding medium.
    'eps_medium': 1.0,

    # Material of the metal sphere.
    # Data files for Olmon and Yang should be put in the Metals directory.
    # Only one of the options below should be uncommented:
    'metal': 'Drude',
    # 'metal': 'Olmon evaporated gold',
    # 'metal': 'Olmon template-stripped gold',
    # 'metal': 'Olmon single-crystal gold',
    # 'metal': 'Yang silver',

    # Whether to enable nonlocality.
    # Only used if nonlocality is enabled or if the metal is 'Drude'.
    'nonlocal': True,

    # Permittivity contribution due to the bound response.
    # Only used if the metal is 'Drude'.
    'eps_inf': 1.0,

    # Plasma frequency in eV.
    # Only used if nonlocality is enabled or if the metal is 'Drude'.
    'hbar omega_p': 8.1,

    # Damping rate in eV.
    'hbar gamma': 0.047,

    # Fermi velocity in metres per second.
    # Only used if nonlocality is enabled.
    'v_F': 1.40e6,

    # Diffusion constant in metres squared per second.
    # Only used if nonlocality is enabled.
    # Setting D to zero uses the Hydrodynamic Theory.
    # Setting D to a higher value uses the Generalized Nonlocal Optical Response.
    'D': 8.62e-4,
}


# GEOMETRY PARAMETERS

geometry = {
    # Radius of the metal sphere.
    'radius': 30,

    # Unit in which the radius is given.
    # Valid values: 'm', 'cm', 'mm', 'um', 'nm', 'A'.
    'unit': 'nm',
}


# DIPOLE PARAMETERS

dipole = {
    # Orientation of the dipole with respect to the sphere.
    # Valid values: 'radial', 'tangential', 'averaged'.
    'orientation': 'averaged',

    # Quantum efficiency of the dipole in vacuum.
    'q_0': 1.0,
}


# DIPOLE-SPHERE DISTANCE PARAMETERS

distance = {
    # Smallest distance between the dipole and the surface of the sphere.
    'min': 1.0,

    # Largest distance between the dipole and the surface of the sphere.
    'max': 10.0,

    # Number of points between the smallest and largest distances.
    # If it is 1, only the smallest distance is considered.
    'n': 20,

    # Unit in which the distances are given.
    # Valid values: 'm', 'cm', 'mm', 'um', 'nm', 'A'.
    'unit': 'nm',
}


# EMISSION PARAMETERS

emission = {
    # Smallest value of the emission parameter of the dipole.
    'min': 1.0,

    # Largest value of the emission parameter of the dipole.
    'max': 4.0,

    # Number of points between the smallest and largest emission parameter values.
    # If it is 1, only the smallest emission parameter value is considered.
    'n': 20,

    # Type of the given emission parameter.
    # Valid values: 'omega', 'hbar omega (J)', 'hbar omega (eV)', 'frequency (Hz)',
    #               'wavelength (m)', 'wavelength (nm)'.
    'label': 'hbar omega (eV)',
}

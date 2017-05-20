import numpy as np
from scipy import constants
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

import computations
from parameters import *


def processing(save, show, n_max,
               eps_medium, metal, hbar_omega_p, hbar_gamma,
               radius, radius_unit, orientation, q_0,
               distance_min, distance_max, distance_n, distance_unit,
               emission_min, emission_max, emission_n, emission_label):
    r = convert_units(radius, radius_unit)
    distance = np.linspace(distance_min, distance_max, num=distance_n)
    d = convert_units(distance, distance_unit)
    emission = np.linspace(emission_min, emission_max, num=emission_n)
    omega = convert_emission_to_omega(emission, emission_label)
    eps_metal = permittivity(omega, metal, hbar_omega_p, hbar_gamma)
    gamma_tot, gamma_r, gamma_nr = computations.decay_rates_vectorized(n_max, eps_medium, eps_metal, omega, r, d, orientation)
    q = computations.quantum_efficiency(gamma_tot, gamma_r, q_0)
    if save:
        distance_grid, emission_grid = np.meshgrid(distance, emission)
        X = map(np.ravel, (distance_grid, emission_grid, gamma_tot, gamma_r, gamma_nr, q))
        columns = ('distance (' + distance_unit + ')',
                   emission_label,
                   'normalized total decay rate',
                   'normalized radiative decay rate',
                   'normalized nonradiative decay rate',
                   'quantum efficiency')
        np.savetxt('results.txt', np.stack(X, axis=1), header=', '.join(columns))
    if show:
        make_plot(distance, distance_n, distance_unit,
                  emission, emission_n, emission_label,
                  gamma_tot, gamma_r, gamma_nr, q)


def convert_units(x, x_unit):
    factors = {'m': 1e0,
               'cm': 1e-2,
               'mm': 1e-3,
               'um': 1e-6,
               'nm': 1e-9,
               'A': 1e-10}
    return x * factors[x_unit]


def convert_emission_to_omega(x, x_label):
    if x_label == 'omega':
        result = x
    elif x_label == 'hbar omega (J)':
        result = x / constants.hbar
    elif x_label == 'hbar omega (eV)':
        result = convert_eV_to_Hz(x)
    elif x_label == 'frequency (Hz)':
        result = 2.0 * constants.pi * x
    elif x_label == 'wavelength (m)':
        result = 2.0 * constants.pi * constants.c / x
    elif x_label == 'wavelength (nm)':
        result = 2.0 * constants.pi * constants.c / (x*1.0e-9)
    else:
        result = np.nan
    return result


def permittivity(omega, metal, hbar_omega_p, hbar_gamma):
    if metal == 'Drude':
        omega_p = convert_eV_to_Hz(hbar_omega_p)
        gamma = convert_eV_to_Hz(hbar_gamma)
        eps = 1.0 - (omega_p**2.0)/(omega*(omega + 1j*gamma))
    else:
        eps = np.nan
    return eps


def convert_eV_to_Hz(x_eV):
    return x_eV / constants.hbar * constants.eV


def make_plot(distance, distance_n, distance_unit,
              emission, emission_n, emission_label,
              gamma_tot, gamma_r, gamma_nr, q):
    plt.figure(figsize=(15, 3))
    labels = {'omega': r'$\omega$',
              'hbar omega (J)': r'$\hbar \omega$',
              'hbar omega (eV)': r'$\hbar \omega$ (eV)',
              'frequency (Hz)': r'$\nu$',
              'wavelength (m)': r'$\lambda$',
              'wavelength (nm)': r'$\lambda$ (nm)',
              'gamma_sp': r'$\gamma_\mathrm{sp} / \gamma_0$',
              'gamma_r': r'$\gamma_\mathrm{r} / \gamma_0$',
              'gamma_nr': r'$\gamma_\mathrm{nr} / \gamma_0$',
              'q': r'$q$'}
    plot_params = (
        (gamma_tot, 'gamma_sp', 'log'),
        (gamma_r, 'gamma_r', 'log'),
        (gamma_nr, 'gamma_nr', 'log'),
        (q, 'q', 'linear')
        )
    if distance_n > 1 and emission_n > 1:
        X, Y = np.meshgrid(distance, emission)
        for i, (Z, Z_label, Z_scale) in enumerate(plot_params, start=1):
            if Z_scale == 'log':
                Z_norm = LogNorm(vmin=Z.min(), vmax=Z.max())
            else:
                Z_norm=None
            plt.subplot(1, len(plot_params), i)
            plt.imshow(Z, aspect='auto', interpolation='bilinear',
                       norm=Z_norm, origin='lower',
                       extent=[X.min(), X.max(), Y.min(), Y.max()])
            plt.xlabel('distance (' + distance_unit + ')')
            plt.ylabel(labels[emission_label])
            plt.title(labels[Z_label])
            plt.colorbar()
    else:
        if emission_n == 1:
            x = distance
            x_label = 'distance (' + distance_unit + ')'
        else:
            x = emission
            x_label = emission_label
        for i, (y, y_label, y_scale) in enumerate(plot_params, start=1):
            plt.subplot(1, len(plot_params), i)
            plt.plot(x, np.ravel(y))
            plt.xlabel(labels.get(x_label, x_label))
            plt.xlim(x[0], x[-1])
            plt.ylabel(labels.get(y_label, y_label))
            plt.yscale(y_scale)
    plt.tight_layout()
    plt.show()
    plt.close()


processing(save, show, n_max,
           eps_medium, metal, hbar_omega_p, hbar_gamma,
           radius, radius_unit, orientation, q_0,
           distance_min, distance_max, distance_n, distance_unit,
           emission_min, emission_max, emission_n, emission_label)

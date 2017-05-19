import numpy as np
from scipy import constants
import matplotlib.pyplot as plt
import computations
from parameters import *


def processing(save, show, n_max,
               a, d, orientation, q_0,
               x_min, x_max, x_n, x_label,
               eps_medium, metal, hbar_omega_p, hbar_gamma):
    x = np.linspace(x_min, x_max, x_n)
    omega = convert_x_to_omega(x, x_label)
    eps_metal = permittivity(omega, metal, hbar_omega_p, hbar_gamma)
    (gamma_tot, gamma_r, gamma_nr) = computations.decay_rates_vectorized(n_max, eps_medium, eps_metal, omega, a, d, orientation)
    q = computations.quantum_efficiency(gamma_tot, gamma_r, q_0)
    if show:
        make_figure(x, x_label,
                    (gamma_tot, 'gamma_sp', 'log'),
                    (gamma_r, 'gamma_r', 'log'),
                    (gamma_nr, 'gamma_nr', 'log'),
                    (q, 'q', 'linear'))


def convert_x_to_omega(x, x_label):
    if x_label == 'omega':
        return x
    elif x_label == 'hbar omega (J)':
        return x / constants.hbar
    elif x_label == 'hbar omega (eV)':
        return convert_eV_to_Hz(x)
    elif x_label == 'frequency (Hz)':
        return 2.0 * constants.pi * x
    elif x_label == 'wavelength (m)':
        return 2.0 * constants.pi * constants.c / x
    elif x_label == 'wavelength (nm)':
        return 2.0 * constants.pi * constants.c / (x*1.0e-9)
    else:
        return np.nan


def permittivity(omega, metal, hbar_omega_p, hbar_gamma):
    eps = np.nan
    if metal == 'Drude':
        omega_p = convert_eV_to_Hz(hbar_omega_p)
        gamma = convert_eV_to_Hz(hbar_gamma)
        eps = 1.0 - (omega_p**2.0)/(omega*(omega + 1j*gamma))
    return eps


def convert_eV_to_Hz(x_eV):
    return x_eV / constants.hbar * constants.eV


def make_figure(x, x_label, *args):
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
    for (i, (y, y_label, y_scale)) in enumerate(args, start=1):
        plt.subplot(1, len(args), i)
        plt.plot(x, y)
        plt.xlabel(labels[x_label])
        plt.ylabel(labels[y_label])
        plt.yscale(y_scale)
    plt.show()


processing(save, show, n_max,
           a, d, orientation, q_0,
           x_min, x_max, x_n, x_label,
           eps_medium, metal, hbar_omega_p, hbar_gamma)

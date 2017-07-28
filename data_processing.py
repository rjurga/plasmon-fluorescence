import numpy as np
from scipy import constants
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

import computations


def processing(save, show, n_max,
               eps_medium, metal, nonlocal,
               hbar_omega_p, hbar_gamma, v_F, D,
               radius, radius_unit, orientation, q_0,
               distance_min, distance_max, distance_n, distance_unit,
               emission_min, emission_max, emission_n, emission_label):
    """Convert parameters, get decay rates, save and plot results."""
    r = convert_units(radius, radius_unit)
    distance = np.linspace(distance_min, distance_max, num=distance_n)
    d = convert_units(distance, distance_unit)
    emission = np.linspace(emission_min, emission_max, num=emission_n)
    omega = convert_emission_to_omega(emission, emission_label)
    eps_metal = permittivity(omega, metal, eps_medium, hbar_omega_p, hbar_gamma)
    eps_inf = bound_response(eps_metal, omega, hbar_omega_p, hbar_gamma)
    omega_p = convert_eV_to_Hz(hbar_omega_p)
    xi = np.sqrt(3.0/5.0*np.square(v_F) - 1j*omega*D)
    gamma_tot, gamma_r = computations.decay_rates_vectorized(n_max, nonlocal, eps_medium, eps_metal, eps_inf, omega_p, xi, omega, r, d, orientation)
    gamma_nr = computations.nonradiative_decay_rate(gamma_tot, gamma_r)
    q = computations.quantum_efficiency(gamma_tot, gamma_r, q_0)
    if save:
        save_data(distance, emission, gamma_tot, gamma_r, gamma_nr, q)
    if show:
        make_plot(distance, distance_n, distance_unit,
                  emission, emission_n, emission_label,
                  gamma_tot, gamma_r, gamma_nr, q)


def convergence(n_max, eps_medium, metal, hbar_omega_p, hbar_gamma,
                radius, radius_unit, orientation, q_0,
                distance_min, distance_unit, emission_min, emission_label):
    """Plot decay rates as a function of the max angular mode order."""
    r = convert_units(radius, radius_unit)
    d = convert_units(np.array([distance_min]), distance_unit)
    omega = convert_emission_to_omega(np.array([emission_min]), emission_label)
    eps_metal = permittivity(omega, metal, eps_medium, hbar_omega_p, hbar_gamma)
    eps_inf = bound_response(eps_metal, omega, hbar_omega_p, hbar_gamma)
    xi = np.sqrt(3.0/5.0*np.square(v_F) - 1j*omega*D)
    gamma_tot = np.empty(n_max)
    gamma_r = np.empty(n_max)
    for i, n in enumerate(range(1, n_max+1)):
        gamma_tot[i], gamma_r[i] = computations.decay_rates_vectorized(n, nonlocal, eps_medium, eps_metal, eps_inf, omega_p, xi, omega, r, d, orientation)
    plot_params = (
        (gamma_tot, r'$\gamma_\mathrm{sp} / \gamma_0$', 'linear'),
        (gamma_r, r'$\gamma_\mathrm{r} / \gamma_0$', 'linear'),
        )
    make_1d_plot(range(1, n_max+1), r'$n_\mathrm{max}$', plot_params, style='.')


def convert_units(x, x_unit):
    """Return length converted to metres."""
    factors = {'m': 1e0,
               'cm': 1e-2,
               'mm': 1e-3,
               'um': 1e-6,
               'nm': 1e-9,
               'A': 1e-10}
    return x * factors[x_unit]


def convert_emission_to_omega(x, x_label):
    """Convert emission parameter to radian per second and return omega."""
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


def permittivity(omega, metal, eps_inf, hbar_omega_p, hbar_gamma):
    """Return the permittivity at omega for the specified metal."""
    if metal == 'Drude':
        omega_p = convert_eV_to_Hz(hbar_omega_p)
        gamma = convert_eV_to_Hz(hbar_gamma)
        eps = eps_inf - np.square(omega_p)/(omega*(omega + 1j*gamma))
    elif (('Olmon' in metal) and ('gold' in metal)) or (metal == 'Yang silver'):
        params = {'Olmon evaporated gold': ('Metals/Olmon_PRB2012_EV.dat', None, 2),
                  'Olmon template-stripped gold': ('Metals/Olmon_PRB2012_TS.dat', None, 2),
                  'Olmon single-crystal gold': ('Metals/Olmon_PRB2012_SC.dat', None, 2),
                  'Yang silver': ('Metals/Ag_C_corrected.csv', ',', 1)}
        fname, d, s = params[metal]
        data = np.loadtxt(fname, delimiter=d, skiprows=s, usecols=(0,2,3))
        data = np.flipud(data)
        omega_data = convert_eV_to_Hz(data[:, 0])
        re_eps = np.interp(omega, omega_data, data[:, 1], left=np.nan, right=np.nan)
        im_eps = np.interp(omega, omega_data, data[:, 2], left=np.nan, right=np.nan)
        eps = re_eps + 1j*im_eps
    else:
        eps = np.nan
    return eps


def bound_response(eps, omega, hbar_omega_p, hbar_gamma):
    """Return the bound response at omega."""
    omega_p = convert_eV_to_Hz(hbar_omega_p)
    gamma = convert_eV_to_Hz(hbar_gamma)
    return eps + np.square(omega_p)/(omega*(omega + 1j*gamma))


def convert_eV_to_Hz(x_eV):
    """Return input converted from eV to Hz."""
    return x_eV / constants.hbar * constants.eV


def save_data(distance, emission, gamma_tot, gamma_r, gamma_nr, q):
    """Save the decay rates and quantum efficiency in results.txt."""
    distance_grid, emission_grid = np.meshgrid(distance, emission)
    X = map(np.ravel, (distance_grid, emission_grid, gamma_tot, gamma_r, gamma_nr, q))
    columns = ('distance (' + distance_unit + ')',
                emission_label,
                'normalized total decay rate',
                'normalized radiative decay rate',
                'normalized nonradiative decay rate',
                'quantum efficiency')
    np.savetxt('results.txt', np.stack(X, axis=1), header=', '.join(columns))


def make_plot(distance, distance_n, distance_unit,
              emission, emission_n, emission_label,
              gamma_tot, gamma_r, gamma_nr, q):
    """Plot the decay rates and quantum efficiency."""
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
        (gamma_tot, labels['gamma_sp'], 'log'),
        (gamma_r, labels['gamma_r'], 'log'),
        (gamma_nr, labels['gamma_nr'], 'log'),
        (q, labels['q'], 'linear')
        )
    if distance_n > 1 and emission_n > 1:
        x_label = 'distance (' + distance_unit + ')'
        y_label = labels[emission_label]
        make_2d_plot(distance, x_label, emission, y_label, plot_params)
    else:
        if emission_n == 1:
            x = distance
            x_label = 'distance (' + distance_unit + ')'
        else:
            x = emission
            x_label = labels[emission_label]
        make_1d_plot(x, x_label, plot_params)


def make_2d_plot(x, x_label, y, y_label, plot_params):
    """Make a 2d map plot of the decay rates and quantum efficiency."""
    plt.figure(figsize=(4*len(plot_params), 3))
    X, Y = np.meshgrid(x, y)
    for i, (Z, Z_label, Z_scale) in enumerate(plot_params, start=1):
        if Z_scale == 'log':
            Z_norm = LogNorm(vmin=Z.min(), vmax=Z.max())
        else:
            Z_norm=None
        plt.subplot(1, len(plot_params), i)
        plt.imshow(Z, aspect='auto', interpolation='bilinear',
                    norm=Z_norm, origin='lower',
                    extent=[X.min(), X.max(), Y.min(), Y.max()])
        plt.xlabel(x_label)
        plt.ylabel(y_label)
        plt.title(Z_label)
        plt.colorbar()
    plt.tight_layout()
    plt.show()
    plt.close()


def make_1d_plot(x, x_label, plot_params, style='-'):
    """Make a 1d plot of the decay rates and quantum efficiency."""
    plt.figure(figsize=(4*len(plot_params), 3))
    for i, (y, y_label, y_scale) in enumerate(plot_params, start=1):
        plt.subplot(1, len(plot_params), i)
        plt.plot(x, np.ravel(y), style)
        plt.xlabel(x_label)
        plt.xlim(x[0], x[-1])
        plt.ylabel(y_label)
        plt.yscale(y_scale)
    plt.tight_layout()
    plt.show()
    plt.close()


if __name__ == "__main__":
    from parameters import *
    if save or show:
        processing(save, show, n_max,
                eps_medium, metal, nonlocal,
                hbar_omega_p, hbar_gamma, v_F, D,
                radius, radius_unit, orientation, q_0,
                distance_min, distance_max, distance_n, distance_unit,
                emission_min, emission_max, emission_n, emission_label)
    if show_convergence:
        convergence(n_max, eps_medium, metal, hbar_omega_p, hbar_gamma,
                    radius, radius_unit, orientation, q_0,
                    distance_min, distance_unit, emission_min, emission_label)

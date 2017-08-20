import numpy as np
from scipy import constants
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

import computations


def processing(params, materials, geometry, dipole, distance, emission):
    """Convert parameters, get decay rates, save and plot results."""
    r = convert_units(geometry['radius'], geometry['unit'])
    d_orig = np.linspace(distance['min'], distance['max'], num=distance['n'])
    d = convert_units(d_orig, distance['unit'])
    em_orig = np.linspace(emission['min'], emission['max'], num=emission['n'])
    omega = convert_emission_to_omega(em_orig, emission['label'])
    materials['omega_p'] = convert_eV_to_Hz(materials['hbar omega_p'])
    materials['gamma'] = convert_eV_to_Hz(materials['hbar gamma'])
    eps_metal = permittivity(omega, materials)
    eps_inf = bound_response(omega, eps_metal, materials)
    gamma_tot, gamma_r = computations.decay_rates_vectorized(params['n_max'], materials['nonlocal'], materials['eps_medium'], eps_metal, eps_inf, materials['omega_p'], materials['gamma'], materials['v_F'], materials['D'], omega, r, d, dipole['orientation'])
    gamma_nr = computations.nonradiative_decay_rate(gamma_tot, gamma_r)
    q = computations.quantum_efficiency(gamma_tot, gamma_r, dipole['q_0'])
    if params['save results']:
        save_data(d_orig, distance['unit'], em_orig, emission['label'], gamma_tot, gamma_r, gamma_nr, q)
    if params['show results']:
        make_plot(d_orig, distance['n'], distance['unit'],
                  em_orig, emission['n'], emission['label'],
                  gamma_tot, gamma_r, gamma_nr, q)


def convergence(params, materials, geometry, dipole, distance, emission):
    """Plot decay rates as a function of the max angular mode order."""
    r = convert_units(geometry['radius'], geometry['unit'])
    d = convert_units(np.array([distance['min']]), distance['unit'])
    omega = convert_emission_to_omega(np.array([emission['min']]), emission['label'])
    materials['omega_p'] = convert_eV_to_Hz(materials['hbar omega_p'])
    materials['gamma'] = convert_eV_to_Hz(materials['hbar gamma'])
    eps_metal = permittivity(omega, materials)
    eps_inf = bound_response(omega, eps_metal, materials)
    gamma_tot = np.empty(params['n_max'])
    gamma_r = np.empty(params['n_max'])
    for i, n in enumerate(range(1, params['n_max']+1)):
        gamma_tot[i], gamma_r[i] = computations.decay_rates_vectorized(n, materials['nonlocal'], materials['eps_medium'], eps_metal, eps_inf, materials['omega_p'], materials['gamma'], materials['v_F'], materials['D'], omega, r, d, dipole['orientation'])
    plot_params = (
        (gamma_tot, r'$\gamma_\mathrm{sp} / \gamma_0$', 'linear'),
        (gamma_r, r'$\gamma_\mathrm{r} / \gamma_0$', 'linear'),
    )
    make_1d_plot(range(1, params['n_max']+1), r'$n_\mathrm{max}$', plot_params, style='.')


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


def permittivity(omega, materials):
    """Return the permittivity at omega for the specified metal."""

    params_Olmon_Yang = {
        # 'metal': (path to file, column delimiter, rows to skip),
        'Olmon evaporated gold': ('Metals/Olmon_PRB2012_EV.dat', None, 2),
        'Olmon template-stripped gold': ('Metals/Olmon_PRB2012_TS.dat', None, 2),
        'Olmon single-crystal gold': ('Metals/Olmon_PRB2012_SC.dat', None, 2),
        'Yang silver': ('Metals/Ag_C_corrected.csv', ',', 1),
    }

    if materials['metal'] == 'Drude':
        eps = materials['eps_inf']
        eps += free_response(omega, materials['omega_p'], materials['gamma'])

    elif materials['metal'] in params_Olmon_Yang.keys():
        fname, d, s = params_Olmon_Yang[materials['metal']]
        data = np.loadtxt(fname, delimiter=d, skiprows=s, usecols=(0,2,3))
        # flip columns such that omega is increasing
        data = np.flipud(data)
        eps = interp_permittivity(omega, data)

    # USER IMPLEMENTED PERMITTIVITY
    # fill and uncomment the block below
    # make sure to convert the frequency to the proper units
    # (omega should be in Hz, you can use convert_emission_to_omega)
    # make sure to properly order the data for interpolation
    # (omega should be increasing)

    # elif materials['metal'] == '':
    #     fname = 'Metals/'
    #     d = None
    #     s = 0
    #     cols = (0,1,2)
    #     data = np.loadtxt(fname, delimiter=d, skiprows=s, usecols=cols)
    #     eps = interp_permittivity(omega, data)

    else:
        eps = np.nan

    return eps


def interp_permittivity(omega, data):
    omega_data = convert_eV_to_Hz(data[:, 0])
    re_eps = np.interp(omega, omega_data, data[:, 1], left=np.nan, right=np.nan)
    im_eps = np.interp(omega, omega_data, data[:, 2], left=np.nan, right=np.nan)
    return re_eps + 1j*im_eps


def bound_response(omega, eps_metal, materials):
    """Return the bound response at omega."""
    eps_inf = np.copy(eps_metal)
    eps_inf -= free_response(omega, materials['omega_p'], materials['gamma'])
    return eps_inf


def free_response(omega, omega_p, gamma):
    """Return the Drude free electrons response at omega."""
    return - np.square(omega_p)/(omega*(omega + 1j*gamma))


def convert_eV_to_Hz(x_eV):
    """Return input converted from eV to Hz."""
    return x_eV / constants.hbar * constants.eV


def save_data(distance, distance_unit, emission, emission_label, gamma_tot, gamma_r, gamma_nr, q):
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
    if params['save results'] or params['show results']:
        processing(params, materials, geometry, dipole, distance, emission)
    if params['show convergence']:
        convergence(params, materials, geometry, dipole, distance, emission)

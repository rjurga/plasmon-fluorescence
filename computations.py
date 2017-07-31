import numpy as np
import scipy.special
from scipy import constants


def decay_rates_vectorized(n_max, nonlocal, eps1, eps2, eps_inf, omega_p, gamma, v_F, D, omega, r, d, orientation):
    """Return the normalized total and radiative decay rates.
    
    This function is vectorized: it iterates over omega and d.
    """
    gamma_tot = np.empty((omega.size, d.size))
    gamma_r = np.empty((omega.size, d.size))
    k1 = np.sqrt(eps1) * omega / constants.c
    k2 = np.sqrt(eps2) * omega / constants.c
    k2_nl = k_longitudinal(nonlocal, eps2, eps_inf, omega_p, gamma, v_F, D, omega)
    for i in range(omega.size):
        an, bn = mie_coefficients(n_max, nonlocal, k1[i]*r, k2[i]*r, k2_nl[i]*r, eps1, eps2[i], eps_inf[i])
        for j in range(d.size):
            gamma_tot[i, j], gamma_r[i, j] = decay_rates(n_max, an, bn, k1[i], r, d[j], orientation)
    return (gamma_tot, gamma_r)


def k_longitudinal(nonlocal, eps, eps_inf, omega_p, gamma, v_F, D, omega):
    """Return the longitudinal wavevector due to the nonlocal response."""
    if nonlocal:
        xi = np.sqrt(3.0/5.0*np.square(v_F) + D*(gamma - 1j*omega))
        k_nl = omega_p / xi * np.square(eps / eps_inf / (eps_inf-eps))
    else:
        k_nl = np.empty(omega.size) * np.nan
    return k_nl


def decay_rates(n_max, an, bn, k1, r, d, orientation):
    """Return the normalized total and radiative decay rates."""
    n = range(1, n_max + 1)
    y1 = k1 * (r+d)
    jn = scipy.special.spherical_jn(n, y1).astype(np.clongdouble)
    hn = spherical_hankel(n, y1, jn)
    if orientation == 'radial':
        return decay_rates_radial(n, bn, jn, hn, y1)
    elif orientation == 'tangential':
        return decay_rates_tangential(n, an, bn, jn, hn, y1)
    elif orientation == 'averaged':
        gamma_tot_radial, gamma_r_radial = decay_rates_radial(n, bn, jn, hn, y1)
        gamma_tot_tangential, gamma_r_tangential = decay_rates_tangential(n, an, bn, jn, hn, y1)
        gamma_tot = (gamma_tot_radial + 2.0*gamma_tot_tangential) / 3.0
        gamma_r = (gamma_r_radial + 2.0*gamma_r_tangential) / 3.0
        return (gamma_tot, gamma_r)


def decay_rates_radial(n, bn, jn, hn, y1):
    """Return the normalized decay rates for a radial dipole."""
    gamma_tot = total_decay_rate_radial(n, bn, hn, y1)
    gamma_r = radiative_decay_rate_radial(n, bn, jn, hn, y1)
    return (gamma_tot, gamma_r)


def decay_rates_tangential(n, an, bn, jn, hn, y1):
    """Return the normalized decay rates for a tangential dipole."""
    jnprime = scipy.special.spherical_jn(n, y1, derivative=True).astype(np.clongdouble)
    psinprime = psi_n_prime(y1, jn, jnprime)
    zetanprime = zeta_n_prime(n, y1, jnprime, hn)
    gamma_tot = total_decay_rate_tangential(n, an, bn, zetanprime, hn, y1)
    gamma_r = radiative_decay_rate_tangential(n, an, bn, jn, hn, psinprime, zetanprime, y1)
    return (gamma_tot, gamma_r)


def total_decay_rate_radial(n, bn, hn, y1):
    """Return the normalized total decay rate for a radial dipole.
    
    The result is the same as the Purcell factor.
    """
    terms = np.array([i * (i+1) * (2*i + 1) for i in n])
    terms = terms * bn * np.square(hn/y1)
    return 1.0 + 1.5*np.real(np.sum(terms))


def total_decay_rate_tangential(n, an, bn, zetanprime, hn, y1):
    """Return the normalized total decay rate for a tangential dipole.
    
    The result is the same as the Purcell factor.
    """
    terms = np.array([i + 0.5 for i in n])
    terms = terms * (bn*np.square(zetanprime/y1) + an*np.square(hn))
    return 1.0 + 1.5*np.real(np.sum(terms))


def radiative_decay_rate_radial(n, bn, jn, hn, y1):
    """Return the normalized radiative decay rate for a radial dipole."""
    terms = np.array([i * (i+1) * (2*i + 1) for i in n])
    terms = terms * np.square(np.absolute((jn+bn*hn)/y1))
    return 1.5 * np.sum(terms)


def radiative_decay_rate_tangential(n, an, bn, jn, hn, psinprime, zetanprime, y1):
    """Return the normalized radiative decay rate for a tangential dipole."""
    terms = np.array([2*i + 1 for i in n])
    terms = terms * (np.square(np.absolute(jn+an*hn)) + np.square(np.absolute((psinprime+bn*zetanprime)/y1)))
    return 0.75 * np.sum(terms)


def nonradiative_decay_rate(gamma_tot, gamma_r):
    """Return the normalized nonradiative decay rate."""
    return gamma_tot - gamma_r


def quantum_efficiency(gamma_tot, gamma_r, q_0):
    """Return the quantum efficiency."""
    gamma_int_0 = 1.0/q_0 - 1.0
    q = gamma_r / (gamma_tot + gamma_int_0)
    return q

def mie_coefficients(n_max, nonlocal, rho1, rho2, rho2_nl, eps1, eps2, eps_inf):
    """Return the a and b Mie coefficients."""
    n = range(1, n_max + 1)
    jn1 = scipy.special.spherical_jn(n, rho1).astype(np.clongdouble)
    jn2 = scipy.special.spherical_jn(n, rho2).astype(np.clongdouble)
    jnprime1 = scipy.special.spherical_jn(n, rho1, derivative=True).astype(np.clongdouble)
    jnprime2 = scipy.special.spherical_jn(n, rho2, derivative=True).astype(np.clongdouble)
    hn1 = spherical_hankel(n, rho1, jn1)
    psinprime1 = psi_n_prime(rho1, jn1, jnprime1)
    psinprime2 = psi_n_prime(rho2, jn2, jnprime2)
    zetanprime1 = zeta_n_prime(n, rho1, jnprime1, hn1)
    deltan = delta_n(n, rho2_nl, eps2, eps_inf, jn2) if nonlocal else 0.0
    an = mie_bn(1.0, 1.0, jn1, jn2, hn1, psinprime1, psinprime2, zetanprime1, 0.0)
    bn = mie_bn(eps1, eps2, jn1, jn2, hn1, psinprime1, psinprime2, zetanprime1, deltan)
    return (an, bn)


def mie_bn(eps1, eps2, jn1, jn2, hn1, psinprime1, psinprime2, zetanprime1, deltan):
    """Return the b Mie coefficient."""
    numerator = eps1*jn1*(psinprime2+deltan) - eps2*jn2*psinprime1
    denominator = eps2*jn2*zetanprime1 - eps1*hn1*(psinprime2+deltan)
    return numerator / denominator


def spherical_hankel(n, z, jn, derivative=False):
    """Return the spherical Hankel function of the first kind."""
    yn = scipy.special.spherical_yn(n, z, derivative=derivative).astype(np.clongdouble)
    return jn + 1j*yn


def psi_n_prime(z, jn, jnprime):
    """Return the derivative of psi."""
    return jn + z*jnprime


def zeta_n_prime(n, z, jnprime, hn):
    """Return the derivative of zeta."""
    hnprime = spherical_hankel(n, z, jnprime, derivative=True)
    return hn + z*hnprime

def delta_n(n, rho2_nl, eps2, eps_inf, jn2):
    """Return the nonlocal correction delta."""
    n = np.array(n)
    jn2_nl = scipy.special.spherical_jn(n, rho2_nl).astype(np.clongdouble)
    jnprime2_nl = scipy.special.spherical_jn(n, rho2_nl, derivative=True).astype(np.clongdouble)
    return n * (n+1) * jn2 * (eps2-eps_inf) / eps_inf * jn2_nl / (rho2_nl * jnprime2_nl)

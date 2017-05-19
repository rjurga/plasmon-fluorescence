import numpy as np
import scipy.special
from scipy import constants


def decay_rates_vectorized(n_max, eps1, eps2, omega, a, d, orientation):
    gamma_tot = np.empty(omega.shape)
    gamma_r = np.empty(omega.shape)
    for i in range(omega.size):
        (gamma_tot[i], gamma_r[i]) = decay_rates(n_max, eps1, eps2[i], omega[i], a, d, orientation)
    gamma_nr = nonradiative_decay_rate(gamma_tot, gamma_r)
    return (gamma_tot, gamma_r, gamma_nr)


def decay_rates(n_max, eps1, eps2, omega, a, d, orientation):
    n = range(1, n_max + 1)
    k1 = np.sqrt(eps1) * omega / constants.c
    k2 = np.sqrt(eps2) * omega / constants.c
    an, bn = mie_coefficients(n, k1*a, k2*a, eps1, eps2)
    kd = k1*(a+d)
    jn = scipy.special.spherical_jn(n, kd)
    hn = spherical_hankel(n, kd, jn)
    if orientation == 'radial':
        return decay_rates_radial(n, bn, jn, hn, kd)
    elif orientation == 'tangential':
        return decay_rates_tangential(n, an, bn, jn, hn, kd)
    elif orientation == 'averaged':
        gamma_tot_radial, gamma_r_radial = decay_rates_radial(n, bn, jn, hn, kd)
        gamma_tot_tangential, gamma_r_tangential = decay_rates_tangential(n, an, bn, jn, hn, kd)
        gamma_tot = (gamma_tot_radial + 2.0*gamma_tot_tangential) / 3.0
        gamma_r = (gamma_r_radial + 2.0*gamma_r_tangential) / 3.0
        return (gamma_tot, gamma_r)


def decay_rates_radial(n, bn, jn, hn, kd):
    gamma_tot = total_decay_rate_radial(n, bn, hn, kd)
    gamma_r = radiative_decay_rate_radial(n, bn, jn, hn, kd)
    return (gamma_tot, gamma_r)


def decay_rates_tangential(n, an, bn, jn, hn, kd):
    jnprime = scipy.special.spherical_jn(n, kd, derivative=True)
    psinprime = psi_n_prime(kd, jn, jnprime)
    zetanprime = zeta_n_prime(n, kd, jnprime, hn)
    gamma_tot = total_decay_rate_tangential(n, an, bn, zetanprime, hn, kd)
    gamma_r = radiative_decay_rate_tangential(n, an, bn, jn, hn, psinprime, zetanprime, kd)
    return (gamma_tot, gamma_r)


def total_decay_rate_radial(n, bn, hn, kd):
    terms = np.array([i * (i+1) * (2*i + 1) for i in n])
    terms = terms * bn * np.square(hn/kd)
    return 1.0 + 1.5*np.real(np.sum(terms))


def total_decay_rate_tangential(n, an, bn, zetanprime, hn, kd):
    terms = np.array([i + 0.5 for i in n])
    terms = terms * (bn*np.square(zetanprime/kd) + an*np.square(hn))
    return 1.0 + 1.5*np.real(np.sum(terms))


def radiative_decay_rate_radial(n, bn, jn, hn, kd):
    terms = np.array([i * (i+1) * (2*i + 1) for i in n])
    terms = terms * np.square((jn+bn*hn)/kd)
    return 1.5 * np.real(np.sum(terms))


def radiative_decay_rate_tangential(n, an, bn, jn, hn, psinprime, zetanprime, kd):
    terms = np.array([2*i + 1 for i in n])
    terms = terms * (np.square(jn+an*hn) + np.square((psinprime+bn*zetanprime)/kd))
    return 0.75 * np.real(np.sum(terms))


def nonradiative_decay_rate(gamma_tot, gamma_r):
    return gamma_tot - gamma_r


def quantum_efficiency(gamma_tot, gamma_r, q_0):
    gamma_int_0 = 1.0/q_0 - 1.0
    q = gamma_r / (gamma_tot + gamma_int_0)
    return q

def mie_coefficients(n, rho1, rho2, eps1, eps2):
    jn1 = scipy.special.spherical_jn(n, rho1)
    jn2 = scipy.special.spherical_jn(n, rho2)
    jnprime1 = scipy.special.spherical_jn(n, rho1, derivative=True)
    jnprime2 = scipy.special.spherical_jn(n, rho2, derivative=True)
    hn1 = spherical_hankel(n, rho1, jn1)
    psinprime1 = psi_n_prime(rho1, jn1, jnprime1)
    psinprime2 = psi_n_prime(rho2, jn2, jnprime2)
    zetanprime1 = zeta_n_prime(n, rho1, jnprime1, hn1)
    an = mie_bn(1.0, 1.0, jn1, jn2, hn1, psinprime1, psinprime2, zetanprime1)
    bn = mie_bn(eps1, eps2, jn1, jn2, hn1, psinprime1, psinprime2, zetanprime1)
    return (an, bn)


def mie_bn(eps1, eps2, jn1, jn2, hn1, psinprime1, psinprime2, zetanprime1):
    numerator = eps1*jn1*psinprime2 - eps2*jn2*psinprime1
    denominator = eps2*jn2*zetanprime1 - eps1*hn1*psinprime2
    return numerator / denominator


def spherical_hankel(n, z, jn, derivative=False):
    yn = scipy.special.spherical_yn(n, z, derivative=derivative)
    return jn + 1j*yn


def psi_n_prime(z, jn, jnprime):
    return jn + z*jnprime


def zeta_n_prime(n, z, jnprime, hn):
    hnprime = spherical_hankel(n, z, jnprime, derivative=True)
    return hn + z*hnprime

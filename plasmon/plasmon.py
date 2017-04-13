import numpy as np
import scipy.special


def total_decay_rate_radial(n_max, bn, hn, kd):
    terms = np.array((n * (n+1) * (2*n + 1) for n in range(1, n_max + 1)))
    terms = terms * bn * np.square(hn/kd)
    return 1.0 + 1.5*np.real(np.sum(terms))


def total_decay_rate_tangential(n_max, an, bn, zetanprime, hn, kd):
    terms = np.array((n + 0.5 for n in range(1, n_max + 1)))
    terms = terms * (bn*np.square(zetanprime/kd) + an*np.square(hn))
    return 1.0 + 1.5*np.real(np.sum(terms))


def radiative_decay_rate_radial(n_max, bn, jn, hn, kd):
    terms = np.array((n * (n+1) * (2*n + 1) for n in range(1, n_max + 1)))
    terms = terms * np.square((jn+bn*hn)/kd)
    return 1.5 * np.real(np.sum(terms))


def radiative_decay_rate_tangential(n_max, an, bn, jn, hn, psinprime, zetanprime, kd):
    terms = np.array((2*n + 1 for n in range(1, n_max + 1)))
    terms = terms * (np.square(jn+an*hn) + np.square((psinprime+bn*zetanprime)/kd))
    return 0.75 * np.real(np.sum(terms))


def nonradiative_decay_rate(gamma_tot, gamma_r):
    return gamma_tot - gamma_r


def mie_coefficients(n_max, rho1, rho2, eps1, eps2):
    jn1 = spherical_bessel(n_max, rho1)
    jn2 = spherical_bessel(n_max, rho2)
    jn1prime = spherical_bessel(n_max, rho1, d=True)
    jn2prime = spherical_bessel(n_max, rho2, d=True)
    hn1 = spherical_hankel(n_max, rho1, jn1)
    psinprime1 = psi_n_prime(rho1, jn1, jn1prime)
    psinprime2 = psi_n_prime(rho2, jn2, jn2prime)
    zetanprime1 = zeta_n_prime(n_max, rho1, jn1prime, hn1)
    an = mie_bn(1.0, 1.0, jn1, jn2, hn1, psinprime1, psinprime2, zetanprime1)
    bn = mie_bn(eps1, eps2, jn1, jn2, hn1, psinprime1, psinprime2, zetanprime1)
    return an, bn


def mie_bn(eps1, eps2, jn1, jn2, hn1, psinprime1, psinprime2, zetanprime1):
    numerator = eps1*jn1*psinprime2 - eps2*jn2*psinprime1
    denominator = eps2*jn2*zetanprime1 - eps1*hn1*psinprime2
    return numerator / denominator


def spherical_bessel(n_max, z, d=False):
    n = np.ndarray(1, n_max + 1)
    return scipy.special.spherical_jn(n, z, derivative=d)


def spherical_hankel(n_max, z, jn, d=False):
    n = np.ndarray(1, n_max + 1)
    yn = scipy.special.spherical_yn(n, z, derivative=d)
    return jn + 1j*yn


def psi_n_prime(z, jn, jnprime):
    return jn + z*jnprime


def zeta_n_prime(n_max, z, jnprime, hn):
    hnprime = spherical_hankel(n_max, z, jnprime, d=True)
    return hn + z*hnprime

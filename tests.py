import numpy as np
import scipy.special
from scipy import constants

import computations
import data_processing


def test_mie():
    """Test the functions involved in the Mie coefficients computation.

    Compare to values obtained with Wolfram Alpha.
    """
    n = 1
    a = 10e-9
    omega = 2.48 / constants.hbar * constants.eV
    eps1 = 1.0
    eps2 = -2.377 + 1j*2.856
    k1 = np.sqrt(eps1) * omega / constants.c
    k2 = np.sqrt(eps2) * omega / constants.c
    rho1 = k1*a
    rho2 = k2*a
    jn1 = scipy.special.spherical_jn(n, rho1)
    assert np.isclose(jn1, 0.041827106236)
    jn2 = scipy.special.spherical_jn(n, rho2)
    assert np.isclose(np.real(jn2), 0.034734565982)
    assert np.isclose(np.imag(jn2), 0.073239296279)
    hn1 = computations.spherical_hankel(n, rho1, jn1)
    assert np.isclose(np.real(hn1), 0.0418271062361)
    assert np.isclose(np.imag(hn1), -63.8076276033)
    jnprime1 = scipy.special.spherical_jn(n, rho1, derivative=True)
    assert np.isclose(jnprime1, 0.331755278538607)
    jnprime2 = scipy.special.spherical_jn(n, rho2, derivative=True)
    assert np.isclose(np.real(jnprime2), 0.337084148426481)
    assert np.isclose(np.imag(jnprime2), -0.004531343067615)
    psinprime1 = computations.psi_n_prime(rho1, jn1, jnprime1)
    assert np.isclose(psinprime1, 0.0835220176842857)
    psinprime2 = computations.psi_n_prime(rho2, jn2, jnprime2)
    assert np.isclose(np.real(psinprime2), 0.070389454025903)
    assert np.isclose(np.imag(psinprime2), 0.146716099221274)
    zetanprime1 = computations.zeta_n_prime(n, rho1, jnprime1, hn1)
    assert np.isclose(np.real(zetanprime1), 0.0835220176842857)
    assert np.isclose(np.imag(zetanprime1), 62.8155149095646)
    an = computations.mie_bn(1.0, 1.0, jn1, jn2, hn1, psinprime1, psinprime2, zetanprime1)
    assert np.isclose(np.real(an), -1.96544240e-6)
    assert np.isclose(np.imag(an), -2.34432221e-6)
    bn = computations.mie_bn(eps1, eps2, jn1, jn2, hn1, psinprime1, psinprime2, zetanprime1)
    assert np.isclose(np.real(bn), -0.00140178)
    assert np.isclose(np.imag(bn), 0.00149810)
    return 'Tests pass: Mie'


def finite_element_method_tests():
    r = data_processing.convert_units(30, 'nm')
    d = data_processing.convert_units(np.array([5]), 'nm')
    emission = np.linspace(1.0, 4.0, num=10)
    omega = data_processing.convert_emission_to_omega(emission, 'hbar omega (eV)')
    metal = 'Drude'
    hbar_omega_p = 8.1
    hbar_gamma = constants.hbar / (14e-15 * constants.eV)
    eps_medium = 1.0
    eps_metal = data_processing.permittivity(omega, metal, eps_medium, hbar_omega_p, hbar_gamma)
    n_max = 111
    orientation = 'radial'
    gamma_tot, gamma_r = computations.decay_rates_vectorized(n_max, eps_medium, eps_metal, omega, r, d, orientation)
    fem_data = np.loadtxt('Tests/FEM_01.txt', skiprows=10)
    assert np.allclose(emission, fem_data[:, 0])
    assert np.allclose(np.transpose(gamma_tot), fem_data[:, 1], atol=0.0, rtol=1.0e-3)
    orientation = 'tangential'
    gamma_tot, gamma_r = computations.decay_rates_vectorized(n_max, eps_medium, eps_metal, omega, r, d, orientation)
    fem_data = np.loadtxt('Tests/FEM_02.txt', skiprows=10)
    assert np.allclose(emission, fem_data[:, 0])
    assert np.allclose(np.transpose(gamma_tot), fem_data[:, 1], atol=0.0, rtol=1.0e-2)
    eps_medium = 2.0
    eps_metal = data_processing.permittivity(omega, metal, eps_medium, hbar_omega_p, hbar_gamma)
    orientation = 'radial'
    gamma_tot, gamma_r = computations.decay_rates_vectorized(n_max, eps_medium, eps_metal, omega, r, d, orientation)
    fem_data = np.loadtxt('Tests/FEM_03.txt', skiprows=10)
    assert np.allclose(emission, fem_data[:, 0])
    assert np.allclose(np.transpose(gamma_tot), fem_data[:, 1], atol=0.0, rtol=1.0e-3)
    orientation = 'tangential'
    gamma_tot, gamma_r = computations.decay_rates_vectorized(n_max, eps_medium, eps_metal, omega, r, d, orientation)
    fem_data = np.loadtxt('Tests/FEM_04.txt', skiprows=10)
    assert np.allclose(emission, fem_data[:, 0])
    assert np.allclose(np.transpose(gamma_tot), fem_data[:, 1], atol=0.0, rtol=2.0e-2)
    distance = np.linspace(1.0, 10.0, num=10)
    d = data_processing.convert_units(distance, 'nm')
    emission = 2.5
    omega = data_processing.convert_emission_to_omega(np.array([emission]), 'hbar omega (eV)')
    eps_medium = 1.0
    eps_metal = data_processing.permittivity(omega, metal, eps_medium, hbar_omega_p, hbar_gamma)
    orientation = 'radial'
    gamma_tot, gamma_r = computations.decay_rates_vectorized(n_max, eps_medium, eps_metal, omega, r, d, orientation)
    fem_data = np.loadtxt('Tests/FEM_05.txt', skiprows=10)
    assert np.allclose(d, fem_data[:, 0])
    assert np.allclose(gamma_tot, fem_data[:, 1], atol=0.0, rtol=3.0e-2)
    orientation = 'tangential'
    gamma_tot, gamma_r = computations.decay_rates_vectorized(n_max, eps_medium, eps_metal, omega, r, d, orientation)
    fem_data = np.loadtxt('Tests/FEM_06.txt', skiprows=10)
    assert np.allclose(d, fem_data[:, 0])
    assert np.allclose(gamma_tot, fem_data[:, 1], atol=0.0, rtol=3.0e-2)
    return 'Tests pass: FEM comparison'


if __name__ == "__main__":
    print(test_mie())
    print(finite_element_method_tests())

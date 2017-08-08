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
    an = computations.mie_bn(1.0, 1.0, jn1, jn2, hn1, psinprime1, psinprime2, zetanprime1, 0.0)
    assert np.isclose(np.real(an), -1.96544240e-6)
    assert np.isclose(np.imag(an), -2.34432221e-6)
    bn = computations.mie_bn(eps1, eps2, jn1, jn2, hn1, psinprime1, psinprime2, zetanprime1, 0.0)
    assert np.isclose(np.real(bn), -0.00140178)
    assert np.isclose(np.imag(bn), 0.00149810)
    return 'Tests pass: Mie local'


def test_mie_nonlocal():
    """Test the functions involved in the nonlocal Mie coefficients computation.

    Compare to values obtained with Wolfram Alpha.
    """
    n = 1
    a = 10e-9
    omega = 2.48 / constants.hbar * constants.eV
    omega_p = 8.1 / constants.hbar * constants.eV
    gamma = 0.047 / constants.hbar * constants.eV
    v_F = 1.40e6
    D = 8.62e-4
    eps1 = 1.0
    eps2 = -2.377 + 1j*2.856
    eps_inf = 1.0
    k1 = np.sqrt(eps1) * omega / constants.c
    k2 = np.sqrt(eps2) * omega / constants.c
    k2_nl = computations.k_longitudinal(True, eps2, eps_inf, omega_p, gamma, v_F, D, omega)
    assert np.isclose(np.real(k2_nl), -2.98397e9, atol=0.0, rtol=1.0e-3)
    assert np.isclose(np.imag(k2_nl), 5.25972e9, atol=0.0, rtol=1.0e-3)
    rho1 = k1*a
    rho2 = k2*a
    rho2_nl = k2_nl*a
    jn1 = scipy.special.spherical_jn(n, rho1)
    jn2 = scipy.special.spherical_jn(n, rho2)
    hn1 = computations.spherical_hankel(n, rho1, jn1)
    jnprime1 = scipy.special.spherical_jn(n, rho1, derivative=True)
    jnprime2 = scipy.special.spherical_jn(n, rho2, derivative=True)
    psinprime1 = computations.psi_n_prime(rho1, jn1, jnprime1)
    psinprime2 = computations.psi_n_prime(rho2, jn2, jnprime2)
    zetanprime1 = computations.zeta_n_prime(n, rho1, jnprime1, hn1)
    jn2_nl = scipy.special.spherical_jn(n, rho2_nl)
    assert np.isclose(np.real(jn2_nl), 5.119867079e20)
    assert np.isclose(np.imag(jn2_nl), -2.786558951e20)
    jnprime2_nl = scipy.special.spherical_jn(n, rho2_nl, derivative=True)
    assert np.isclose(np.real(jnprime2_nl), -2.7063597583597e20)
    assert np.isclose(np.imag(jnprime2_nl), -5.0690425090715e20)
    deltan = computations.delta_n(n, rho2_nl, eps2, eps_inf, jn2)
    assert np.isclose(np.real(deltan), -0.0119636)
    assert np.isclose(np.imag(deltan), 0.00117763)
    an = computations.mie_bn(1.0, 1.0, jn1, jn2, hn1, psinprime1, psinprime2, zetanprime1, 0.0)
    assert np.isclose(np.real(an), -1.96544240e-6)
    assert np.isclose(np.imag(an), -2.34432221e-6)
    bn = computations.mie_bn(eps1, eps2, jn1, jn2, hn1, psinprime1, psinprime2, zetanprime1, 0.0)
    assert np.isclose(np.real(bn), -0.00140178)
    assert np.isclose(np.imag(bn), 0.00149810)
    return 'Tests pass: Mie nonlocal'


def test_decay_rates():
    """Test the functions involved in the decay rates computation.

    Compare to values obtained with finite element methods.
    """
    r = data_processing.convert_units(30, 'nm')
    metal = 'Drude'
    eps_inf = 1.0
    hbar_omega_p = 8.1
    omega_p = data_processing.convert_eV_to_Hz(hbar_omega_p)
    hbar_gamma = constants.hbar / (14.0e-15 * constants.eV)
    gamma = data_processing.convert_eV_to_Hz(hbar_gamma)
    n_max = 111
    test_fem_local_emission_air(r, metal, eps_inf, hbar_omega_p, omega_p, hbar_gamma, gamma, n_max)
    test_fem_local_emission_dielectric(r, metal, eps_inf, hbar_omega_p, omega_p, hbar_gamma, gamma, n_max)
    test_fem_local_distance_air(r, metal, eps_inf, hbar_omega_p, omega_p, hbar_gamma, gamma, n_max)
    test_fem_nonlocal_emission_air(r, metal, eps_inf, hbar_omega_p, omega_p, hbar_gamma, gamma, n_max)
    metal = 'Olmon single-crystal gold'
    test_fem_emission_exp_eps(r, metal, eps_inf, hbar_omega_p, omega_p, hbar_gamma, gamma, n_max)
    return 'Tests pass: FEM comparison'


def test_fem_local_emission_air(r, metal, eps_inf, hbar_omega_p, omega_p, hbar_gamma, gamma, n_max):
    """Compare decay rates with FEM calculations.
    
    The comparison is for a varying emission frequency in air.
    Only the relative difference is evaluated.
    """
    nonlocal = False
    v_F = 0.0
    D = 0.0
    d = data_processing.convert_units(np.array([5]), 'nm')
    emission = np.linspace(1.0, 4.0, num=10)
    omega = data_processing.convert_emission_to_omega(emission, 'hbar omega (eV)')
    eps_medium = 1.0
    eps_metal = data_processing.permittivity(omega, metal, eps_inf, hbar_omega_p, hbar_gamma)
    eps_inf = data_processing.bound_response(eps_metal, omega, hbar_omega_p, hbar_gamma)
    orientation = 'radial'
    gamma_tot, gamma_r = computations.decay_rates_vectorized(n_max, nonlocal, eps_medium, eps_metal, eps_inf, omega_p, gamma, v_F, D, omega, r, d, orientation)
    fem_data = np.loadtxt('Tests/FEM_01.txt', skiprows=17)
    assert np.allclose(emission, fem_data[:, 0])
    assert np.allclose(np.transpose(gamma_tot), fem_data[:, 1], atol=0.0, rtol=1.0e-3)
    orientation = 'tangential'
    gamma_tot, gamma_r = computations.decay_rates_vectorized(n_max, nonlocal, eps_medium, eps_metal, eps_inf, omega_p, gamma, v_F, D, omega, r, d, orientation)
    fem_data = np.loadtxt('Tests/FEM_02.txt', skiprows=17)
    assert np.allclose(emission, fem_data[:, 0])
    assert np.allclose(np.transpose(gamma_tot), fem_data[:, 1], atol=0.0, rtol=1.0e-2)
    return 'Tests pass: FEM comparison for changing emission parameter in air with local Drude metal'


def test_fem_local_emission_dielectric(r, metal, eps_inf, hbar_omega_p, omega_p, hbar_gamma, gamma, n_max):
    """Compare decay rates with FEM calculations.
    
    The comparison is for a varying emission frequency in a dielectric medium.
    Only the relative difference is evaluated.
    """
    nonlocal = False
    v_F = 0.0
    D = 0.0
    d = data_processing.convert_units(np.array([5]), 'nm')
    emission = np.linspace(1.0, 4.0, num=10)
    omega = data_processing.convert_emission_to_omega(emission, 'hbar omega (eV)')
    gamma = data_processing.convert_eV_to_Hz(hbar_gamma)
    eps_medium = 2.0
    eps_metal = data_processing.permittivity(omega, metal, eps_inf, hbar_omega_p, hbar_gamma)
    eps_local = data_processing.bound_response(eps_metal, omega, hbar_omega_p, hbar_gamma)
    orientation = 'radial'
    gamma_tot, gamma_r = computations.decay_rates_vectorized(n_max, nonlocal, eps_medium, eps_metal, eps_local, omega_p, gamma, v_F, D, omega, r, d, orientation)
    fem_data = np.loadtxt('Tests/FEM_03.txt', skiprows=17)
    assert np.allclose(emission, fem_data[:, 0])
    assert np.allclose(np.transpose(gamma_tot), fem_data[:, 1], atol=0.0, rtol=1.0e-3)
    orientation = 'tangential'
    gamma_tot, gamma_r = computations.decay_rates_vectorized(n_max, nonlocal, eps_medium, eps_metal, eps_local, omega_p, gamma, v_F, D, omega, r, d, orientation)
    fem_data = np.loadtxt('Tests/FEM_04.txt', skiprows=17)
    assert np.allclose(emission, fem_data[:, 0])
    assert np.allclose(np.transpose(gamma_tot), fem_data[:, 1], atol=0.0, rtol=2.0e-2)
    return 'Tests pass: FEM comparison for changing emission parameter in dielectric with local Drude metal'


def test_fem_local_distance_air(r, metal, eps_inf, hbar_omega_p, omega_p, hbar_gamma, gamma, n_max):
    """Compare decay rates with FEM calculations.
    
    The comparison is for a varying distance frequency in air.
    Only the relative difference is evaluated.
    """
    nonlocal = False
    v_F = 0.0
    D = 0.0
    distance = np.linspace(1.0, 10.0, num=10)
    d = data_processing.convert_units(distance, 'nm')
    emission = 2.5
    omega = data_processing.convert_emission_to_omega(np.array([emission]), 'hbar omega (eV)')
    gamma = data_processing.convert_eV_to_Hz(hbar_gamma)
    eps_medium = 1.0
    eps_metal = data_processing.permittivity(omega, metal, eps_inf, hbar_omega_p, hbar_gamma)
    eps_local = data_processing.bound_response(eps_metal, omega, hbar_omega_p, hbar_gamma)
    orientation = 'radial'
    gamma_tot, gamma_r = computations.decay_rates_vectorized(n_max, nonlocal, eps_medium, eps_metal, eps_local, omega_p, gamma, v_F, D, omega, r, d, orientation)
    fem_data = np.loadtxt('Tests/FEM_05.txt', skiprows=17)
    assert np.allclose(d, fem_data[:, 0])
    assert np.allclose(gamma_tot, fem_data[:, 1], atol=0.0, rtol=3.0e-2)
    orientation = 'tangential'
    gamma_tot, gamma_r = computations.decay_rates_vectorized(n_max, nonlocal, eps_medium, eps_metal, eps_local, omega_p, gamma, v_F, D, omega, r, d, orientation)
    fem_data = np.loadtxt('Tests/FEM_06.txt', skiprows=17)
    assert np.allclose(d, fem_data[:, 0])
    assert np.allclose(gamma_tot, fem_data[:, 1], atol=0.0, rtol=3.0e-2)
    return 'Tests pass: FEM comparison for changing distance in air with local Drude metal'


def test_fem_nonlocal_emission_air(r, metal, eps_inf, hbar_omega_p, omega_p, hbar_gamma, gamma, n_max):
    """Compare decay rates with FEM calculations.
    
    The comparison is for a Drude metal with a nonlocal response.
    Only the relative difference is evaluated.
    """
    nonlocal = True
    v_F = 1.40e6
    D = 0.0
    d = data_processing.convert_units(np.array([5]), 'nm')
    emission = np.linspace(1.0, 4.0, num=10)
    omega = data_processing.convert_emission_to_omega(emission, 'hbar omega (eV)')
    eps_medium = 1.0
    eps_metal = data_processing.permittivity(omega, metal, eps_inf, hbar_omega_p, hbar_gamma)
    eps_local = data_processing.bound_response(eps_metal, omega, hbar_omega_p, hbar_gamma)
    orientation = 'radial'
    gamma_tot, gamma_r = computations.decay_rates_vectorized(n_max, nonlocal, eps_medium, eps_metal, eps_local, omega_p, gamma, v_F, D, omega, r, d, orientation)
    fem_data = np.loadtxt('Tests/FEM_07.txt', skiprows=17)
    assert np.allclose(emission, fem_data[:, 0])
    assert np.allclose(np.transpose(gamma_tot), fem_data[:, 1], atol=0.0, rtol=2.0e-2)
    D = 8.62e-4
    gamma_tot, gamma_r = computations.decay_rates_vectorized(n_max, nonlocal, eps_medium, eps_metal, eps_local, omega_p, gamma, v_F, D, omega, r, d, orientation)
    fem_data = np.loadtxt('Tests/FEM_08.txt', skiprows=17)
    assert np.allclose(emission, fem_data[:, 0])
    assert np.allclose(np.transpose(gamma_tot), fem_data[:, 1], atol=0.0, rtol=1.0e-2)
    return 'Tests pass: FEM comparison for changing emission parameter in air with nonlocal Drude metal'


def test_fem_emission_exp_eps(r, metal, eps_inf, hbar_omega_p, omega_p, hbar_gamma, gamma, n_max):
    """Compare decay rates with FEM calculations.
    
    The comparison is for a metal with permittivity of gold given by Olmon.
    Only the relative difference is evaluated.
    """
    emission = np.linspace(1.0, 4.0, num=10)
    omega = data_processing.convert_emission_to_omega(emission, 'hbar omega (eV)')
    orientation = 'radial'
    d = data_processing.convert_units(np.array([5]), 'nm')
    nonlocal = False
    v_F = 0.0
    D = 0.0
    eps_medium = 1.0
    eps_metal = data_processing.permittivity(omega, metal, eps_inf, hbar_omega_p, hbar_gamma)
    eps_local = data_processing.bound_response(eps_metal, omega, hbar_omega_p, hbar_gamma)
    gamma_tot, gamma_r = computations.decay_rates_vectorized(n_max, nonlocal, eps_medium, eps_metal, eps_local, omega_p, gamma, v_F, D, omega, r, d, orientation)
    fem_data = np.loadtxt('Tests/FEM_09.txt', skiprows=17)
    assert np.allclose(emission, fem_data[:, 0])
    assert np.allclose(np.transpose(gamma_tot), fem_data[:, 1], atol=0.0, rtol=1.0e-3)
    eps_medium = 2.0
    eps_metal = data_processing.permittivity(omega, metal, eps_inf, hbar_omega_p, hbar_gamma)
    eps_local = data_processing.bound_response(eps_metal, omega, hbar_omega_p, hbar_gamma)
    gamma_tot, gamma_r = computations.decay_rates_vectorized(n_max, nonlocal, eps_medium, eps_metal, eps_local, omega_p, gamma, v_F, D, omega, r, d, orientation)
    fem_data = np.loadtxt('Tests/FEM_10.txt', skiprows=17)
    assert np.allclose(emission, fem_data[:, 0])
    assert np.allclose(np.transpose(gamma_tot), fem_data[:, 1], atol=0.0, rtol=1.0e-3)
    nonlocal = True
    v_F = 1.40e6
    D = 0.0
    eps_medium = 1.0
    eps_metal = data_processing.permittivity(omega, metal, eps_inf, hbar_omega_p, hbar_gamma)
    eps_local = data_processing.bound_response(eps_metal, omega, hbar_omega_p, hbar_gamma)
    gamma_tot, gamma_r = computations.decay_rates_vectorized(n_max, nonlocal, eps_medium, eps_metal, eps_local, omega_p, gamma, v_F, D, omega, r, d, orientation)
    fem_data = np.loadtxt('Tests/FEM_11.txt', skiprows=17)
    assert np.allclose(emission, fem_data[:, 0])
    assert np.allclose(np.transpose(gamma_tot), fem_data[:, 1], atol=0.0, rtol=2.0e-2)
    eps_medium = 2.0
    eps_metal = data_processing.permittivity(omega, metal, eps_inf, hbar_omega_p, hbar_gamma)
    eps_local = data_processing.bound_response(eps_metal, omega, hbar_omega_p, hbar_gamma)
    gamma_tot, gamma_r = computations.decay_rates_vectorized(n_max, nonlocal, eps_medium, eps_metal, eps_local, omega_p, gamma, v_F, D, omega, r, d, orientation)
    fem_data = np.loadtxt('Tests/FEM_12.txt', skiprows=17)
    assert np.allclose(emission, fem_data[:, 0])
    assert np.allclose(np.transpose(gamma_tot), fem_data[:, 1], atol=0.0, rtol=2.0e-2)
    D = 8.62e-4
    eps_medium = 1.0
    eps_metal = data_processing.permittivity(omega, metal, eps_inf, hbar_omega_p, hbar_gamma)
    eps_local = data_processing.bound_response(eps_metal, omega, hbar_omega_p, hbar_gamma)
    gamma_tot, gamma_r = computations.decay_rates_vectorized(n_max, nonlocal, eps_medium, eps_metal, eps_local, omega_p, gamma, v_F, D, omega, r, d, orientation)
    fem_data = np.loadtxt('Tests/FEM_13.txt', skiprows=17)
    assert np.allclose(emission, fem_data[:, 0])
    assert np.allclose(np.transpose(gamma_tot), fem_data[:, 1], atol=0.0, rtol=1.0e-3)
    eps_medium = 2.0
    eps_metal = data_processing.permittivity(omega, metal, eps_inf, hbar_omega_p, hbar_gamma)
    eps_local = data_processing.bound_response(eps_metal, omega, hbar_omega_p, hbar_gamma)
    gamma_tot, gamma_r = computations.decay_rates_vectorized(n_max, nonlocal, eps_medium, eps_metal, eps_local, omega_p, gamma, v_F, D, omega, r, d, orientation)
    fem_data = np.loadtxt('Tests/FEM_14.txt', skiprows=17)
    assert np.allclose(emission, fem_data[:, 0])
    assert np.allclose(np.transpose(gamma_tot), fem_data[:, 1], atol=0.0, rtol=1.0e-3)
    return 'Tests pass: FEM comparison for metal with permittivity of gold given by Olmon'


if __name__ == "__main__":
    print(test_mie())
    print(test_mie_nonlocal())
    print(test_decay_rates())

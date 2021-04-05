# -*- coding: utf-8 -*-
"""

Postprocessing: calculate E-resolved microcanonical rate constants.

@author: Lei Lei
"""

import numpy as np
from scipy.special import gamma

def rho_vib_E(freq_ls, energy_step, max_energy, initial_value=None, calc_step=5):
    '''Calculate the vibrational density of states using Beyer-Swinehart direct count algorithm.
       Inputs are in order the vibrational frequencies of the harmonic oscillator, the energy step
       and the maximum energy for rovibrational density of states output, the initial vectors for
       direct count calculation, which should be the energy vector and density of states vector of
       the rotational density of states, and the energy step for direct count.
       Vibrational frequencies and energy are in the unit of cm-1.
       Output is the density of states at specific energy levels in the unit of 1/cm-1.'''

    # initialization

    E_grid = np.array(range(0, int(max_energy), int(energy_step)) + [max_energy])    # outpyt grids
    rho = np.zeros(len(E_grid[1:]))  # zero energy is not counted

    # main loop
    if initial_value == None:
        calc_grid = np.array(range(0, int(max_energy), int(calc_step)) + [max_energy])                # calculation grids
        calc_rho = np.zeros(int(max_energy / calc_step) + 1)
        calc_rho[0] = 1.
    else:
        calc_grid = initial_value[0]
        calc_rho = initial_value[1]
    # direct count
    for f in freq_ls:
        pos = np.sum(calc_grid < round(f))
        for e in calc_grid[pos:]:
            curr_loc = np.sum(calc_grid < e)
            pri_loc = np.sum(calc_grid < e-round(f))
            calc_rho[curr_loc] += calc_rho[pri_loc]

    # output
    for n, e in enumerate(E_grid[1:]):
        if n == 0:
            pos = 0
        pre_pos = pos
        pos = np.sum(calc_grid <= e)
        rho[n] = np.sum(calc_rho[pre_pos:pos]) / energy_step


    return E_grid[1:], rho

def rho_rot_E(one_D_B_ls, two_D_B_ls, one_D_sigma_ls, two_D_sigma_ls, energy_step, max_energy):
    '''Calculate the rotational density of states for molecules that are unhindered.
       Inputs are in order rotational constants for 1D and 2D rotors, the symmetry
       number for 1D and 2D rotors, the calculation energy step and the maximum energy
       to be considered. Rotational constants and energy are in the unit of cm-1.
       Output is the density of state at specified energy levels in the unit of 1/cm-1.'''
    # initialization
    E_grid = np.array(range(0, int(max_energy), int(energy_step))[1:] + [max_energy])    # outpyt grids
    rho = np.zeros(len(E_grid))

    for n, e in enumerate(E_grid):
        if len(one_D_B_ls) == 0:
            temp = np.prod([np.power(B_i * sigma_i, -1) for (B_i, sigma_i) in zip(two_D_B_ls, two_D_sigma_ls)])
        elif len(two_D_B_ls) == 0:
            temp = np.prod([np.power(B_i, -0.5) / sigma_i for (B_i, sigma_i) in zip(one_D_B_ls, one_D_sigma_ls)])
        else:
            temp_1 = np.prod([np.power(B_i, -0.5) / sigma_i for (B_i, sigma_i) in zip(one_D_B_ls, one_D_sigma_ls)])
            temp_2 = np.prod([np.power(B_i * sigma_i, -1.) for (B_i, sigma_i) in zip(two_D_B_ls, two_D_sigma_ls)])
            temp = temp_1 * temp_2

        rho[n] = np.power(np.pi, len(one_D_B_ls) / 2.) / gamma(len(two_D_B_ls) + len(one_D_B_ls) / 2.)
        rho[n] *= np.power(e, len(two_D_B_ls) + len(one_D_B_ls) / 2. - 1.) * temp

    return E_grid, rho

def rho_rtrans_E(m1, m2, energy_step, max_energy):
    '''Calculate the density of states for the relative translational motions between
       two particles. Inputs are mass of two particles in kg, the energy step and
       the maximum energy to calculate density of states. Output is the density of
       states at specific energy levels in the unit of cm-3 / cm-1.'''

    E_grid = np.array(range(0, int(max_energy), int(energy_step))[1:] + [max_energy])    # outpyt grids
    rho = np.zeros(len(E_grid))

    reduced_mass = (m1 * m2) / (m1 + m2)    # kg
    h = 6.62607004e-34    # Planck's constant, m2 kg s-1
    J_to_cm = 1.e-3 * 6.0221409e+23 / 4.184 * 349.759    # conversion factor from joule to wavenumber

    for n, e in enumerate(E_grid):
        temp = np.pi / 4. * np.power(8. * reduced_mass / h ** 2., 1.5) * np.power(e / J_to_cm, 0.5) # unit: m-3 J-1
        rho[n] = temp / (1.e6 * J_to_cm)                  # convert unit to cm-3 / cm-1

    return E_grid, rho

def convolve_rho_E(mode_1, mode_2):
    '''Convolve the rovibrational density of states with relative translational density of
       states to get the total density of states. Inputs are nested lists containing the
       energy levels and density of states for the modes to be convolved. To get a better
       result, it is optimal to have the same resolution for both modes. Output is the
       convolved energy levels and density of states in the unit of 1/cm-1.'''

    if any(mode_1[0] != mode_2[0]):
        raise Exception('The energy grids do not match.')
    E_grid = mode_1[0]
    delta_E = np.diff(E_grid)[0]
    rho_1 = np.insert(mode_1[1], 0, 1.)
    rho_2 = np.insert(mode_2[1], 0, 1.)
    rho = np.zeros(len(E_grid))

    for n, e in enumerate(E_grid):
        pos = np.sum(E_grid <= e) + 1
        temp = rho_1[:pos] * rho_2[:pos][::-1]
        rho[n] = np.sum(temp[1:] + temp[:-1]) / 2. * delta_E

    return E_grid, rho

def calculate_partition_function(E_grid, rho_E, T):
    '''Calculate the partition function based on the density of states.
       Input are the energy levels, the correspoinding density of
       states and temperature in K. Output is the Boltzmann factor
       wighted integration of density of states, which is the partition function.'''

    delta_E = E_grid[1] - E_grid[0]
    kB_cm = 0.69503476    # Boltzmann constant, cm-1 K-1

    boltz_factor = np.exp(-np.array(E_grid) / (kB_cm * T))

    q = np.sum(np.array(rho_E)[1:] * boltz_factor[1:] + np.array(rho_E)[:-1] * boltz_factor[:-1]) / 2.0 * delta_E
#    q = np.sum(np.array(rho_E) * boltz_factor) * delta_E
    return q

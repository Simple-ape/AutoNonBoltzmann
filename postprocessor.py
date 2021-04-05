# -*- coding: utf-8 -*-
"""

Postprocessing: calculate E-resolved microcanonical rate constants.

@author: Lei Lei
"""

import Mess_executor as ME
import numpy as np
from scipy.stats import norm
import os, sys

class PES():

    def __init__(self, nominal_src, nominal_MESS, conditions, hot_branching, swd, ped_output, new_dir=True):
        self.conditions = conditions  # a list in the order of temperatuer list, pressure list, and energy grid
        self.hot = hot_branching # a list in the order of hot species, staring energy level (kcal/mol), energy sapcing (kcal/mol), and number of energy levels to evaluate
        self.MESS_file = nominal_MESS
        self.ME_model = ME.Mess_Executor(nominal_src, nominal_MESS, conditions)
        self.ME_model.execute_MESS(nominal_MESS, hot_branching, ped_output, new_dir)
        self.swd = self.ME_model.swd + '/'   # species calculation directory
        self.mwd = self.ME_model.mwd  # main directory the contains all the source codes
        # for debugging
        #self.swd = swd
        #self.mwd = '/home/lab-lei/Documents/AutoNonBoltzmann/'

    def extract_hot_branching(self):
        """ Extract hot branching fractions."""
        with open(self.swd + self.MESS_file.split('.')[0]+'.log', 'r') as log:
            lines = log.readlines()
        hot_e_levels = []
        branching_ratio = []
        branching_line = 1e10
        for n, line in enumerate(lines):
            line = line.strip()
            if line.startswith('Hot distribution branching ratios'):
                branching_line = n + 2
            if n >= branching_line:
                if not line.startswith(self.hot[0]):
                    break
                else:
                    hot_e_levels.append(float(line.split()[1]))
                    branching_ratio.append([float(x) for x in line.split()[2:]])
        self.hot_e = np.array(hot_e_levels)
        self.branching = np.array(branching_ratio)
        os.chdir(self.mwd)

    def extract_ped_branching(self, energy_level, target):
        """Extract ped probabilities for target species."""
        with open(self.swd + 'ped.out', 'r') as log:
            lines = log.readlines()
        header = 'Initial well: %sInitial energy[kcal/mol] = %s' %(self.hot[0], energy_level)
        section_flag = False
        ped_e = []
        ped_probability = []
        # locate the section and extract the data
        for line in lines:
            if line.strip() == header:
                section_flag = True
                continue
            if section_flag:
                line = line.strip()
                # determine the location of target species
                if line.startswith('E'):
                    line = line.split()
                    try:
                        ind = line.index(target) - 1
                    except ValueError:
                        sys.exit("No such species %s in the PED output..." %target)
                else:
                    try:
                        line = line.split()
                        ped_e.append(float(line[0]))
                        ped_probability.append(float(line[int(ind)]))
                    except IndexError:
                        section_flag = False
                        break
        return np.array(ped_e), np.array(ped_probability)

    def ped_partition_prob(self, incoming_e, exit_e):
        """Energy partition between bimolecular species.
           Input:
                 incoming_e: list, energies for incoming complex in kcal/mol.
                 exit_e: float, energy for exit fragment in kcal/mol.
           Output:
                 result_p: list, probabilities of complexes with energies in incoming_e to form fragment with exit_e.
        """
        #!!!!! This function needs to be further modify based on statistical mechanical theories.
        #      Currently, the partition from Goldsmith et al. PCI 2015 is useda as a placeholder.

        mu = lambda E: -18.9 + 0.8 * E
        var = lambda E: 5.6 + 0.2 * E

        result_p = []
        for e in incoming_e:
            mu_e = mu(e)
            var_e = var(e)
            result_p.append(norm.pdf(exit_e, mu_e, var_e))

        return np.array(result_p)

    def ped_probability(self, target, exit_e_levels):
        """Calculate the ped probability, which is defined as the probability that complex enters with energy E and fragment exit with energy E'. Output is a second-dimensional tensor with each column being probability for complex at energy E."""

        self.extract_hot_branching()
        self.exit_e_levels = exit_e_levels
        result = np.zeros((len(exit_e_levels), len(self.hot_e)))
        # normailize ped branching
        for j, enter_e in enumerate(self.hot_e):
            ped_e, ped_prob_e = self.extract_ped_branching(enter_e, target)
            ped_prob_e = ped_prob_e / (np.sum(ped_prob_e, axis=0))
            delta_e = abs(ped_e[1] - ped_e[0])
            for i, exit_e in enumerate(exit_e_levels):
                partition_prob = self.ped_partition_prob(ped_e, exit_e)
                temp = ped_prob_e * partition_prob
                result[i,j] = np.sum(temp[1:] + temp[:-1]) / 2. * delta_e
        self.ped_prob = result
        os.chdir(self.mwd)

    def concentration(self, T, P):
        # calculate bath gas concentration in molecule cm-3, P and T should be in torr and K
        return P / (62.36359822 * T) / 1000 * 6.0221409E+23

class entrance_PES():
    """
    Entrance PES is where the non-Boltzmann rate constants are calculated.
    """

    def __init__(self, nominal_src, nominal_MESS, conditions, bimolecular_ke, fragment_density):
        self.conditions = conditions  # a list in the order of temperatuer list, pressure list, and energy grid
        self.MESS_file = nominal_MESS
        self.ME_model = ME.Mess_Executor(nominal_src, nominal_MESS, conditions)
        self.swd = self.ME_model.swd + '/'   # species calculation directory
        self.mwd = self.ME_model.mwd  # main directory the contains all the source codes

        # read user-specified E-resolved microcanonical rate cosntants for bimolecular reaction W + X, and density of states for X
        for f in [bimolecular_ke, fragment_density]:
            try:
                with open(f, 'r') as fhand:
                    lines = fhand.readlines()
                temp_e = [] # assume in cm-1
                temp = []
                for line in lines:
                    line = line.strip().split()
                    temp_e.append(float(line[0]))
                    temp.append(float(line[1]))
            except IOError:
                sys.exit('Cannot find the file: %s...' %f)
            if f == bimolecular_ke:
                self.ke_e = temp_e
                self.ke = temp
            else:
                self.density_e = temp_e
                self.density = temp
        # convert density of states into Boltzmann distribution
        T = self.pes_1.conditions[0] # temperature, K
        k_B = 0.69503476    # Boltzmann constant, cm-1 K-1
        temp = np.array(density) * np.exp(- k_B * T / np.array(density_e))
        self.boltz_dist = temp / np.sum(temp)

    def k_non_boltzmann_E(self, ped_e, ped_prob):
        """
        Calculate the E-resoved microcanonical rate constants for W + X that ultimately form P for W at energy E.
        Currently, semi-microcanonical TST method is applide. To use other method, e.g. effective temperature model,
        modify this function accordingly.
        """
        delta_e = np.abs(ped_e[1] - ped_e[0])
        result = []
        for n,e in enumerate(self.ke_e):
            # slice transition probability
            temp_prob = ped_prob[:, ped_e >= e]
            temp = np.min(len(self.boltz_dist), temp_prob.shape[1])
            result.append(np.sum(temp_prob[:,:temp] * self.boltz_dist[:temp], axis=1) * delta_e * self.ke[n])
        self.ke_non_boltz = result

    def add_fictitious_channels(self, connecting_well, well_mass, rho_w, X_mole_fraction):
        """Add fictitious channels for each non-Boltzmann reaction sequence."""
        os.chdir(self.mwd)
        # read user-specifc density of states for well in entrance well in the form of E (cm-1), rho (1/cm-1)
        try:
            with open(rho_w, 'r') as fhand:
                lines = fhand.readlines()
        except IOError:
            sys.exit('Cannot find file %s...' %rho_w)
        rho_w = []
        for line in lines:
            line = line.strip().split()
            rho_w.append(float(line[1]))
        self.rho_w = np.array(rho_w)

        # create fictitious barrier for each of the non-Botlzmann sequence
        for n in range(self.ke_non_boltz.shape[1]):
            ne = self.effective_Ne(self.ke_non_boltz[:,n], X_mole_fraction)
            ne_file = self.ME_model.add_fictitious_channel(connecting_well, well_mass)
            with open(ne_file, 'w+') as fhand:
                for i in range(len(ne)):
                    fhand.write('%s    %s\n' %(self.ke_e[i], ne[i]))
            self.ME_model.new_ne_file.append(ne_file)
        # run master equation calculation
        self.ME_model.execute_MESS(nominal_MESS, hot_branching, ped_output=False, new_dir=False)

    def effective_Ne(self, k_e, X_mole_fraction):
        """Generate effective number of states for the transition states."""
        M = self.concentration(self.conditions[0], self.conditions[1])   # bath gas concentration in molecule cm-3
        k_e_eff = k_e * M * X_mole_fraction
        h = 6.62607004e-34    # Planck's constant, m2 kg s-1
        temp = min(k_e_eff, self.rho_w)
        ne = h * self.rho_w[:temp] * k_e_eff[:temp]
        return ne


class connecting_PES(PES):
    """
    Connecting PESs are those where:
    (1) hot species R* enters PES1 and breaking into (A+B)*;
    (2) A* enters PES2, reacts with thermalized C and the complex AC eventually breaks into (D+E)*;
    (3) D* exits from PES2.

    For these PESs, we want to calculate the probability of R* entering with E that eventually form D* exiting with E'.
    """

    def __init__(self, PES_1, PES_2, fragment_density, energy_shift):
        # inherit results from the previous calculations, this enables the recursive implementation of this class for arbitrary number of coupled PESs
        self.pes_1 = PES_1
        self.pes_2 = PES_2
        self.conditions = self.pes_1.conditions # calculation conditions, i.e. T/P/E_grid
        self.hot_e = self.pes_1.hot_e # entrance energies
        # for the terminating PES, there is not exit_e_levels attribute
        try:
            self.exit_e_levels = self.pes_2.exit_e_levels # exit energies
        except AttributeError:
            pass

        # read the user-specified density of states, assuming it matches the energy grid of PAPR-MESS calculations
        try:
            with open(fragment_density, 'r') as fhand:
                lines = fhand.readlines()
            density_e = [] # assume in cm-1
            density = []
            for line in lines:
                line = line.strip().split()
                density_e.append(float(line[0]))
                density.append(float(line[1]))
        except IOError:
            sys.exit('Cannot find the file: %s...' %fragment_density)

        # convert density of states into Boltzmann distribution
        T = self.pes_1.conditions[0] # temperature, K
        k_B = 0.69503476    # Boltzmann constant, cm-1 K-1
        temp = np.array(density) * np.exp(- k_B * T / np.array(density_e))
        self.boltz_dist = temp / np.sum(temp)

        # The energy_shift defines the difference between the reference energy of two PESs.
        # For example, if exiting from A + B on PES1 and entering A + B on PES2, then the enrgy
        # shift equals the zero point energy of A + B (relative to the reference energy of PES1)
        # minus that of A + C (relative to the reference energy of PES2).
        if hasattr(self.pes_2, 'hot'):
            self.pes_2.hot[1] += energy_shift
        if hasattr(self.pes_2, 'hot_e'):
            self.pes_2.hot_e += energy_shift

    def single_transition_prob(self, incoming_ind, intermediate_ind, exit_ind):
        """Calculate the probability of a single transition that complex enters the i-th PES with incoming_e, breaks into fragment with intermediate_e, reacts with a thrid body in the (i+1)-th PES, and finally breaks into fragment with exit_e. """
        incoming_prob = self.pes_1.ped_prob[intermediate_ind, incoming_ind]
        incoming_e = self.pes_1.hot_e[incoming_ind]
        intermediate_e = self.pes_1.exit_e_levels[intermediate_ind]
        exit_e = self.pes_2.exit_e_levels[exit_ind]
        delta_e = np.abs(self.pes_1.exit_e_levels[1] - self.pes_1.exit_e_levels[0])
        # slice the right transition matrix
        temp_prob = self.pes_2.ped_prob[exit_ind, :]
        temp_prob = temp_prob[self.pes_2.hot[1] >= intermediate_e]
        temp = min(len(self.boltz_dist), len(temp_prob))
        # numerical integration
        result = np.sum(self.boltz_dist[:temp] * temp_prob[:temp]) * delta_e * incoming_prob
        return result

    def transition_prob(self):
        """Calculate the full specturm of transition probability."""
        delta_e = np.abs(self.pes_1.hot_e[1] - self.pes_1.hot_e[0])
        I, J = self.pes_1.ped_prob.shape # number of intermediate states and number of incoming states
        M, N = self.pes_2.ped_prob.shape # number of exit states and number of intermediate states
        result = np.zeros((M,J))
        for j in range(J):
            for m in range(M):
                temp = []
                for i in range(I):
                    temp.append(self.single_transition_prob(j, i, m))
                result[m,j] = np.sum(temp) * delta_e
        self.ped_prob = result

class terminating_PES(connecting_PES):
    """
    Terminating PES is the one where species R* enters, reacts with thermal R', and finally form thermal P's.
    For this PES, we want to calculate the probability of hot R* at E that fianlly form thermal P's.
    """

    def E_branching(self, incoming_ind, intermediate_ind):
        """The hot branching fractions of A* at E reacting with thermal B to form thermal products P."""
        incoming_prob = self.pes_1.ped_prob[intermediate_ind, incoming_ind]
        incoming_e = self.pes_1.hot_e[incoming_ind]
        intermediate_e = self.pes_1.exit_e_levels[intermediate_ind]
        delta_e = np.abs(self.pes_1.exit_e_levels[1] - self.pes_1.exit_e_levels[0])
        # slice the right transition matrix
        result = []
        for n in range(self.pes_2.branching.shape[1]):
            temp_prob = self.pes_2.branching[:, n]
            temp_prob = temp_prob[self.pes_2.hot_e >= intermediate_e]
            temp = min(len(self.boltz_dist), len(temp_prob))
            # numerical integration
            result.append(np.sum(self.boltz_dist[:temp] * temp_prob[:temp]) * delta_e * incoming_prob)
        return result

    def final_branching(self):
        """Calculate the full E-dependent hot branching."""
        delta_e = np.abs(self.pes_2.hot_e[1] - self.pes_2.hot_e[0])
        I, J = self.pes_1.ped_prob.shape # number of intermediate states and number of incoming states
        M, N = self.pes_2.branching.shape # number of exit states and number of species
        result = np.zeros((N,J))
        for j in range(J):
            temp = []
            for i in range(I):
                temp.append(self.E_branching(j, i))
            result[:,j] = np.sum(temp, axis=0) * delta_e
        self.branching_prob = result

####################################################################################################################################
if __name__ == '__main__':
    # global parameters for PAPR-MESS calculations
    Temperature = 1000   # K
    Pressure = 760   # Torr
    Energy_grid = 100.   # 1/cm
    conditions = [Temperature, Pressure, Energy_grid]

    model_1_src = '/home/lab-lei/Documents/AutoNonBoltzmann/C3H7+O2/'
    model_1_name = 'me1.inp'
    src_1 = '/home/lab-lei/Documents/AutoNonBoltzmann/Non-Boltzmann_calculation/calculation_3/me1/'
    hot_entrance_1 = ['W1', np.linspace(0, 30, 20)]
    hot_exit_1 = ['P1', np.linspace(-20, 40, 30)]
    model1 = PES(model_1_src, model_1_name, conditions, hot_entrance_1, swd=src_1, ped_output=True)
    model1.ped_probability(hot_exit_1[0], hot_exit_1[1])

    model_2_src = '/home/lab-lei/Documents/AutoNonBoltzmann/C3H7O2+O2/'
    model_2_name = 'me2.inp'
    src_2 = '/home/lab-lei/Documents/AutoNonBoltzmann/Non-Boltzmann_calculation/calculation_3/me2/'
    hot_entrance_2 = ['W1', np.linspace(-15, 45, 30)]
    hot_exit_2 = ['P1', np.linspace(-10, 40, 30)]
    model2 = PES(model_2_src, model_2_name, conditions, hot_entrance_2, swd=src_2, ped_output=True, new_dir=False)
    model2.ped_probability(hot_exit_2[0], hot_exit_2[1])

    ME1_ME2 = connecting_PES(model1, model2, 'rho_O2_internal_relative_prob.dat', -5.0)
    ME1_ME2.transition_prob()

    model_3_src = '/home/lab-lei/Documents/AutoNonBoltzmann/C3H7O2+O2/'
    model_3_name = 'me3.inp'
    src_3 = '/home/lab-lei/Documents/AutoNonBoltzmann/Non-Boltzmann_calculation/calculation_3/me3/'
    hot_entrance_3 = ['W1', np.linspace(0, 50, 30)]
    model3 = PES(model_3_src, model_3_name, conditions, hot_entrance_3, swd=src_3, ped_output=False, new_dir=False)
    model3.extract_hot_branching()

    Branching = terminating_PES(ME1_ME2, model3, 'rho_O2_internal_relative_prob.dat', -10.0)
    Branching.final_branching()



# -*- coding: utf-8 -*-
"""

PAPR-MESS Species Classes.

@author: Lei Lei
"""

import sys

class Mess_Input:
    """Define classes of species involved in the PAPR-MESS inputs."""
    # assign name to the class
    def name(self):
        if hasattr(self, 'Well'):
            self.name = self.Well
        elif hasattr(self, 'Bimolecular'):
            self.name = self.Bimolecular
        elif hasattr(self, 'Barrier'):
            self.name = self.Barrier

    # determine if units are attched to the key words
    def hasunit(self, attr):
        """True if unit is associated to the command."""
        value = self.__dict__[attr]
        unit_pool = ['[1/cm]', '[kcal/mol]', '[K]', '[torr]', '[angstrom]', '[amu]', '[au]', '[atm]']
        if type(value) is list:
            if value[-1] in unit_pool:
                return True
            else: return False
        elif type(value) is dict:
            if value['unit'] in unit_pool:
                return True
            else: return False
        else: return False

    def partial_match_key(self, keyword):
        """Partial match the key and get the key."""
        key_ls = self.__dict__.keys()
        for key in key_ls:
            if keyword in key:
                return key
            elif 'Nej' in keyword:
                return 'NejScale'
        return 0

    def change_Energy(self, percentage_diff):
        """Change the energy of specific species by defined percentage."""
        key = self.partial_match_key('Energy')
        org_eng = float(self.__dict__[key][0])
        unit = self.__dict__[key][1]
        if unit == '[kcal/mol]':
            new_eng = org_eng + percentage_diff / 349.759
        elif unit == '[cm-1]':
            new_eng = org_eng + percentage_diff
        else:
            sys.exit("Error: Unrecognizable unit: %s." %unit)
        self.__dict__[key][0] = str(new_eng)

    def change_Vib_Frequency(self, percentage_diff):
        """Change the vibrational frequencies of specific species by defined percentage."""
        key = 'Frequencies'
        org_fre = self.__dict__[key]['value']
        temp = []
        for f in org_fre:
            nf = f * (1. + percentage_diff)
            temp.append(nf)
        self.__dict__[key]['value'] = temp

    def change_Symmetry(self, percentage_diff):
        """Change the symmetry factor of specific species by defined percentage."""
        key = 'SymmetryFactor'
        org_sym = self.__dict__[key][0]
        new_sym = org_sym * (1. + percentage_diff)
        self.__dict__[key][0] = new_sym

    def Hindered_rotor_correction(self):
        """In PAPR-MESS code, for hindered rotor axis and symmetry have to be integers,
           otherwise it causes error."""
        target_list = ['Axis', 'Symmetry']
        for target in target_list:
            if not hasattr(self, target):
                break
            else:
                temp = []
                for x in self.__dict__[target]:
                    temp.append(int(x))
                self.__dict__[target] = temp

class Computation_Cond(Mess_Input):
    """Computational conditions of MESS."""
    def __init__(self):
        self.order = []

    def Pressure_unit(self):
        key = self.partial_match_key('Pressure')
        unit = self.__dict__[key][-1]
        return unit.strip('[]')

    def change_Temperature(self, Temp_list):
        """Change the simulation temperature list."""
        key = 'TemperatureList'
        self.__dict__[key][0] = str(Temp_list)

    def change_Pressure(self, Pres_list):
        """Change the simulation temperature list."""
        key = 'PressureList'
        unit = self.__dict__[key][1]
        if unit == '[torr]':
            self.__dict__[key][0] = str(Pres_list)
        elif unit == '[atm]':
            self.__dict__[key][0] = str([_/760. for _ in Pres_list])
        else:
            sys.exit("Error: Unrecognizable unit: %s." %unit)

    def change_energy_grid(self, energy_grid=100):
        """Cahnge the energy grid (in [1/cm]) for calculation."""
        key = 'EnergyStep'
        if hasattr(self, key):
            self.__dict__[key][0] = str(energy_grid)
        else:
            del(self.__dict__['EnergyStepOverTemperature'])
            self.__dict__[key] = [str(energy_grid), '[1/cm]']
            self.order[self.order.index('EnergyStepOverTemperature')] = 'EnergyStep'

    def hot_reaction(self, species, E_levels):
        """Calculate the hot energy branching fractions."""
        self.__dict__['HotEnergies'] = [str(species), E_levels, '[kcal/mol]']
        self.order.append('HotEnergies')

    def ped(self, output_name):
        """Calculate the PED output."""
        if not hasattr(self, 'PEDOutput'):
            self.__dict__['PEDOutput'] = [output_name]
            self.order.append('PEDOutput')

    def drop_log_output_command(self):
        cmd_ls = ['RateOutput', 'LogOutput', 'PEDSpecies']
        for x in cmd_ls:
            if hasattr(self, x):
                del(self.__dict__[x])
                self.order.pop(self.order.index(x))

# Define collision and relaxation energy transfer
class Collision_Relaxation(Mess_Input):
    """Collisional energy transfer and relaxation energy transfer."""
    def __init__(self):
        self.order = []

class Relaxtion_Exponential(Collision_Relaxation):
    def __init__(self, Factor, Power, ExponentCutoff):
        self.Factor = [str(Factor), '[1/cm]']
        self.Power = Power
        self.ExponentCutoff = ExponentCutoff

class Lennard_Jones(Collision_Relaxation):
    def __init__(self, Epsilons, Sigmas, Masses):
        self.Epsilons = [str(Epsilons), '[1/cm]']
        self.Sigmas = [str(Sigmas), '[angstrom]']
        self.Masses = [str(Masses), '[amu]']

class Bimolecular(Mess_Input):
    """Bimolecular species initilizor."""

    def __init__(self):
        self.order = []

    def get_bimolecular_pair(self):
        """Get the PAPR-MESS bimolecular species paris."""
        key_ls = self.__dict__.keys()
        fra = []
        for key in key_ls:
            if 'Fragment' in key:
                fra.append(self.__dict__[key][0])
        return "+".join(fra)


class Well(Mess_Input):
    """Well initilizor."""

    def __init__(self):
        self.order = []

class Barrier(Mess_Input):
    """Barrier initilizor."""

    def __init__(self):
        self.order = []

    def change_Img_Frequency(self, percentage_diff):
        """Change the tunneling imaginary frequencies by specified percentage."""
        org_fre = float(self.ImaginaryFrequency[0])
        new_fre = org_fre * (1. + percentage_diff)
        self.ImaginaryFrequency[0] = str(new_fre)

    def change_Nej_file(self, percentage_diff):
        """Change the Nej file by specified percentage."""
        nej_file = self.__dict__['File'][0]
        with open(nej_file, 'r') as fhand:
            base = fhand.readlines()
        for n, x in enumerate(base):
            if x.strip():
                line = x.strip().split()
                line[1] = str(float(line[1]) * (1. + percentage_diff))
                base[n] = '    '.join(line) + '\n'
        with open(nej_file, 'w') as fhand:
            fhand.writelines(base)



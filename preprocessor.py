# -*- coding: utf-8 -*-
"""
Preprocessor for PAPR-MESS input files.

Input: nominal PAPR-MESS input files
Output: species classes for stationary points on the PES

@author: Lei Lei
"""

from input_cleaner import file_cleaner
from class_generator import class_generator

class Preprocessor:
    """Parent class for preprocessing, used to obtain the species classes."""

    def __init__(self, nominal_file):
        self.nominal_file = nominal_file

    def clean_input(self):
        print("Cleaning input file for %s..." %self.nominal_file.split(".")[0])
        self.cleaned_file = file_cleaner(self.nominal_file)
        return 1

    def generate_species_classes(self):
        print("Generating PAPR-MESS classes for %s..." %self.nominal_file.split(".")[0])
        self.species_classes, self.section_order, self.files_to_copy = class_generator(self.cleaned_file)
        return 1

    def add_fictitious_channel(self, connecting_well, well_mass):
        """For non-Boltzmann reaction sequence, add a fictitious barrier that use Ne file to calculate microcanonical rate constants."""
        n = 1
        # determine the barrier number
        while True:
            fict_b = "DummyB%d" %n
            if fict_b in section_oder:
                n += 1
                continue
            else:
                self.section_oder.append(fict_b)
                break
        # add dummy Barrier and Products to the PAPR-MESS class
        self.species_classes['DP%d'%n] = {'Bimolecular': [''DummyP%d' %n']}
        self.species_classes['DP%d'%n]['Dummy'] = ['']
        self.species_classes['DP%d'%n]['order'] = ['Bimolecular', 'Dummy']
        # for dummy barrier, use ne file to calculate microcanonical rate constants
        self.species_classes["DB%d"%n] = {'Barrier':[fict_b, connecting_well, 'DummyP%d' %n]}
        self.species_classes["DB%d"%n]['Read'] = ['']
        ne_file = 'Ne_NonBoltzmann_Dummy%d.dat' %n
        self.species_classes["DB%d"%n]['File'] = ne_file
        self.species_classes["DB%d"%n]['Mass'] = [str(well_mass), '[amu]']
        self.species_classes["DB%d"%n]['SymmetryFactor'] = [1.0]
        self.species_classes["DB%d"%n]['GroundEnergy'] = self.species_classes[connecting_well].ZeroEnergy
        self.species_classes["DB%d"%n]['End'] = ['']
        self.species_classes["DB%d"%n]['order'] = ['Barrier', 'Read', 'File', 'Mass', 'SymmetryFactor', 'GroundEnergy', 'End', 'End']
        return ne_file

if __name__ == '__main__':
    nominal_species = Preprocessor('c2h2.inp')
    # nominal_species = Preprocessor('test.inp')
    nominal_species.clean_input()
    nominal_species.generate_species_classes()

# -*- coding: utf-8 -*-
"""

Write Perturbed PAPR-MESS input files.

@author: Lei Lei
"""

import preprocessor, io, os, sys, copy, shutil

class Mess_Executor:
    """Parent class for PAPR-MESS simulations."""

    def __init__(self, nominal_src, nominal_file, conditions_ls):
        self.src = nominal_src
        self.input_name = nominal_file
        self.nominal_model = preprocessor.Preprocessor(nominal_src + nominal_file)
        self.nominal_model.clean_input()
        self.nominal_model.generate_species_classes()
        self.Temp_ls = [conditions_ls[0]]          # should be a list
        self.Pres_ls = [conditions_ls[1]]          # should be a list
        self.Energy_grid = conditions_ls[2]      # should be a float
        self.new_ne_file = []

    def make_directory(self):
        """Declare the related working directories related to a calculation."""
        self.mwd = os.getcwd()    # main directory the contains all the source codes
        self.cwd = self.mwd + '/Non-Boltzmann_calculation'  # calculation working directory

    def write_class(self, output_name, model, tar_class):
        """Write a Mess_input class to file based on the order attribute."""
        # define the model and input class to be written
        curr_model = self.__dict__[model]
        curr_class = curr_model.species_classes[tar_class]
        curr_class.Hindered_rotor_correction()
        order = curr_class.__dict__["order"]
        # append a sepcific part into the input file
        fhand = io.open(output_name, 'ab')
        spec_list = ['Frequencies', 'ElectronicLevels', 'WellDepth', 'FourierExpansion']

        for k in order:
            if type(curr_class.__dict__[k]) is list:
                if k.split(' ')[0] in spec_list:
                    unit = curr_class.__dict__[k][-1]
                    value = curr_class.__dict__[k][0].split(',')
                    for x in value:
                        line = "%s%s \t\t\t\t %s\n" %(k.split(' ')[0], unit, x.strip('[]'))
                        fhand.write(line)
                elif k.split(' ')[0] == 'HotEnergies':
                    unit = curr_class.__dict__[k][-1]
                    e_levels = curr_class.__dict__[k][-2]
                    species = curr_class.__dict__[k][-3]
                    line = "%s%s \t\t\t\t %s\n" %(k.split(' ')[0], unit, len(e_levels))
                    fhand.write(line)
                    for _ in e_levels:
                        fhand.write("%s \t\t %s\n" %(species,_))
                elif curr_class.hasunit(k):
                    unit = curr_class.__dict__[k][-1]
                    value = curr_class.__dict__[k][0].replace(',', '  ')
                    line = "%s%s \t\t\t\t %s\n" %(k.split(' ')[0], unit, value.strip('[]'))
                    fhand.write(line)
                else:
                    value = str(curr_class.__dict__[k]).replace(',', '  ')
                    value = value.replace('\'', '  ')
                    temp_key = k.split(' ')[0]
                    if temp_key in ['Group', 'Axis', 'Symmetry', 'MassExpansionSize', 'PotentialExpansionSize', 'HamiltonSizeMin', 'HamiltonSizeMax', 'GridSize']:
                        temp_value = curr_class.__dict__[k]
                        line = "%s \t\t\t\t %s\n" %(temp_key, '  '.join([str(int(_)) for _ in temp_value]).replace('\'', '  '))
                    else:
                        line = "%s \t\t\t\t %s\n" %(temp_key, value.strip('[]'))
                    fhand.write(line)

            elif type(curr_class.__dict__[k]) is dict:
                if k.split(' ')[0] not in spec_list:
                    temp = copy.deepcopy(curr_class.__dict__[k])
                    unit = temp.pop('unit')
                    order = temp.pop('order')
                    value = temp
                    line = "%s%s \t\t\t\t %s\n" %(k.split(' ')[0], unit, len(value))
                    fhand.write(line)
                    for atom in order:
                        geo = str(value[atom]).replace(',', '    ')
                        line = "%s \t\t\t\t %s\n" %(atom.split(' ')[1], geo.strip('[]'))
                        fhand.write(line)
                else:
                    temp = copy.deepcopy(curr_class.__dict__[k])
                    unit = temp.pop('unit')
                    value = temp['value']
                    line = "%s%s \t\t\t\t %s\n" %(k.split(' ')[0], unit, len(value))
                    fhand.write(line)
                    for x in value:
                        tup = type(x) is tuple
                        value_x = str(x).replace(',', '   ')
                        fhand.write(value_x.strip('[]()') + '    ' + '\n' * tup)
                    fhand.write('\n' * (k.split(' ')[0] == 'Frequencies'))
        fhand.write('\n')
        fhand.close()

    def write_file(self, output_name, model, class_order):
        """Write the entire PAPR-MESS input file."""
        # initialize the file
        fhand = io.open(output_name, 'wb')
        fhand.close()
        # write classes in the given order
        for tar_class in class_order:
            self.write_class(output_name, model, tar_class)

    def new_directory(self, new_dir=True):
        """Create directory for a new calculation."""
        self.make_directory()
        # change path to the PAPR-MESS calculation directory
        if not os.path.exists(self.cwd):
            os.makedirs(self.cwd)

        # create a working directory for each calculation
        n = 1
        while True:
            trial_directory = self.cwd + '/calculation_%s' %n
            if not os.path.exists(trial_directory):
                if new_dir:
                    os.makedirs(trial_directory)
                    self.twd = trial_directory # trial working directory
                else:
                    self.twd = self.cwd + '/calculation_%s' %(n-1)
                break
            n += 1

    def new_system_file(self, filename, hot_branching, ped_output=False):
        # create the system working directory
        system_directory = self.twd + '/%s' %self.input_name.split('.')[0]
        if not os.path.exists(system_directory):
            os.makedirs(system_directory)
        os.chdir(system_directory)
        self.swd = os.getcwd()
        condition_class = self.nominal_model.species_classes['condition']
        condition_class.change_Temperature(self.Temp_ls)
        condition_class.change_Pressure(self.Pres_ls)
        condition_class.change_energy_grid(self.Energy_grid)
        condition_class.drop_log_output_command()
        if ped_output:
            condition_class.ped('ped.out')
        if hot_branching != []:
            condition_class.hot_reaction(hot_branching[0], hot_branching[1])
        self.write_file(filename, "nominal_model", self.nominal_model.section_order)

    def execute_MESS(self, filename, hot_branching=[], ped_output=False, new_dir=True):
        """Execute the PAPR-MESS."""
        self.new_directory(new_dir)

        self.new_system_file(filename, hot_branching, ped_output)
        # copy all the necesseary files
        if len(self.nominal_model.files_to_copy) > 0:
            for f in self.nominal_model.files_to_copy:
                try:
                    shutil.copy('%s/%s' %(self.src, f), '%s' %self.swd)
                except IOError:
                    sys.exit('Error: Cannot find %s in directory: %s' %(f,self.src))

        if len(self.new_ne_file) >0:
            for f in self.new_ne_file:
                try:
                    shutil.copy('%s/%s' %(self.mwd, f), '%s' %self.swd)
                    os.system('rm %s/%s '%(self.mwd, f))
                except IOError:
                    sys.exit('Error: Cannot find %s in directory: %s' %(f,self.src))

        # execute MESS
        print("Running PAPR-MESS for system %s ..." %self.input_name.split(".")[0])
        os.chdir(self.swd)
        os.system('mess %s ' %filename)

if __name__ == '__main__':
    Temperature_list = [1000]   # K
    Pressure_list = [760]   # Torr
    Energy_grid = 100.   # 1/cm
    hot_branching = ['W1', -40.0, 1, 10]   # hot species, starting energy (kcal/mol), energy spacing (kcal/mol), number of energy levels

    conditions = [Temperature_list, Pressure_list, Energy_grid]
    model = Mess_Executor('/home/lab-lei/Documents/AutoNonBoltzmann/C2H3+O2/', 'test.inp', conditions)
    model.execute_MESS('temp.inp', hot_branching, ped_output=False)



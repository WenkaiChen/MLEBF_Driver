#! /usr/bin/env python3

import os
import numpy as np
from scipy.optimize import minimize
from datetime import datetime

banner = ("######################################################\n"
          "######################################################\n"
          "###                                                ###\n"
          "### This is the python script for the optimization ###\n"
          "###     geometry structures using MLEBF method     ###\n"
          "###   (Multi-Layer Energy-Based Fragment Method)   ###\n"
          "###             written by Chen, Wen-Kai           ###\n"
          "###                                                ###\n"
          "###     This script if suitable for MBE_driver     ###\n"
          "###            Package written by Wen-Kai          ###\n"
          "###                                                ###\n"
          "###                  ATTENTION!!!                  ###\n"
          "###          THIS IS NOT A BLACK BOX!!!            ###\n"
          "###   MAKE SURE THAT YOU KNOW WHAT YOU ARE DOING   ###\n"
          "###                                                ###\n"
          "###          Email: wenkaichen@mail.bnu.edu.cn     ###\n"
          "###                                                ###\n"
          "######################################################\n"
          "######################################################\n")

#from Gen_Crdchg_mlebf import from_xyz_to_crdchg
from IO_Data_mlebf import read_coord_from_xyz_file
from IO_Data_mlebf import from_coord_np_to_crdchg
from IO_Data_mlebf import from_coord_np_to_crdchg_freeze
from IO_Data_mlebf import write_coord_np_to_xyz_file



class Optimizer_for_mlebf():
    def __init__(self,mols_args_dict):
        """
        This is a optimizer for utilizing MLEBF program package as the electronic structure method
        to optimize the geometries for a system.
        :param mols_args_dict: A dictionary. Contains many informations. Using these parameters for give parameters for
                                             function "from_coord_np_to_crdchg" and for MLEBF calculations and so on.
        """
        self.xyz_filename=mols_args_dict["init_xyz_file"]
        self.SP_calculation_cmd_0=mols_args_dict["SP_calculation_cmd_0"]
        self.gtol = mols_args_dict['gtol']
        self.mols_args_dict = mols_args_dict
        self.bool_ci = mols_args_dict['bool_ci']
        if self.bool_ci:
            self.SP_calculation_cmd_1 = mols_args_dict['SP_calculation_cmd_1']
            self.sigma = mols_args_dict['sigma']
        self.parameters = {}
        self.parameters['bohr2angstrom'] = 0.5291772

    def run(self):
        if self.mols_args_dict["Optimizer"].lower() == "scipy":
            self.run_opt_scipy()
        elif self.mols_args_dict["Optimizer"].lower() == "ase":
            self.run_opt_ase()
    def run_opt_asse(self):
        print("this function is not implemented.")
    def run_opt_scipy(self):
        init_geo_info = read_coord_from_xyz_file(self.xyz_filename)
        init_coord = init_geo_info["coordinates"] / self.parameters['bohr2angstrom']
        self.ele_list = init_geo_info["element"]

        self.iteration = 0
        self.natom = len(self.ele_list)
        if self.natom != self.mols_args_dict["natoms"]:
            print("Your {xyzfile} is not match to the natom you set in the input dictionary, check and solve it please.".format(xyzfile=self.xyz_filename))
            print("Program exit with code 1.")
            exit(1)
        self.start_time = datetime.now()
        self.last_loop_time = self.start_time
        self.outfile = open("Opt.log", "w", 1)
        self.outfile.write(banner)
        self.outfile.write("\n\nMLEBF Optimizer Runnning... \n")
        self.outfile.write("Quantum Chemistry Software Package: {QM_Package}".format(QM_Package=self.mols_args_dict["QM_Package"]))
        self.outfile.write("STARTING TIME: " + str(self.start_time) + "\n")
        # create the calculation temp dir.
        if "SP_calculation_template_0" in self.mols_args_dict:
            if not os.path.exists('State_0'):
                os.system(
                    "cp -r {SP_calculation_template_0} State_0".format(SP_calculation_template_0=self.mols_args_dict["SP_calculation_template_0"]))
            else:
                os.system('rm -r State_0')
                os.system(
                    "cp -r {SP_calculation_template_0} State_0".format(SP_calculation_template_0=self.mols_args_dict["SP_calculation_template_0"]))
        else:
            print('You need to specific the value "SP_calculation_template_0" in the input dictionary.')
            exit(1)
            '''
            if not os.path.exists('State_0'):
                os.mkdir('State_0')
            os.system('cp {mlebf_External_Scripts_Name_0} State_0/'.format(
                mlebf_External_Scripts_Name_0=self.mols_args_dict['mlebf_External_Scripts_Name_0']))
            '''
        if self.bool_ci:
            if "SP_calculation_template_1" in self.mols_args_dict:
                if not os.path.exists('State_1'):
                    os.system("cp -r {SP_calculation_template_1} State_1".format(
                        SP_calculation_template_1=self.mols_args_dict["SP_calculation_template_1"]))
                else:
                    os.system('rm -r State_1')
                    os.system("cp -r {SP_calculation_template_1} State_1".format(
                        SP_calculation_template_1=self.mols_args_dict["SP_calculation_template_1"]))
            else:
                print('You need to specific the value "SP_calculation_template_1" in the input dictionary.')
                exit(1)
                '''
                if not os.path.exists('State_1'):
                    os.mkdir('State_1')
                os.system('cp {mlebf_External_Scripts_Name_1} State_1/'.format(
                    mlebf_External_Scripts_Name_1=self.mols_args_dict['mlebf_External_Scripts_Name_1']))
                '''
        if self.mols_args_dict['QM_Package'] == 'MLEBF':
            if self.mols_args_dict["bool_single_point"]:
                self.opt_onestep_function_mlebf(init_coord)
            else:
                res = minimize(self.opt_onestep_function_mlebf, init_coord, jac=True,
                               options={'disp': True, 'gtol': self.gtol})
        elif self.mols_args_dict['QM_Package'].lower() == 'gaussian16':
            if self.mols_args_dict["bool_single_point"]:
                self.opt_onestep_function_g16(init_coord)
            else:
                res = minimize(self.opt_onestep_function_g16, init_coord, jac=True,
                               options={'disp': True, 'gtol': self.gtol})
        else:
            print('the QM_Package is not supported currently, please contact the author.\n')
            exit(2)
        self.outfile.write("\n\nELAPSED TIME: " + str(self.last_loop_time - self.start_time) + "\n")
        self.outfile.write("ENDING TIME: " + str(self.last_loop_time) + "\n")
        if self.mols_args_dict["bool_single_point"]:
            self.outfile.write("Single Point Calculation Terminated.\nBye-Bye!")
        else:
            if res.success:
                self.outfile.write("Optimization Terminated Successfully. Good Luck.\nBye-Bye!")
            else:
                self.outfile.write("Optimization Terminated Failed. Sorry.")
        self.outfile.close()

    def opt_onestep_function_mlebf(self, coord_np):
        """

        :param coord_np: In the unit of Bohr !!!! CAUTION!!!!!
        :return:
        """
        # we convert the xyz file format to crdchg file format.
        # convert coordinates from atomic unit to angstrom.
        coord_np = coord_np.reshape([-1,3]) * self.parameters['bohr2angstrom']
        if not self.mols_args_dict["bool_pbc"]:
            if self.mols_args_dict["bool_freeze"]:
                from_coord_np_to_crdchg_freeze(
                    coord_np = coord_np,
                    nmol=self.mols_args_dict["nmol"],natoms=self.mols_args_dict["natoms"],
                    mol_list=self.mols_args_dict["mol_list"],
                    dic_mol_natom=self.mols_args_dict["dic_mol_natom"],
                    dic_mol_element=self.mols_args_dict["dic_mol_element"],
                    dic_mol_chg=self.mols_args_dict["dic_mol_chg"],
                    single_mol_charge_dic=self.mols_args_dict["single_mol_charge_dic"],
                    single_mol_spin_dic=self.mols_args_dict["single_mol_spin_dic"],
                    dic_mol_index=self.mols_args_dict["dic_mol_index"],
                    active_freeze_index_list=self.mols_args_dict["active_freeze_index_list"]
                )
            else:
                from_coord_np_to_crdchg(
                    coord_np = coord_np,
                    nmol=self.mols_args_dict["nmol"],natoms=self.mols_args_dict["natoms"],
                    mol_list=self.mols_args_dict["mol_list"],
                    dic_mol_natom=self.mols_args_dict["dic_mol_natom"],
                    dic_mol_element=self.mols_args_dict["dic_mol_element"],
                    dic_mol_chg=self.mols_args_dict["dic_mol_chg"],
                    single_mol_charge_dic=self.mols_args_dict["single_mol_charge_dic"],
                    single_mol_spin_dic=self.mols_args_dict["single_mol_spin_dic"],
                    dic_mol_index=self.mols_args_dict["dic_mol_index"])
        else:
            if self.mols_args_dict["bool_freeze"]:
                from_coord_np_to_crdchg_freeze(
                    coord_np = coord_np,
                    nmol=self.mols_args_dict["nmol"],natoms=self.mols_args_dict["natoms"],
                    mol_list=self.mols_args_dict["mol_list"],
                    dic_mol_natom=self.mols_args_dict["dic_mol_natom"],
                    dic_mol_element=self.mols_args_dict["dic_mol_element"],
                    dic_mol_chg=self.mols_args_dict["dic_mol_chg"],
                    single_mol_charge_dic=self.mols_args_dict["single_mol_charge_dic"],
                    single_mol_spin_dic=self.mols_args_dict["single_mol_spin_dic"],
                    dic_mol_index=self.mols_args_dict["dic_mol_index"],
                    active_freeze_index_list=self.mols_args_dict["active_freeze_index_list"],
                    pbc=self.mols_args_dict["bool_pbc"],
                    lat_parm=self.mols_args_dict["lat_parm"]
                )
            else:
                from_coord_np_to_crdchg(
                    coord_np = coord_np,
                    nmol=self.mols_args_dict["nmol"],natoms=self.mols_args_dict["natoms"],
                    mol_list=self.mols_args_dict["mol_list"],
                    dic_mol_natom=self.mols_args_dict["dic_mol_natom"],
                    dic_mol_element=self.mols_args_dict["dic_mol_element"],
                    dic_mol_chg=self.mols_args_dict["dic_mol_chg"],
                    single_mol_charge_dic=self.mols_args_dict["single_mol_charge_dic"],
                    single_mol_spin_dic=self.mols_args_dict["single_mol_spin_dic"],
                    dic_mol_index=self.mols_args_dict["dic_mol_index"],
                    pbc=self.mols_args_dict["bool_pbc"],
                    lat_parm=self.mols_args_dict["lat_parm"])
        self.outfile.write("Step : {iteration}\n".format(iteration=self.iteration))
        self.outfile.write("\tTIME starting this step: " + str(self.last_loop_time) + "\n")
        self.outfile.write("\t@@@COORD   : {natoms}\n\t@@@COORD   : {iteration}\n".format(natoms=self.natom,iteration=self.iteration))
        for i in range(self.natom):
            self.outfile.write("\t@@@COORD   : {ele:4s}{x:15.8f}{y:15.8f}{z:15.8f}\n".format(ele=self.ele_list[i],
                                                                                             x=coord_np[i][0],
                                                                                             y=coord_np[i][1],
                                                                                             z=coord_np[i][2]))
        # now we call shell to run mlebf program: MBE_Driver
        os.system('cp cluster_0.crdchg State_0/')
        os.chdir('State_0')
        os.system(self.SP_calculation_cmd_0)
        os.chdir('..')
        if self.bool_ci:
            os.system('cp cluster_0.crdchg State_1/')
            os.chdir('State_1')
            os.system(self.SP_calculation_cmd_1)
            os.chdir('..')
        # read tot_ene, tot_grad.
        ene_0 = np.loadtxt('./State_0/tot_ene.ene')
        grad_0 = np.loadtxt('./State_0/tot_grad.grad')
        self.outfile.write("\t@@@Energy_0: {ene:20.10f} hartree\n".format(ene=np.float64(ene_0)))
        for i in range(self.natom):
            self.outfile.write("\t@@@GRAD_0  : {x:20.10f}{y:20.10f}{z:20.10f}\n".format(x=grad_0[i][0],
                                                                                        y=grad_0[i][1],
                                                                                        z=grad_0[i][2]))
        use_hess_0 = False
        if os.path.isfile('./State_0/tot_hess.hess'):
            hess_0 = np.loadtxt('./State_0/tot_hess.hess')
            use_hess_0 = True
        if self.bool_ci:
            ene_1 = np.loadtxt('./State_1/tot_ene.ene')
            grad_1= np.loadtxt('./State_1/tot_grad.grad')
            self.outfile.write("\t@@@Energy_1: {ene:20.10f} hartree\n".format(ene=np.float64(ene_1)))
            for i in range(self.natom):
                self.outfile.write("\t@@@GRAD_1  : {x:20.10f}{y:20.10f}{z:20.10f}\n".format(x=grad_1[i][0],
                                                                                            y=grad_1[i][1],
                                                                                            z=grad_1[i][2]))
            use_hess_1 = False
            if os.path.isfile('./State_1/tot_hess.hess'):
                hess_1 = np.loadtxt('./State_1/tot_hess.hess')
                use_hess_1 = True
        # get the final energy, gradients maybe hessian.
        if self.bool_ci:
            # Penalty function parameters and calculation
            alpha = 0.02
            e_mean = (ene_0 + ene_1) / 2
            e_diff = ene_1 - ene_0
            g_ij = e_diff ** 2 / (e_diff + alpha)
            ene = e_mean + self.sigma * g_ij
            grad = 0.5 * (grad_1 + grad_0) + self.sigma * \
                   ((e_diff ** 2 + 2 * alpha * e_diff) /
                    (e_diff + alpha) ** 2) * (grad_1 - grad_0)
            use_hess = False
            self.outfile.write("\t@@@Penalty function value: {value:20.10f} hartree\n".format(value=ene))
        else:
            ene = ene_0
            grad = grad_0
            use_hess = use_hess_0
            if use_hess:
                hess = hess_0
        grad = grad.reshape(-1)
        self.iteration += 1
        time = datetime.now()
        self.outfile.write("\tTIME using during this step: " + str(time-self.last_loop_time) + "\n")
        self.last_loop_time = time
        if not use_hess:
            return (ene, grad)
        else:
            return (ene, grad, hess)

    def opt_onestep_function_g16(self, coord_np):
        """

        :param coord_np: In the unit of Bohr !!!! CAUTION!!!!!
        :return:
        """
        # we convert the xyz file format to crdchg file format.
        # convert coordinates from atomic unit to angstrom.
        coord_np = coord_np.reshape([-1,3]) * self.parameters['bohr2angstrom']
        write_coord_np_to_xyz_file(coord_np,self.ele_list,'current_geom.xyz')
        self.outfile.write("Step : {iteration}\n".format(iteration=self.iteration))
        self.outfile.write("\tTIME starting this step: " + str(self.last_loop_time) + "\n")
        self.outfile.write("\t@@@COORD   : {natoms}\n\t@@@COORD   : {iteration}\n".format(natoms=self.natom,iteration=self.iteration))
        for i in range(self.natom):
            self.outfile.write("\t@@@COORD   : {ele:4s}{x:15.8f}{y:15.8f}{z:15.8f}\n".format(ele=self.ele_list[i],
                                                                                             x=coord_np[i][0],
                                                                                             y=coord_np[i][1],
                                                                                             z=coord_np[i][2]))
        # now we call shell to run mlebf program: MBE_Driver
        os.system('cp current_geom.xyz State_0/')
        os.chdir('State_0')
        os.system(self.SP_calculation_cmd_0)
        os.chdir('..')
        if self.bool_ci:
            os.system('cp current_geom.xyz State_1/')
            os.chdir('State_1')
            os.system(self.SP_calculation_cmd_1)
            os.chdir('..')
        # read tot_ene, tot_grad.
        ene_0 = np.loadtxt('./State_0/tot_ene.ene')
        grad_0 = np.loadtxt('./State_0/tot_grad.grad')
        self.outfile.write("\t@@@Energy_0: {ene:20.10f} hartree\n".format(ene=np.float64(ene_0)))
        for i in range(self.natom):
            self.outfile.write("\t@@@GRAD_0  : {x:20.10f}{y:20.10f}{z:20.10f}\n".format(x=grad_0[i][0],
                                                                                        y=grad_0[i][1],
                                                                                        z=grad_0[i][2]))
        use_hess_0 = False
        if os.path.isfile('./State_0/tot_hess.hess'):
            hess_0 = np.loadtxt('./State_0/tot_hess.hess')
            use_hess_0 = True
        if self.bool_ci:
            ene_1 = np.loadtxt('./State_1/tot_ene.ene')
            grad_1= np.loadtxt('./State_1/tot_grad.grad')
            self.outfile.write("\t@@@Energy_1: {ene:20.10f} hartree\n".format(ene=np.float64(ene_1)))
            for i in range(self.natom):
                self.outfile.write("\t@@@GRAD_1  : {x:20.10f}{y:20.10f}{z:20.10f}\n".format(x=grad_1[i][0],
                                                                                            y=grad_1[i][1],
                                                                                            z=grad_1[i][2]))
            use_hess_1 = False
            if os.path.isfile('./State_1/tot_hess.hess'):
                hess_1 = np.loadtxt('./State_1/tot_hess.hess')
                use_hess_1 = True
        # get the final energy, gradients maybe hessian.
        if self.bool_ci:
            # Penalty function parameters and calculation
            alpha = 0.02
            e_mean = (ene_0 + ene_1) / 2
            e_diff = ene_1 - ene_0
            g_ij = e_diff ** 2 / (e_diff + alpha)
            ene = e_mean + self.sigma * g_ij
            grad = 0.5 * (grad_1 + grad_0) + self.sigma * \
                   ((e_diff ** 2 + 2 * alpha * e_diff) /
                    (e_diff + alpha) ** 2) * (grad_1 - grad_0)
            use_hess = False
            self.outfile.write("\t@@@Penalty function value: {value:20.10f} hartree\n".format(value=ene))
        else:
            ene = ene_0
            grad = grad_0
            use_hess = use_hess_0
            if use_hess:
                hess = hess_0
        grad = grad.reshape(-1)
        self.iteration += 1
        time = datetime.now()
        self.outfile.write("\tTIME using during this step: " + str(time-self.last_loop_time) + "\n")
        self.last_loop_time = time
        if not use_hess:
            return (ene, grad)
        else:
            return (ene, grad, hess)

    def opt_onestep_function_bagel(self, coord_np):
        coord_np = coord_np.reshape([-1,3]) * self.parameters['bohr2angstrom']
        write_coord_np_to_xyz_file(coord_np,self.ele_list,'current_geom.xyz')
        self.outfile.write("Step : {iteration}\n".format(iteration=self.iteration))
        self.outfile.write("\tTIME starting this step: " + str(self.last_loop_time) + "\n")
        self.outfile.write("\t@@@COORD   : {natoms}\n\t@@@COORD   : {iteration}\n".format(natoms=self.natom,iteration=self.iteration))
        for i in range(self.natom):
            self.outfile.write("\t@@@COORD   : {ele:4s}{x:15.8f}{y:15.8f}{z:15.8f}\n".format(ele=self.ele_list[i],
                                                                                             x=coord_np[i][0],
                                                                                             y=coord_np[i][1],
                                                                                             z=coord_np[i][2]))
        # now we call shell to run mlebf program: MBE_Driver
        os.system('cp current_geom.xyz State_0/bagel.xyz')
        os.chdir('State_0')
        os.system(self.SP_calculation_cmd_0)
        os.chdir('..')
        if self.bool_ci:
            # different energy and gradient can be calculated at one time
            pass
        # read tot_ene, tot_grad.
        ene_0 = np.loadtxt('./State_0/tot_ene.ene')
        grad_0 = np.loadtxt('./State_0/tot_grad.grad')
        self.outfile.write("\t@@@Energy_0: {ene:20.10f} hartree\n".format(ene=np.float64(ene_0)))
        for i in range(self.natom):
            self.outfile.write("\t@@@GRAD_0  : {x:20.10f}{y:20.10f}{z:20.10f}\n".format(x=grad_0[i][0],
                                                                                        y=grad_0[i][1],
                                                                                        z=grad_0[i][2]))
        use_hess_0 = False
        if os.path.isfile('./State_0/tot_hess.hess'):
            hess_0 = np.loadtxt('./State_0/tot_hess.hess')
            use_hess_0 = True
        if self.bool_ci:
            ene_1 = np.loadtxt('./State_0/tot_ene_2.ene')
            grad_1= np.loadtxt('./State_0/tot_grad_2.grad')
            self.outfile.write("\t@@@Energy_1: {ene:20.10f} hartree\n".format(ene=np.float64(ene_1)))
            for i in range(self.natom):
                self.outfile.write("\t@@@GRAD_1  : {x:20.10f}{y:20.10f}{z:20.10f}\n".format(x=grad_1[i][0],
                                                                                            y=grad_1[i][1],
                                                                                            z=grad_1[i][2]))
            use_hess_1 = False
            if os.path.isfile('./State_0/tot_hess_2.hess'):
                hess_1 = np.loadtxt('./State_0/tot_hess_2.hess')
                use_hess_1 = True
        # get the final energy, gradients maybe hessian.
        if self.bool_ci:
            # Penalty function parameters and calculation
            alpha = 0.02
            e_mean = (ene_0 + ene_1) / 2
            e_diff = ene_1 - ene_0
            g_ij = e_diff ** 2 / (e_diff + alpha)
            ene = e_mean + self.sigma * g_ij
            grad = 0.5 * (grad_1 + grad_0) + self.sigma * \
                   ((e_diff ** 2 + 2 * alpha * e_diff) /
                    (e_diff + alpha) ** 2) * (grad_1 - grad_0)
            use_hess = False
            self.outfile.write("\t@@@Penalty function value: {value:20.10f} hartree\n".format(value=ene))
        else:
            ene = ene_0
            grad = grad_0
            use_hess = use_hess_0
            if use_hess:
                hess = hess_0
        grad = grad.reshape(-1)
        self.iteration += 1
        time = datetime.now()
        self.outfile.write("\tTIME using during this step: " + str(time-self.last_loop_time) + "\n")
        self.last_loop_time = time
        if not use_hess:
            return (ene, grad)
        else:
            return (ene, grad, hess)


if __name__ == "__main__":
    # initialize the parameters for generating crdchg file.
    mols_args_dict = {
        "init_xyz_file": "coord.xyz", # A string. the filename of xya file format
        "QM_Package": "MLEBF", # MLEBF or Gaussian16 for current version
        "Optimizer": "Scipy", # Scipy or ASE
        "SP_calculation_cmd_0": "echo -n > MBE.log; MBE_Driver -c cluster_0.crdchg -o 3 -r 20 --tot_state 2 --current_state 1 >> MBE.log",
        "mlebf_External_Scripts_Name_0": "MBE_External", # the MBE_External File for running the State_0 MLEBF program.
        "SP_calculation_template_0": "MLEBF_template_0",
        "nmol": 3, # A int, how many fragments/monomers do you have
        "natoms": 9, # A int, total number of atoms of the full system
        "mol_list": ["Wat","Wat","Wat"], # A list, you need to list all the molecular with a name. save kind of mol have the same name.
        "dic_mol_natom": {"Wat":3}, # A dictionary, how many atoms does every kind of monomer have. (Water have 3 atoms)
        "dic_mol_element": {"Wat":['O','H','H']}, # A dictionary, every kind of monomer has an element_list. Consisting with the sort of xyz file format!!!
        "dic_mol_index": {"Wat":1}, # A dictionary, every kind of monomer has an index. Writing this index after single_mol_charge & single_mol_spin
        "dic_mol_chg": {"Wat": [-0.834, 0.417, 0.417]}, # A dictionary, every atom in every kind of monomer has a charge using as the background charge.
        "single_mol_charge_dic": {"Wat": 0}, # A dictionary, total charge of the every kind of monomer. (total charge of water is 0)
        "single_mol_spin_dic": {"Wat": 1}, # A dictionary, total spin of the every kind of monomer. (total spin of water is 0)
        "bool_single_point": False, # True for single point calculation while False for Optimization
        "bool_ci": False,
        "sigma": 3.5,
        "SP_calculation_cmd_1": "echo -n > MBE.log; MBE_Driver -c cluster_0.crdchg -o 3 -r 20 --tot_state 2 --current_state 2 >> MBE.log", # only useful when bool_ci == True
        "mlebf_External_Scripts_Name_1": "MBE_External", # the MBE_External File for running the State_1 MLEBF program. only useful when bool_ci == True
        "SP_calculation_template_1": "MLEBF_template_1",
        "gtol": 1e-5, # Gradient norm must be less than gtol before successful termination.
        "bool_freeze": False, # whether we need to freeze part of the fragment or not
        "active_freeze_index_list": ['A', 'A', 'F'], # if bool_freeze == True, then you need to list all the molecular need freeze or not. 'F' means freeze while 'A' means active.
        "bool_pbc": False, # PBC_MLEBF or not
        "lat_parm": "6.0940 0.0 0.0 0.0 6.0940 0.0 0.0 0.0 6.0940", # A string, recording the lattice vector. (A string 9 float)
                                                                    # Caution: if you use this optimizer to optimizer Gaussian, 9 numbers is needed.
    }
    optimizer = Optimizer_for_mlebf(mols_args_dict=mols_args_dict)
    optimizer.run()



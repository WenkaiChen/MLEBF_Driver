#! /usr/bin/env python3
import os
import numpy as np
import sys
from datetime import datetime

banner = ("######################################################\n"
          "######################################################\n"
          "###                                                ###\n"
          "### This is the python script for the calculating  ###\n"
          "###   approximated non-adiabatic coupling vector   ###\n"
          "###                using MLEBF method              ###\n"
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

from IO_Data_mlebf import read_coord_from_xyz_file
from IO_Data_mlebf import from_coord_np_to_crdchg

def calculate_nac_SS(ene_0, ene_1, grad_0_2d, grad_1_2d, hess_0, hess_1):
    """
    Calculate the non-adiabatic coupling vectors between two singlets.
    :param ene_0: A float. Energy of the state with lower energy.
    :param ene_1: A float. Energy of the state with higher energy.
    :param grad_0_2d: A np.array. Gradients of the lower state. shape: (natom,3)
    :param grad_1_2d: A np.array. Gradients of the higher state. shape: (natom,3)
    :param hess_0: A np.array. Hessian matrix of the lower state. shape: (3*natom, 3*natom)
    :param hess_1: A np.array. Hessian matrix of the higher state. shape: (3*natom, 3*natom)
    :return: nac: A np.array. non-adiabatic coupling vector. (natom,3)
    """
    natom = grad_0_2d.shape[0]
    ene_diff = ene_1 - ene_0
    grad_0 = grad_0_2d.reshape(-1,1)
    grad_1 = grad_1_2d.reshape(-1,1)
    if ene_diff == 0:
        ene_diff = 0.0000000001
    gi_gi = np.dot(grad_0,grad_0.transpose())
    gi_gj = np.dot(grad_0,grad_1.transpose())
    gj_gj = np.dot(grad_1,grad_1.transpose())
    #gj_gi = np.dot(grad_1,grad_0.transpose())

    grad_diff = 0.5*(grad_1-grad_0)
    grad_diff2= np.dot(grad_diff,grad_diff.transpose())

    d2E2_d2r = 2*(ene_diff * (hess_1-hess_0) + gi_gi+gj_gj-2*gi_gj)
    magnitude = (d2E2_d2r/8) - grad_diff2

    # SVD
    u, s, vh = np.linalg.svd(magnitude)
    ev = vh[0]
    # get one phase
    e = max(ev[0:2].min(),ev[0:2].max(),key=abs)
    if e < 0:
        ev = -ev
    ew = s[0]
    count = 0
    hopping_direction = np.zeros(shape=(natom,3))
    for iatom in range(natom):
        for xyz in range(3):
            hopping_direction[iatom][xyz] = ev[count]
            count += 1
    hopping_magnitude = np.sqrt(ew)/ene_diff
    nac = hopping_direction * hopping_magnitude

    return nac

def calculate_nac_TT(ene_0, ene_1, grad_0_2d, grad_1_2d, hess_0, hess_1):
    """
    Calculate the non-adiabatic coupling vectors between two triplets.
    :param ene_0: A float. Energy of the state with lower energy.
    :param ene_1: A float. Energy of the state with higher energy.
    :param grad_0_2d: A np.array. Gradients of the lower state. shape: (natom,3)
    :param grad_1_2d: A np.array. Gradients of the higher state. shape: (natom,3)
    :param hess_0: A np.array. Hessian matrix of the lower state. shape: (3*natom, 3*natom)
    :param hess_1: A np.array. Hessian matrix of the higher state. shape: (3*natom, 3*natom)
    :return: nac: A np.array. non-adiabatic coupling vector. (natom*3,3*3)
    """
    natom = grad_0_2d.shape[0]
    ene_diff = ene_1 - ene_0
    grad_0 = grad_0_2d.reshape(-1,1)
    grad_1 = grad_1_2d.reshape(-1,1)
    if ene_diff == 0:
        ene_diff = 0.0000000001
    gi_gi = np.dot(grad_0,grad_0.transpose())
    gi_gj = np.dot(grad_0,grad_1.transpose())
    gj_gj = np.dot(grad_1,grad_1.transpose())
    #gj_gi = np.dot(grad_1,grad_0.transpose())

    grad_diff = 0.5*(grad_1-grad_0)
    grad_diff2= np.dot(grad_diff,grad_diff.transpose())

    d2E2_d2r = 2*(ene_diff * (hess_1-hess_0) + gi_gi+gj_gj-2*gi_gj)
    magnitude = (d2E2_d2r/8) - grad_diff2

    # SVD
    u, s, vh = np.linalg.svd(magnitude)
    ev = vh[0]
    # get one phase
    e = max(ev[0:2].min(),ev[0:2].max(),key=abs)
    if e < 0:
        ev = -ev
    ew = s[0]
    count = 0
    hopping_direction = np.zeros(shape=(natom,3))
    for iatom in range(natom):
        for xyz in range(3):
            hopping_direction[iatom][xyz] = ev[count]
            count +=1
    hopping_magnitude = np.sqrt(ew)/ene_diff
    nac_ = hopping_direction * hopping_magnitude
    nac = np.tile(nac_,[3,3])

    return nac

def _main_SP(mols_args_dict):
    print(banner)
    xyz_filename=mols_args_dict["init_xyz_file"]
    # read the input coordinate.
    init_geo_info = read_coord_from_xyz_file(xyz_filename)
    init_coord = init_geo_info["coordinates"] # coordinates in unit of angstrom
    coord_np = init_coord.reshape([-1,3])
    from_coord_np_to_crdchg(coord_np = coord_np,
                            nmol=mols_args_dict["nmol"],
                            natoms=mols_args_dict["natoms"],
                            mol_list=mols_args_dict["mol_list"],
                            dic_mol_natom=mols_args_dict["dic_mol_natom"],
                            dic_mol_element=mols_args_dict["dic_mol_element"],
                            dic_mol_chg=mols_args_dict["dic_mol_chg"],
                            single_mol_charge_dic=mols_args_dict["single_mol_charge_dic"],
                            single_mol_spin_dic=mols_args_dict["single_mol_spin_dic"],
                            dic_mol_index=mols_args_dict["dic_mol_index"])
    # create the calculation directories. for both lower and higher state.
    if "SP_calculation_template_0" in mols_args_dict:
        if not os.path.exists('State_0'):
            os.system(
                "cp -r {SP_calculation_template_0} State_0".format(SP_calculation_template_0=mols_args_dict["SP_calculation_template_0"]))
        else:
            pass
            #os.system('rm -r State_0')
            #os.system(
            #    "cp -r {SP_calculation_template_0} State_0".format(SP_calculation_template_0=mols_args_dict["SP_calculation_template_0"]))
    else:
        print('You need to specific the value "SP_calculation_template_0" in the input dictionary.')
        exit(1)
    if "SP_calculation_template_1" in mols_args_dict:
        if not os.path.exists('State_1'):
            os.system(
                "cp -r {SP_calculation_template_1} State_1".format(SP_calculation_template_1=mols_args_dict["SP_calculation_template_1"]))
        else:
            pass
            #os.system('rm -r State_1')
            #os.system(
            #    "cp -r {SP_calculation_template_1} State_1".format(SP_calculation_template_1=mols_args_dict["SP_calculation_template_1"]))
    else:
        print('You need to specific the value "SP_calculation_template_1" in the input dictionary.')
        exit(1)


    mlebf_cmd_0 = mols_args_dict['SP_calculation_cmd_0']
    mlebf_cmd_1 = mols_args_dict['SP_calculation_cmd_1']

def main_cal_NAC(state='Singlet',State_0_dir='State_0',State_1_dir='State_1'):
    ene_0 = np.loadtxt('./{State_0_dir}/tot_ene.ene'.format(State_0_dir=State_0_dir))
    ene_1 = np.loadtxt('./{State_1_dir}/tot_ene.ene'.format(State_1_dir=State_1_dir))
    grad_0 = np.loadtxt('./{State_0_dir}/tot_grad.grad'.format(State_0_dir=State_0_dir))
    grad_1 = np.loadtxt('./{State_1_dir}/tot_grad.grad'.format(State_1_dir=State_1_dir))
    hess_0 = np.loadtxt('./{State_0_dir}/tot_hess.hess'.format(State_0_dir=State_0_dir))
    hess_1 = np.loadtxt('./{State_1_dir}/tot_hess.hess'.format(State_1_dir=State_1_dir))
    if state.lower() == 'singlet':
        nac = calculate_nac_SS(ene_0,ene_1,grad_0,grad_1,hess_0,hess_1)
        return nac
    elif state.lower() == 'triplet':
        nac = calculate_nac_TT(ene_0,ene_1,grad_0,grad_1,hess_0,hess_1)
        return nac


if __name__ == '__main__':
    state = sys.argv[0]
    State_0_dir = sys.argv[1]
    State_1_dir = sys.argv[2]
    nac_file_name = sys.argv[3]
    nac = main_cal_NAC(state,State_0_dir,State_1_dir)
    np.savetxt(nac_file_name,X=nac,fmt='%20.10f')
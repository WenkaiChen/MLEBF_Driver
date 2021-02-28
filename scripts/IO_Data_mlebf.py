#! /usr/bin/env python3
#-*- coding : utf-8-*-
# coding:unicode_escape

# This file is consisting of functions for reading the output data from external
# Quantum Chemistry Packages, such as Gaussian.

import math
import numpy as np
import os
import sys
from math import sqrt, exp, log

from ase.data import chemical_symbols, atomic_numbers
from ase.units import  Hartree
from ase.calculators.calculator import Calculator, FileIOCalculator, Parameters, kpts2mp,ReadError
from datetime import datetime
from ase import Atoms
from ase.optimize import LBFGS
from ase.optimize import BFGS
from ase.constraints import FixAtoms
from ase.constraints import FixBondLength
from ase.constraints import FixInternals  #This class allows to fix an arbitrary number of bond lengths, angles and dihedral angles.
from ase.io import read, write, iread
from math import pi

atmcnum2ele_dic = {
    1:  'H'  ,
    2:  'He' ,
    3:  'Li' ,
    4:  'Be' ,
    5:  'B'  ,
    6:  'C'  ,
    7:  'N'  ,
    8:  'O'  ,
    9:  'F'  ,
    10:  'Ne' ,
    11:  'Na' ,
    12:  'Mg' ,
    13:  'Al' ,
    14:  'Si' ,
    15:  'P'  ,
    16:  'S'  ,
    17:  'Cl' ,
    18:  'Ar' ,
    19:  'K'  ,
    20:  'Ca' ,
    21:  'Sc' ,
    22:  'Ti' ,
    23:  'V'  ,
    24:  'Cr' ,
    25:  'Mn' ,
    26:  'Fe' ,
    27:  'Co' ,
    28:  'Ni' ,
    29:  'Cu' ,
    30:  'Zn' ,
    31:  'Ga' ,
    32:  'Ge' ,
    33:  'As' ,
    34:  'Se' ,
    35:  'Br' ,
    36:  'Kr' ,
    37:  'Rb' ,
    38:  'Sr' ,
    39:  'Y'  ,
    40:  'Zr' ,
    41:  'Nb' ,
    42:  'Mo' ,
    43:  'Tc' ,
    44:  'Ru' ,
    45:  'Rh' ,
    46:  'Pd' ,
    47:  'Ag' ,
    48:  'Cd' ,
    49:  'In' ,
    50:  'Sn' ,
    51:  'Sb' ,
    52:  'Te' ,
    53:  'I'  ,
    54:  'Xe' ,
    55:  'Cs' ,
    56:  'Ba' ,
    57:  'La' ,
    58:  'Ce' ,
    59:  'Pr' ,
    60:  'Nd' ,
    61:  'Pm' ,
    62:  'Sm' ,
    63:  'Eu' ,
    64:  'Gd' ,
    65:  'Tb' ,
    66:  'Dy' ,
    67:  'Ho' ,
    68:  'Er' ,
    69:  'Tm' ,
    70:  'Yb' ,
    71:  'Lu' ,
    72:  'Hf' ,
    73:  'Ta' ,
    74:  'W'  ,
    75:  'Re' ,
    76:  'Os' ,
    77:  'Ir' ,
    78:  'Pt' ,
    79:  'Au' ,
    80:  'Hg' ,
    81:  'Tl' ,
    82:  'Pb' ,
    83:  'Bi' ,
    84:  'Po' ,
    85:  'At' ,
    86:  'Rn' ,
    87:  'Fr' ,
    88:  'Ra' ,
    89:  'Ac' ,
    90:  'Th' ,
    91:  'Pa' ,
    92:  'U'  ,
    93:  'Np' ,
    94:  'Pu' ,
    95:  'Am' ,
    96:  'Cm' ,
    97:  'Bk' ,
    98:  'Cf' ,
    99:  'Es' ,
    100:  'Fm' ,
    101:  'Md' ,
    102:  'No' ,
    103:  'Lr' ,
    104:  'Rf' ,
    105:  'Db' ,
    106:  'Sg' ,
    107:  'Bh' ,
    108:  'Hs' ,
    109:  'Mt' ,
    110:  'Ds' ,
    111:  'Rg' ,
    112:  'Cn' ,
    114:  'Uuq',
    116:  'Uuh'}


def read_hessian_gaussian16_fchk(fchk_name):
    # get the number of atoms of the system and number of hessians
    os.system("grep 'Number of atoms' {fchk_name} | awk '{{ print $5 }}' > _tmp_natom_system".format(
        fchk_name=fchk_name))
    natom = np.loadtxt("_tmp_natom_system").astype(np.int32)
    os.system("rm _tmp_natom_system")
    nhess_fchk = (natom*3+1) * natom*3 / 2
    nhess_fchk_lines = math.ceil(nhess_fchk / 5)

    # read hessian
    os.system(
        "grep -A{nhess_fchk_lines} 'Cartesian Force Constants' {fchk_name}  | "
        "tail -n {nhess_fchk_lines} | xargs > _tmp_hess_origindata".format(
            fchk_name=fchk_name, nhess_fchk_lines=nhess_fchk_lines))
    hess_origin = np.loadtxt("_tmp_hess_origindata").reshape(-1)

    # reshape the hessian from a triangle vector to a matrix
    hess_full = generate_full_mat_from_lowtriangle_vector(hess_origin, n_dim = natom*3)
    os.system("rm _tmp_hess_origindata")
    return hess_full

def read_hessian_openmolcas_hdf5(hdf5_name):
    import h5py
    with h5py.File(hdf5_name, "r") as f:
        hess_origin = f["HESSIAN"][()]
        natom = len(f["CENTER_MASSES"][()])
    hess_full = generate_full_mat_from_lowtriangle_vector(hess_origin,natom*3)
    return hess_full

def read_dipole_moment_gaussian_fchk(fchk_name):
    os.system('grep -A1 "Dipole Moment" {fchk_name} | tail -n 1 > _tmp_dipole_moment'.format(fchk_name=fchk_name))
    dipole_moment = np.loadtxt('_tmp_dipole_moment')
    os.system('rm _tmp_dipole_moment')
    return dipole_moment.reshape(1,-1)

def read_dipole_derivatives_gaussian_fchk(fchk_name):
    # get the number of atoms of the system and number of hessians
    os.system("grep 'Number of atoms' {fchk_name} | awk '{{ print $5 }}' > _tmp_natom_system".format(
        fchk_name=fchk_name))
    natom = np.loadtxt("_tmp_natom_system").astype(np.int32)
    os.system("rm _tmp_natom_system")
    ndp_der_fchk = natom * 9
    ndp_der_fchk_lines = math.ceil(ndp_der_fchk / 5)
    # read dipole moment derivatives
    os.system(
        "grep -A{ndp_der_fchk_lines} 'Dipole Derivatives' {fchk_name}  | "
        "tail -n {ndp_der_fchk_lines} | xargs > _tmp_dp_der_origindata".format(
            fchk_name=fchk_name, ndp_der_fchk_lines=ndp_der_fchk_lines))
    dp_der = np.loadtxt('_tmp_dp_der_origindata').reshape(-1,9)
    return dp_der

def generate_full_mat_from_uptriangle_vector(vec_origin, n_dim):
    """
    This function is used for generating a full matrix from a vector containing the element of the triangle matrix
    :param vec_origin: A 1D np.array ( or a numpy vector ). the element of the triangle matrix.
    :param n_dim: A int. the dimension of the matrix
    :return: A 2D np.array ( or a matrix). the full matrix.
    """
    mat_uptri = np.zeros((n_dim,n_dim))
    count = 0
    for i in range(n_dim):
        for j in range(i,n_dim):
            mat_uptri[i][j] = vec_origin[count]
            count += 1
    mat_full = mat_uptri + mat_uptri.T
    np.fill_diagonal(mat_full,np.diag(mat_uptri))
    return mat_full

def generate_full_mat_from_lowtriangle_vector(vec_origin, n_dim):
    """
    This function is used for generating a full matrix from a vector containing the element of the triangle matrix
    :param vec_origin: A 1D np.array ( or a numpy vector ). the element of the triangle matrix.
    :param n_dim: A int. the dimension of the matrix
    :return: A 2D np.array ( or a matrix). the full matrix.
    """
    mat_lowtri = np.zeros((n_dim,n_dim))
    count = 0
    for i in range(n_dim):
        for j in range(i+1):
            mat_lowtri[i][j] = vec_origin[count]
            count += 1
    mat_full = mat_lowtri + mat_lowtri.T
    np.fill_diagonal(mat_full,np.diag(mat_lowtri))
    return mat_full

def read_coord_from_xyz_file(xyzfile_name):
    """
    read the coordinates from the xyz file
    :param xyzfile_name: file name of the xyzfile
    :return: A dictionary: { "element": A list, "coordinates": A 2D numpy.array }
    """
    with open(xyzfile_name,'r') as f:
        lines = f.readlines()
    natom = int(lines[0].split()[0])
    ele_list = []
    coord_list = []
    for i in range(2,natom+2):
        line = lines[i]
        ele_list.append(line.split()[0])
        coord_list.append([float(line.split()[1]),float(line.split()[2]),float(line.split()[3])])
    coord_np = np.array(coord_list)
    return {"element": ele_list,
            "coordinates":coord_np}

def read_coord_from_crdchg_file(crdchg_name):
    with open(crdchg_name, 'r') as f:
        lines = f.readlines()
    nmol = int(lines[0].split()[0])
    natom = int(lines[0].split()[1])
    ele_list = []
    coord_list = []
    for i in range(3,natom+nmol+2):
        line = lines[i]
        if len(line.split()) == 5:
            ele_list.append(line.split()[0])
            coord_list.append([float(line.split()[1]),float(line.split()[2]),float(line.split()[3])])
    coord_np = np.array(coord_list)
    return {"element": ele_list,
            "coordinates":coord_np}

def write_coord_np_to_xyz_file(coord_np, ele_list, xyzfile_name):
    """

    :param coord_np: A numpy array with shape (natom,3)
    :param ele_list: A list.
    :param xyzfile_name: A string.
    :return:
    """
    natom = coord_np.shape[0]
    with open(xyzfile_name, 'w') as f:
        f.write("{natom}\n\n".format(natom=natom))
        for j in range(natom):
            f.write("{ele:4s}{cx:16.8f}{cy:16.8f}{cz:16.8f}\n".format(ele=ele_list[j],
                                                                      cx=coord_np[j][0],
                                                                      cy=coord_np[j][1],
                                                                      cz=coord_np[j][2]))

def from_raw_to_xyzfile(filename, ele_list, xyz_prefix='coord', nmol=None):
    """
    convert DeepMD raw file format to xyz files.
    :param filename: A string. the filename of the coord.raw file
    :param ele_list: A list. The element of the molecular. such as ['O', 'H', 'H']
    :param nmol: A int. or None. number of geoms need to write.
    :return: nothing. write the coordinate xyz files.
    """
    coord_np = np.loadtxt(filename)
    if not nmol:
        nmol = coord_np.shape[0]
    natom_per_mol = int(coord_np.shape[1]/3)
    for i in range(1,nmol+1):
        print(i)
        with open(xyz_prefix+"_{}.xyz".format(i), 'w') as f:
            f.write("{natom}\n\n".format(natom=natom_per_mol))
            for j in range(natom_per_mol):
                f.write("{ele:4s}{cx:16.8f}{cy:16.8f}{cz:16.8f}\n".format(ele=ele_list[j],
                                                                          cx=coord_np[i][j*3],
                                                                          cy=coord_np[i][j*3+1],
                                                                          cz=coord_np[i][j*3+2]))

def generate_EANN_SP_1_from_xyz_file(xyzfile_name):
    coord_np = read_coord_from_xyz_file(xyzfile_name=xyzfile_name)["coordinates"]
    np.savetxt('1', X=coord_np, fmt="%20.10f", header="point=       1 ")

def from_coord_np_to_crdchg(coord_np,nmol, natoms, mol_list, dic_mol_natom,
                            dic_mol_element,dic_mol_chg,
                            single_mol_charge_dic,single_mol_spin_dic,
                            dic_mol_index,
                            pbc=False, lat_parm=None):
    current_mol = 0
    mol_num = []

    mol_num = [0]
    for _,i in enumerate(mol_list):
        if _ == 0:
            mol_num.append(dic_mol_natom[i])
        else:
            mol_num.append(dic_mol_natom[i]+mol_num[_])
    #print(tot_mol)
    resnum_old = 0
    for current_t_mol in range(1):
        if not pbc:
            snap_shot = "{nmol}  {natoms}\n-----\n".format(nmol=nmol,natoms=natoms)
        else:
            snap_shot = "{nmol}  {natoms}  {lat_parm}\n-----\n".format(nmol=nmol,natoms=natoms,lat_parm=lat_parm)
        current_mol = 0
        for index,line in enumerate(coord_np[current_t_mol*natoms:(current_t_mol+1)*natoms,:]):
            while True:
                if index < mol_num[current_mol+1]: # this is the correct molecular.
                    atom_index = index - mol_num[current_mol]
                    break
                else:
                    current_mol += 1
            mol_name = mol_list[current_mol]
            if atom_index == 0:
                snap_shot += '{single_mol_charge} {single_mol_spin} {index}\n'.format(single_mol_charge=single_mol_charge_dic[mol_name],
                                                                                      single_mol_spin=single_mol_spin_dic[mol_name],
                                                                                      index=dic_mol_index[mol_name])

            snap_shot += "{ele:4s}{x:15.8f}{y:15.8f}{z:15.8f}{chg:11.5f}\n".format(ele=dic_mol_element[mol_name][atom_index],
                                                                                   x=float(line[0]),
                                                                                   y=float(line[1]),
                                                                                   z=float(line[2]),
                                                                                   chg=dic_mol_chg[mol_name][atom_index])
        with open('cluster_{}.crdchg'.format(current_t_mol), 'w') as f:
            f.write(snap_shot)

def from_coord_np_to_crdchg_freeze(coord_np,nmol, natoms, mol_list, dic_mol_natom,
                                   dic_mol_element,dic_mol_chg,
                                   single_mol_charge_dic,single_mol_spin_dic,
                                   dic_mol_index, active_freeze_index_list,
                                   pbc=False, lat_parm=None):
    current_mol = 0
    mol_num = []

    mol_num = [0]
    for _,i in enumerate(mol_list):
        if _ == 0:
            mol_num.append(dic_mol_natom[i])
        else:
            mol_num.append(dic_mol_natom[i]+mol_num[_])
    #print(tot_mol)
    resnum_old = 0
    for current_t_mol in range(1):
        if not pbc:
            snap_shot = "{nmol}  {natoms}\n-----\n".format(nmol=nmol,natoms=natoms)
        else:
            snap_shot = "{nmol}  {natoms}  {lat_parm}\n-----\n".format(nmol=nmol,natoms=natoms,lat_parm=lat_parm)
        current_mol = 0
        for index,line in enumerate(coord_np[current_t_mol*natoms:(current_t_mol+1)*natoms,:]):
            while True:
                if index < mol_num[current_mol+1]: # this is the correct molecular.
                    atom_index = index - mol_num[current_mol]
                    break
                else:
                    current_mol += 1
            mol_name = mol_list[current_mol]
            if atom_index == 0:
                snap_shot += '{single_mol_charge} {single_mol_spin} {index} {af}\n'.format(single_mol_charge=single_mol_charge_dic[mol_name],
                                                                                           single_mol_spin=single_mol_spin_dic[mol_name],
                                                                                           index=dic_mol_index[mol_name],
                                                                                           af=active_freeze_index_list[current_mol])

            snap_shot += "{ele:4s}{x:15.8f}{y:15.8f}{z:15.8f}{chg:11.5f}\n".format(ele=dic_mol_element[mol_name][atom_index],
                                                                                   x=float(line[0]),
                                                                                   y=float(line[1]),
                                                                                   z=float(line[2]),
                                                                                   chg=dic_mol_chg[mol_name][atom_index])
        with open('cluster_{}.crdchg'.format(current_t_mol), 'w') as f:
            f.write(snap_shot)

def from_g16_external_to_xyz(g16_name,xyzfile_name):
    """
    the function to convert g16 external input file to xyz file.
    :param g16_name: A string.
    :param xyzfile_name: A string.
    """
    with open(g16_name, 'r') as f:
        lines = f.readlines()
    line0 = lines[0]
    natoms= int(line0.split()[0])
    deriv = int(line0.split()[1])
    charge= int(line0.split()[2])
    spin  = int(line0.split()[3])

    coords = np.zeros((natoms,3))
    ele_list = []

    for i, line in enumerate(lines[1:1+natoms]):

        tokens = line.split()
        a = atmcnum2ele_dic[int(tokens[0])]

        c = np.array([float(tokens[1]), float(tokens[2]),float(tokens[3])])*0.52917721092

        coords[i,0] = c[0]
        coords[i,1] = c[1]
        coords[i,2] = c[2]

        ele_list.append(a)
    write_coord_np_to_xyz_file(coord_np=coords,ele_list=ele_list,xyzfile_name=xyzfile_name)

def from_ene_force_np_to_g16_output(output_filename,energy, dipole_moment=None, force=None):
    if dipole_moment is None:
        dipole_moment = np.array([0.0,0.0,0.0])
    with open(output_filename,'w') as f:
        f.write('{ene:20.12f}{dpx:20.12f}{dpy:20.12f}{dpz:20.12f}\n'.format(ene=energy,
                                                                          dpx=dipole_moment[0],
                                                                          dpy=dipole_moment[1],
                                                                          dpz=dipole_moment[2]))
        if force is not None:
            for i in range(force.shape[0]):
                f.write('{x:20.12f}{y:20.12f}{z:20.12f}\n'.format(x=force[i][0],
                                                                  y=force[i][1],
                                                                  z=force[i][2]))



if __name__ == "__main__":
    filename = sys.argv[1]
    data_name = sys.argv[2]
    file_suffix = filename.split('.')[-1]
    if ((data_name.lower() == 'hess') or (data_name.lower() == 'hessian')):
        if file_suffix.lower() == 'fchk': # Gaussian
            hess = read_hessian_gaussian16_fchk(filename)
        if file_suffix.lower() == 'h5':
            hess = read_hessian_openmolcas_hdf5(filename)
        save_hess_name = filename.split('.')[0] + '.hess'
        np.savetxt(save_hess_name, X=hess, fmt="%18.9e")
    elif (data_name.lower() == 'eann_1'):
        if file_suffix.lower == 'xyz':
            generate_EANN_SP_1_from_xyz_file(filename)
    elif ((data_name.lower() == 'dp') or (data_name.lower() == 'dipole_moment')):
        if file_suffix.lower() == 'fchk': # Gaussian
            dipole_moment = read_dipole_moment_gaussian_fchk(filename)
        save_dp_name = filename.split('.')[0] + '.dp'
        np.savetxt(save_dp_name, X=dipole_moment, fmt='%18.9e')
    elif ((data_name.lower() == 'dp_der') or (data_name.lower() == 'dipole_moment_derivatives')):
        if file_suffix.lower() == 'fchk': # Gaussian
            dp_der = read_dipole_derivatives_gaussian_fchk(filename)
        save_dp_name = filename.split('.')[0] + '.dp_der'
        np.savetxt(save_dp_name, X=dp_der, fmt='%18.9e')
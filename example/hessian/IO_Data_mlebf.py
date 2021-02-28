#! /usr/bin/env python3
#-*- coding : utf-8-*-
# coding:unicode_escape

# This file is consisting of functions for reading the output data from external
# Quantum Chemistry Packages, such as Gaussian.

import math
import numpy as np
import os
import sys

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
        ele_list.append(line.split()[1])
        coord_list.append([float(line.split()[1]),float(line.split()[2]),float(line.split()[3])])
    coord_np = np.array(coord_list)
    return {"element": ele_list,
            "coordinates":coord_np}


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

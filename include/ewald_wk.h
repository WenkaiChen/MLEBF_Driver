#ifndef MBE_DRIVER_EWALD_WK_H
#define MBE_DRIVER_EWALD_WK_H
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <complex>
#include <cstring>
#include <string.h>
#include <iomanip>
#include "element.h"
#include "pbc_tools.h"
#include "classes.h"
#define M_PI 3.14159265358979323846 // only useful in VC++

void get_crd_chg_mol(const Pbc_cell_smols & uc,
                     vector<vector<double > > & coord,
                     vector<double > & charge,
                     vector<vector<int > > & sc_index,
                     vector<int > & mol_index,
                     vector<double > & masses);

void get_crd_chg_mol_uc(const Pbc_cell_smols & uc,
                        vector<vector<double > > & coord,
                        vector<double > & charge,
                        vector<int > & mol_index,
                        vector<double > & masses);

void put_in_cell(vector<vector<double > > & coord,
                 const vector<vector<double > > & lv );

void get_PBC(double &RX, double &RY, double &RZ,
             const double & half_L, const double & half_B, const double & half_H,
             const double & L, const double & B, const double & H);

void ewald_ene_grad(const Pbc_cell_smols & unit_cell,
                    const Pbc_cell_smols & super_cell,
                    const vector<vector<double > > & lattice_vec,
                    const vector<vector<double > > & rec_lac_vec,
                    double & ewald_ene,
                    vector<vector<double > > & ewald_grad);

void uc_ee_in_sc_ene_grad(const Pbc_cell_smols & unit_cell,
                          const Pbc_cell_smols & super_cell,
                          const vector<vector<double > > & lattice_vec,
                          double & sc_ee_ene,
                          vector<vector<double > > & sc_ee_grad); 

void write_ewald_ene_grad(
        const double & ewald_ene,
        const vector<vector<double > > & ewald_grad);
#endif //MBE_DRIVER_EWALD_WK_H

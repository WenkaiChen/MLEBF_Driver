#ifndef MBE_DRIVER_BUILD_MONOMOER_H
#define MBE_DRIVER_BUILD_MONOMOER_H
#include <map>
#include <iomanip>
#include <sstream>
#include "../include/classes.h"
#include "../include/MathUtilities.h"


void
build_monomers(vector<Monomer > & monomers,
               const vector<SingleMol > & smols,
               const vector<vector<double > > & dis_mat,
               const double & rcut,
               map<string, int> & map_monomer);

void
build_monomers(vector<Monomer > & monomers,
               Pbc_cell_smols & uc_smols,
               Pbc_cell_smols & sc_smols,
               const int & stt_index,
               vector<int > & num_uc,
               vector<vector<double > > & dis_mat,
               const double & rcut,
               map<string, int> & map_monomer); 

void
build_dimers(vector<Dimer > & dimers,
             const vector<SingleMol > & smols,
             const vector<vector<double > > & dis_mat,
             const double & rcut_qm,
             const double & rcut_sr,
             map<string, int> & map_dimer);

void
build_dimers(vector<Dimer > & dimers,
             Pbc_cell_smols & uc_smols,
             Pbc_cell_smols & sc_smols,
             const int & stt_index,
             vector<int > & num_uc,
             const vector<vector<double > > & dis_mat,
             const double & rcut_qm,
             const double & rcut_sr,
             map<string, int> & map_dimer);

void
build_trimers(vector<Trimer > & trimers,
              const vector<SingleMol > & smols,
              const vector<vector<double > > & dis_mat,
              const double & rcut_qm,
              const double & rcut_sr,
              map<string, int> & map_trimer,
              bool reduce_order,
              vector<int> ro_index);

void
build_trimers(vector<Trimer > & trimers,
              Pbc_cell_smols & uc_smols,
              Pbc_cell_smols & sc_smols,
              const int & stt_index,
              vector<int > & num_uc,
              const vector<vector<double > > & dis_mat,
              const double & rcut_qm,
              const double & rcut_sr,
              map<string, int> & map_trimer,
              bool reduce_order,
              vector<int > ro_index);

#endif //MBE_DRIVER_BUILD_MONOMOER_H

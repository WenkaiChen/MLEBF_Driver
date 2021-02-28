#ifndef MBE_DRIVER_CAL_DIP_MNT_DER_H
#define MBE_DRIVER_CAL_DIP_MNT_DER_H
#include "classes.h"
#include <vector>
#include <sstream>
#include <iostream>
#include <fstream>
#include "vector_helper.h"
#include "read_crdchg.h" 

void
cal_dipole_moment_der(
                  vector<vector<double > > & dip_mnt_der,
                  vector<Monomer > & monomers,
                  vector<Dimer > & dimers,
                  vector<Trimer > & trimers,
                  const double & rcut_qm,
                  const double & rcut_sr,
                  const double & rcut_lr,
                  const map<string, int> & map_monomer,
                  const map<string, int> & map_dimer,
                  const map<string, int> & map_trimer,
                  const bool & directlyDelta,
                  const vector<int> & dd_index,
                  const bool & reduce_order,
                  const vector<int> & ro_index,
                  const bool & partial_nobgc,
                  const bool & mono_bgc,
                  const bool & cal_grad,
                  const int & order=2,
                  const bool & use_state=false,
                  const int & tot_state=0,
                  const int & current_state=0,
                  const bool & pbc=false,
                  const int & stt_index=-1);

void
read_dp_der_mono(vector<Monomer > & monomers);

void
read_dp_der_dimer(vector<Dimer > & dimers,
              vector<Monomer > & monomers,
              const double & rcut_qm,
              const double & rcut_sr,
              const double & rcut_lr,
              const map<string, int> & map_monomer,
              const map<string, int> & map_dimer,
              const bool & directlyDelta,
              const vector<int> & dd_index,
              const bool & reduce_order,
              const vector<int> &  ro_index,
              const bool & partial_nobgc,
              const bool & mono_bgc,
              const bool & pbc);

void
read_dp_der_trimer(vector<Trimer > & trimers,
               const vector<Dimer > & dimers,
               const vector<Monomer > & monomers,
               const double & rcut_qm,
               const double & rcut_sr,
               const double & rcut_lr,
               const map<string, int> & map_monomer,
               const map<string, int> & map_dimer,
               const map<string, int> & map_trimer,
               const bool & directlyDelta,
               const vector<int> & dd_index,
               const bool & pbc);

vector<vector<double > > mbe_dipole_moment_der(
        const vector<Monomer > & monomers,
        vector<Dimer > & dimers,
        const double & rcut_qm,
        const double & rcut_sr,
        const double & rcut_lr,
        const int & stt_index);

vector<vector<double > > mbe_dipole_moment_der(
        const vector<Monomer > & monomers,
        vector<Dimer > & dimers,
        const vector<Trimer > & trimers,
        const double & rcut_qm,
        const double & rcut_sr,
        const double & rcut_lr,
        const bool & directlyDelta,
        const vector<int> & dd_index,
        const int & stt_index);

void
cal_mono_dp_der(vector<vector<double > > & dip_mnt_der,
            const vector<Monomer > & monomers,
            const map<int, int> & natom_set,
            const map<int, int> & natom_cum);

void
cal_dimer_dp_der(vector<vector<double > > & dip_mnt_der,
             vector<Dimer > & dimers,
             const double & rcut_qm,
             const double & rcut_sr,
             const double & rcut_lr,
             const map<int, int> & natom_set,
             const map<int, int> & natom_cum,
             const int & stt_index);

void
cal_trimer_dp_der(vector<vector<double > > & dip_mnt_der,
              const vector<Monomer > & monomers,
              const vector<Trimer > & trimers,
              const double & rcut_qm,
              const double & rcut_sr,
              const double & rcut_lr,
              const map<int, int> & natom_set,
              const map<int, int> & natom_cum,
              const bool & directlyDelta,
              const vector<int> & dd_index,
              const int & stt_index);

vector<vector<double > >
read_mono_dp_mnt_der_tmp(const stringstream & dp_der_name, const int & natom);

void
write_output_dp_der(const vector<vector<double > > & dip_mnt_der);

#endif //MBE_DRIVER_CAL_DIP_MNT_DER_H

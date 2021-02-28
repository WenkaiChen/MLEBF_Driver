#ifndef MBE_DRIVER_CAL_TRANSITION_DIP_MNT_H
#define MBE_DRIVER_CAL_TRANSITION_DIP_MNT_H
#include "classes.h"
#include <vector>
#include <sstream>
#include <iostream>
#include <fstream>
#include "vector_helper.h"
#include "read_crdchg.h" 

void
cal_trans_dip_mnt(double & trans_dip_mnt,
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
                  const vector<int> &  ro_index,
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
read_trans_dp_mono(vector<Monomer > & monomers);

void
read_trans_dp_dimer(vector<Dimer > & dimers,
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
read_trans_dp_trimer(vector<Trimer > & trimers,
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

double mbe_trans_dip_mnt(
        const vector<Monomer > & monomers,
        vector<Dimer > & dimers,
        const double & rcut_qm,
        const double & rcut_sr,
        const double & rcut_lr,
        const int & stt_index);

double mbe_trans_dip_mnt(
        const vector<Monomer > & monomers,
        vector<Dimer > & dimers,
        const vector<Trimer > & trimers,
        const double & rcut_qm,
        const double & rcut_sr,
        const double & rcut_lr,
        const bool & directlyDelta,
        const vector<int> & dd_index,
        const int & stt_index);

double
read_mono_trans_dpmnt_tmp(const stringstream & tdpname);

void
write_output_trans_dp(const double & trans_dip_mnt);
#endif //MBE_DRIVER_CAL_TRANSITION_DIP_MNT_H

#ifndef MBE_DRIVER_CAL_EANDG_H
#define MBE_DRIVER_CAL_EANDG_H

#include "classes.h"
#include "read_crdchg.h" // want to use function split()
#include "vector_helper.h" // want to use function is_element_in_vector
#include "set_freeze_mol.h"
#include <iostream>
#include <iomanip>
#include <sstream>
#include <map>
#include <cstdlib> // system to exec shell script

/* old version
void
cal_eandg(double & ene,
          vector<vector<double > > & grad,
          vector<Monomer > & monomers,
          vector<Dimer > & dimers,
          const double & rcut_qm,
          const double & rcut_sr,
          const double & rcut_lr);
*/

void
cal_eandg(double & ene,
          vector<vector<double > > & grad,
          vector<Monomer > & monomers,
          vector<Dimer > & dimers,
          vector<Trimer > & trimers,
          const bool & has_freeze_mol,
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
          const int & stt_index=-1); // stt_index = -1 means no-pbc.

void
read_eg_mono(vector<Monomer > & monomers);

void
read_eg_dimer(vector<Dimer > & dimers,
              const double & rcut_qm,
              const double & rcut_sr,
              const double & rcut_lr);

void
read_eg_dimer(vector<Dimer > & dimers,
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
read_eg_trimer(vector<Trimer > & trimers,
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

double mbe_energy(const vector<Monomer > & monomers,
                  const vector<Dimer > & dimers,
                  const double & rcut_qm,
                  const double & rcut_sr,
                  const double & rcut_lr,
                  const int & stt_index);

vector<vector<double > > mbe_gradient(
        const vector<Monomer > & monomers,
        vector<Dimer > & dimers,
        const double & rcut_qm,
        const double & rcut_sr,
        const double & rcut_lr,
        const int & stt_index);

double mbe_energy(const vector<Monomer > & monomers,
                  const vector<Dimer > & dimers,
                  const vector<Trimer > & trimers,
                  const double & rcut_qm,
                  const double & rcut_sr,
                  const double & rcut_lr,
                  const bool & directlyDelta,
                  const vector<int> & dd_index,
                  const int & stt_index);

vector<vector<double > > mbe_gradient(
        const vector<Monomer > & monomers,
        vector<Dimer > & dimers,
        const vector<Trimer > & trimers,
        const double & rcut_qm,
        const double & rcut_sr,
        const double & rcut_lr,
        const bool & directlyDelta,
        const vector<int> & dd_index,
        const int & stt_index);

double cal_dimer_esr(const Dimer & dimer); // useless now

double cal_dimer_elr(const Dimer & dimer); // useless now

double cal_dimer_emm(const Dimer & dimer, const double & prefactor);

void
cal_dimer_gmm(Dimer & dimer, const double & prefactor);

void
cal_mono_grad(vector<vector<double > > & grad,
        const vector<Monomer > & monomers,
        const map<int, int> & natom_set,
        const map<int, int> & natom_cum);

void
cal_dimer_grad(vector<vector<double > > & grad,
               vector<Dimer > & dimers,
               const double & rcut_qm,
               const double & rcut_sr,
               const double & rcut_lr,
               const map<int, int> & natom_set,
               const map<int, int> & natom_cum,
               const int & stt_index);

void
cal_trimer_grad(vector<vector<double > > & grad,
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
read_mono_grad_tmp(const stringstream & gradname, const int & natom);

double
read_mono_ene_tmp(const stringstream & enename);

vector<vector<double > >
read_mono_grad_from_dimer(const stringstream & gradname,
                          const int & natom,
                          const int & start_index,
                          const int & final_index);

void
read_only_e_mono(vector<Monomer > & monomers);

void
read_only_e_dimer(vector<Dimer > & dimers,
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
                  const bool & mono_bgc);

void
read_only_e_trimer(vector<Trimer > & trimers,
                   const vector<Dimer > & dimers,
                   const vector<Monomer > & monomers,
                   const double & rcut_qm,
                   const double & rcut_sr,
                   const double & rcut_lr,
                   const map<string, int> & map_monomer,
                   const map<string, int> & map_dimer,
                   const map<string, int> & map_trimer,
                   const bool & directlyDelta,
                   const vector<int> & dd_index);

vector<vector<double > >
read_dimer_grad_not_in_dimers(const int & index, const int & natom);

double
read_dimer_ene_not_in_dimers(const int & index);

vector<vector<double > >
read_monomer_grad_not_in_monomers(const int & index, const int & natom);

double
read_monomer_ene_not_in_monomers(const int & index);

#endif //MBE_DRIVER_CAL_EANDG_H

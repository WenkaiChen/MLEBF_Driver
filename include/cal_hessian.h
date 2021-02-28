#ifndef MBE_DRIVER_CAL_HESSIAN_H
#define MBE_DRIVER_CAL_HESSIAN_H

#include "classes.h"
#include <vector>
#include <sstream>
#include <iostream>
#include <fstream>
#include "vector_helper.h"
#include "read_crdchg.h" 
#include "set_freeze_mol.h"

void
cal_hessian(vector<vector<double > > & hess,
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
            const int & stt_index=-1);

void
read_hessian_mono(vector<Monomer > & monomers);

void
read_hessian_dimer(vector<Dimer > & dimers,
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
read_hessian_trimer(vector<Trimer > & trimers,
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

vector<vector<double > > mbe_hessian(
        const vector<Monomer > & monomers,
        vector<Dimer > & dimers,
        const double & rcut_qm,
        const double & rcut_sr,
        const double & rcut_lr,
        const int & stt_index);

vector<vector<double > > mbe_hessian(
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
cal_mono_hess(vector<vector<double > > & hess,
              const vector<Monomer > & monomers,
              const map<int, int> & natom_set,
              const map<int, int> & natom_cum);

void
cal_dimer_hess(vector<vector<double > > & hess,
               vector<Dimer > & dimers,
               const double & rcut_qm,
               const double & rcut_sr,
               const double & rcut_lr,
               const map<int, int> & natom_set,
               const map<int, int> & natom_cum,
               const int & stt_index);

void
cal_trimer_hess(vector<vector<double > > & hess,
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
read_mono_hess_from_dimer(const stringstream & hessname,
                          const int & nh, // dimension of hessian of monomer. (3*Natom)
                          const int & start_index,
                          const int & final_index);

vector<vector<double > >
read_mono_hess_tmp(const stringstream & hessname, const int & nh);

void
write_output_hessian(const vector<vector<double > > & hessian);
#endif //MBE_DRIVER_CAL_HESSIAN_H

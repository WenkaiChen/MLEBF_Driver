#ifndef MBE_DRIVER_QM_SP_DRIVER_H
#define MBE_DRIVER_QM_SP_DRIVER_H

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <unistd.h>
#include <cstring>
#include "element.h"
#include "classes.h"


void
write_monomer(const vector<Monomer > & monomers,
              const int & n_mono,
              const bool & reduce_order,
              const vector<int> & ro_index,
              const bool & partial_nobgc,
              const bool & mono_bgc);

void
write_dimer(const vector<Dimer > & dimers,
            const int & n_dimer,
            const double & rcut_qm,
            const bool & reduce_order,
            const vector<int> & ro_index,
            const bool & partial_nobgc,
            const bool & pbc);

void
write_trimer(const vector<Trimer > & trimers,
             const int & n_trimer,
             const double & rcut_qm,
             const bool & pbc,
             map<string, int> & map_monomer,
             map<string, int> & map_dimer,
             map<string, int> & map_trimer);

void
qm_sp_driver(const int & method,
             const vector<Monomer > & monomers,
             const vector<Dimer > & dimers,
             const vector<Trimer > & trimers,
             map<string, int> & map_monomer,
             map<string, int> & map_dimer,
             map<string, int> & map_trimer,
             const double & rcut_qm,
             const bool & reduce_order,
             const vector<int> & ro_index,
             const bool & partial_nobgc,
             const bool & mono_bgc,
             int order=2,
             bool pbc= false);

void
write_output(const double & tot_ene,
             vector<vector<double > > & tot_grad,
             const bool & cal_grad);

#endif //MBE_DRIVER_QM_SP_DRIVER_H

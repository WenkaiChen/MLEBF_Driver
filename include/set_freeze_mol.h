#ifndef MBE_DRIVER_SET_FREEZE_MOL_H
#define MBE_DRIVER_SET_FREEZE_MOL_H

#include <vector>
#include "classes.h"
#include "vector_helper.h"

void
set_freeze_grad(const vector<Monomer > & monomers,
                vector<vector<double > > & tot_grad);

void
set_freeze_hess(const vector<Monomer > & monomers,
                vector<vector<double > > & tot_hess);


#endif //MBE_DRIVER_SET_FREEZE_MOL_H

#ifndef MBE_DRIVER_PBC_TOOLS_H
#define MBE_DRIVER_PBC_TOOLS_H
#include <vector>
#include <cmath>
#include "../include/vector_helper.h"

using namespace std;

void
cell2cellpar(const vector<vector<double > > & lattice_vector,
             vector<double > & pbc_box_info);

void
car2frac(const vector<vector<double > > & coord_cartesian,
         const vector<vector<double > > & lattice_vector,
         vector<vector<double > > & coord_fraction);

void
frac2car(const vector<vector<double > > & coord_fraction,
         const vector<vector<double > > & lattice_vector,
         vector<vector<double > > & coord_cartesian);

vector<vector<double > >
reciprical_lattice_vectors(const vector<vector<double > > & lat_vec);

double
compute_volume_cell(const vector<vector<double > > & lattice_vector);

#endif //MBE_DRIVER_PBC_TOOLS_H

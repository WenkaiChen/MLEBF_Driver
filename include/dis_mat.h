#ifndef MBE_DRIVER_DIS_MAT_H
#define MBE_DRIVER_DIS_MAT_H

#include <string>
#include <vector>
#include "../include/classes.h"
#include "../include/MathUtilities.h"

vector<vector<double > >
compute_dis_mat(const vector<SingleMol > & smols);

void
compute_dis_mat(vector<vector<double > > & dis_mat,
                const vector<SingleMol > & smols,
                const string & fn_box);

#endif //MBE_DRIVER_DIS_MAT_H

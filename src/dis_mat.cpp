//
// Created by wenkai on 18-12-25.
//

#include "../include/dis_mat.h"

vector<vector<double > >
compute_dis_mat(const vector<SingleMol > & smols){
    int nmol = smols.size();
    vector<vector<double > > dis_mat(nmol, vector<double> (nmol,0.0));

    for (int i=0;i<nmol;i++){
        for (int j=i+1;j<nmol;j++){
            double x = smols[i].comass[0] - smols[j].comass[0];
            double y = smols[i].comass[1] - smols[j].comass[1];
            double z = smols[i].comass[2] - smols[j].comass[2];
            double dist = sqrt(x*x + y*y + z*z);
            dis_mat[i][j] = dist;
            dis_mat[j][i] = dist;
        }
    }
    return dis_mat;
}

void
compute_dis_mat(vector<vector<double > > & dis_mat,
        const vector<SingleMol > & smols,
        const string & fn_box){
    int nmol = smols.size();
    if (fn_box == " "){ // no PBC
        dis_mat = compute_dis_mat(smols);
    }
    else {

    }
}
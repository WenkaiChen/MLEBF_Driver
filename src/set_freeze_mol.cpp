#include "../include/set_freeze_mol.h"

void
set_freeze_grad(const vector<Monomer > & monomers,
                vector<vector<double > > & tot_grad) {
    map<int, int> natom_set;
    map<int, int> natom_cum; 
    int ntot_atom = 0;
    for (int i=0;i<monomers.size();i++){
        natom_set[i]=monomers[i].ele.size();
        ntot_atom += monomers[i].ele.size();
    }
    natom_cum[monomers[0].index] = 0;
    for (int i=1;i<monomers.size();i++){
        natom_cum[i] = natom_cum[i-1]+natom_set.at(i-1);
    }
    natom_cum[monomers.size()] = ntot_atom;
    for (int i=0; i<monomers.size();i++){
        if (monomers[i].is_frozen()){
            set_part_matrix_zero(tot_grad,natom_cum[i],natom_cum[i+1],0,3);
        }
    }
}

void
set_freeze_hess(const vector<Monomer > & monomers,
                vector<vector<double > > & tot_hess) {
    map<int, int> natom_set;
    map<int, int> natom_cum; 
    int ntot_atom = 0;
    for (int i=0;i<monomers.size();i++){
        natom_set[i]=monomers[i].ele.size();
        ntot_atom += monomers[i].ele.size();
    }
    natom_cum[monomers[0].index] = 0;
    for (int i=1;i<monomers.size();i++){
        natom_cum[i] = natom_cum[i-1]+natom_set.at(i-1);
    }
    natom_cum[monomers.size()] = ntot_atom;
    for (int i=0; i<monomers.size();i++){
        if (monomers[i].is_frozen()){
            set_part_matrix_zero(tot_hess,natom_cum[i]*3,natom_cum[i+1]*3,0,ntot_atom*3);
        }
        else {
            for (int j=0; j<monomers.size();j++){
                if (monomers[j].is_frozen()){
                    set_part_matrix_zero(tot_hess,natom_cum[i]*3,natom_cum[i+1]*3,natom_cum[j]*3,natom_cum[j+1]*3);
                }
            }
        }
    }
}

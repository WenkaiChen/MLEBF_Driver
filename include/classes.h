#ifndef MBE_DRIVER_STRUCTS_H
#define MBE_DRIVER_STRUCTS_H

#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include "element.h"
#include "vector_helper.h"
#include "pbc_tools.h"
#include <iomanip>

using namespace std;
class SingleMol{
public:
    int  species;
    int  tot_charge;
    int  tot_spin;
    vector<int >  sc_index;
    int  mol_index;
    vector<int > ele;
    vector<vector<double > > coord;
    vector<double > charge;
    vector<double > comass; 
    SingleMol(int s, int tc, int ts, vector<int> e, vector<vector<double > > crd, vector<double > chg, double rb)
            :species(s),tot_charge(tc),tot_spin(ts),ele(e),coord(crd),charge(chg){
        compute_centor_of_mass();
        if ( rb == 0.0 ) {frozen = false;}
        else{
            double r = sqrt(comass[0]*comass[0]+comass[1]*comass[1]+comass[2]*comass[2]);
            if (r<rb) frozen=false;
            else frozen=true;
        }
    }
    SingleMol(int s, int tc, int ts, vector<int> e, vector<vector<double > > crd, vector<double > chg, bool freeze)
            :species(s),tot_charge(tc),tot_spin(ts),ele(e),coord(crd),charge(chg){
        compute_centor_of_mass();
        frozen = freeze;
    }
    SingleMol(vector<int >  ci, int mi, int s, int tc, int ts, vector<int> e, vector<vector<double > > crd, vector<double > chg, double rb)
            :sc_index(ci), mol_index(mi), species(s),tot_charge(tc),tot_spin(ts),ele(e),coord(crd),charge(chg){
        compute_centor_of_mass();
        if ( rb == 0.0 ) {frozen = false;}
        else{
            double r = sqrt(comass[0]*comass[0]+comass[1]*comass[1]+comass[2]*comass[2]);
            if (r<rb) frozen=false;
            else frozen=true;
        }
    }
    SingleMol(vector<vector<double > > crd):coord(crd){}
    bool is_frozen(void) const { return frozen; }

private:
    bool frozen;
    void
    compute_centor_of_mass(void){
        comass.clear();
        comass.push_back(0.0);
        comass.push_back(0.0);
        comass.push_back(0.0);
        double tot_mass = 0;
        int natom = ele.size();
        for (int ii = 0; ii < natom; ii++){
            tot_mass += ele2mass.at(ele[ii]);
            comass[0] += ele2mass.at(ele[ii]) * coord[ii][0];
            comass[1] += ele2mass.at(ele[ii]) * coord[ii][1];
            comass[2] += ele2mass.at(ele[ii]) * coord[ii][2];
        }
        comass[0] = comass[0]/tot_mass;
        comass[1] = comass[1]/tot_mass;
        comass[2] = comass[2]/tot_mass;
    }

};

class Pbc_cell{
public:
    vector<double > pbc_box_info;
    vector<vector<double > > lattice_vector;
    vector<vector<double > > coord_cartesian;
    vector<vector<double > > coord_fraction;
    Pbc_cell(vector<double > pbi,vector<vector<double > > lv, vector<vector<double > > crd, bool cartesian){
        pbc_box_info = pbi;
        lattice_vector = lv;
        if (cartesian){
            coord_cartesian = crd;
            car2frac(coord_cartesian,lattice_vector,coord_fraction);
        }
        else{
            coord_fraction = crd;
            frac2car(coord_fraction,lattice_vector,coord_cartesian);
        }
    }
};

class Pbc_cell_smols{
public:
    vector<double > pbc_box_info;
    vector<vector<double > > lattice_vector;
    vector<SingleMol > smols;
    vector<double > dipole_momemt;
    vector<double > external_charges_chg;
    vector<vector<double > > external_charges_crd;
    Pbc_cell_smols(vector<double> & pbi, vector<vector<double > > & lv, vector<SingleMol > & sms, bool is_sc)
            :pbc_box_info(pbi),lattice_vector(lv),smols(sms){
        if (not is_sc) {
            int n_smols = smols.size();
            for (int i = 0; i < n_smols; i++) {
                smols[i].mol_index = i;
                smols[i].sc_index.clear();
                smols[i].sc_index.reserve(3);
                smols[i].sc_index.resize(3);
                for (int j = 0; j < 3; j++){
                    smols[i].sc_index[j] = 0;
                }
            }
        }
    }
    Pbc_cell_smols make_supercell(const vector<int> & num_uc){
        vector<SingleMol >  sc_smols;
        vector<double > sc_pbi;
        vector<vector<double > > sc_lv;
        vector<double >sc_ecchg;
        vector<vector<double > > sc_eccrd;
        sc_pbi.resize(6);
        sc_lv.resize(3);
        int nmol = smols.size();
        int stt;
        for (int a_mult=-num_uc[0];a_mult<num_uc[0]+1;a_mult++){
            for (int b_mult=-num_uc[1];b_mult<num_uc[1]+1;b_mult++){
                for (int c_mult=-num_uc[2];c_mult<num_uc[2]+1;c_mult++){
        //for (int a_mult=0;a_mult<num_uc[0]+1;a_mult++){
        //    for (int b_mult=0;b_mult<num_uc[1]+1;b_mult++){
        //        for (int c_mult=0;c_mult<num_uc[2]+1;c_mult++){

                    vector<double > trans_vec = double(a_mult)*lattice_vector[0] + double(b_mult)*lattice_vector[1]+
                                                double(c_mult)*lattice_vector[2];
                    vector<int > cell_index;
                    cell_index.resize(3);
                    cell_index[0] = a_mult; cell_index[1] = b_mult; cell_index[2] = c_mult;
                    for (int i_mol=0;i_mol<nmol;i_mol++){
                        vector<vector<double > > coord_new = smols[i_mol].coord + trans_vec;
                        SingleMol _singlemol(cell_index,i_mol,smols[i_mol].species,smols[i_mol].tot_charge,smols[i_mol].tot_spin,
                                             smols[i_mol].ele,coord_new,smols[i_mol].charge,0.0);
                        sc_smols.push_back(_singlemol);
                    }
                    if (a_mult == -num_uc[0]){
                        sc_ecchg.push_back(external_charges_chg[0]);
                        sc_eccrd.push_back(external_charges_crd[0] + trans_vec);
                    }
                    if (a_mult == num_uc[0]){
                        sc_ecchg.push_back(external_charges_chg[1]);
                        sc_eccrd.push_back(external_charges_crd[1] + trans_vec);
                    }
                    if (b_mult == -num_uc[1]){
                        sc_ecchg.push_back(external_charges_chg[2]);
                        sc_eccrd.push_back(external_charges_crd[2] + trans_vec);
                    }
                    if (b_mult == num_uc[1]){
                        sc_ecchg.push_back(external_charges_chg[3]);
                        sc_eccrd.push_back(external_charges_crd[3] + trans_vec);
                    }
                    if (c_mult == -num_uc[2]){
                        sc_ecchg.push_back(external_charges_chg[4]);
                        sc_eccrd.push_back(external_charges_crd[4] + trans_vec);
                    }
                    if (c_mult == num_uc[2]){
                        sc_ecchg.push_back(external_charges_chg[5]);
                        sc_eccrd.push_back(external_charges_crd[5] + trans_vec);
                    }
                }
            }
        }
        for (int i=0; i<3; i++){
            sc_pbi[i] = pbc_box_info[i] * (2*num_uc[i] + 1);
            sc_lv[i].resize(3);
            for (int j=0; j<3; j++){
                sc_lv[i][j] = lattice_vector[i][j] * (2*num_uc[i] + 1);
            }
        }
        sc_pbi[4] = pbc_box_info[4];
        sc_pbi[5] = pbc_box_info[5];
        sc_pbi[6] = pbc_box_info[6];
        Pbc_cell_smols sc(sc_pbi,sc_lv,sc_smols,true);
        sc.set_external_charges(sc_ecchg, sc_eccrd);
        return sc;
    }
    void calculate_dipole_moment(void){
        dipole_momemt.clear();
        dipole_momemt.resize(3);
        center_uc.resize(3);
        center_uc[0] = lattice_vector[0][0] * 0.5 + lattice_vector[1][0] * 0.5 + lattice_vector[2][0] * 0.5;
        center_uc[1] = lattice_vector[0][1] * 0.5 + lattice_vector[1][1] * 0.5 + lattice_vector[2][1] * 0.5;
        center_uc[2] = lattice_vector[0][2] * 0.5 + lattice_vector[1][2] * 0.5 + lattice_vector[2][2] * 0.5;
        int nmol = smols.size();
        for (int i=0; i<nmol;i++){
            int natom = smols[i].ele.size();
            for (int j=0;j<natom;j++) {
                vector<double> _coord_tmp = smols[i].coord[j] - center_uc;
                dipole_momemt = dipole_momemt + _coord_tmp * smols[i].charge[j];
            }
        }
    }
    void set_external_charges(void){
        calculate_dipole_moment();
        external_charges_chg.clear();
        external_charges_chg.reserve(6);
        external_charges_crd.clear();
        external_charges_crd.reserve(6);
        for (int i=0;i<3;i++){
            vector<double > coord1;
            vector<double > coord2;
            double charge;
            coord1 = center_uc;
            coord2 = center_uc;
            coord1[i] = coord1[i] - center_uc[i];
            coord2[i] = coord2[i] + center_uc[i];
            charge = dipole_momemt[i] / (coord2[i] - coord1[i]);
            external_charges_chg.push_back(charge);
            external_charges_crd.push_back(coord1);
            external_charges_chg.push_back(-charge);
            external_charges_crd.push_back(coord2);
        }
    }
    void set_external_charges(vector<double > & chg, vector<vector<double > > & crd){
        external_charges_crd = crd;
        external_charges_chg = chg;
    }
private:
    vector<double > center_uc;

};

class Monomer{
public:
    int index;
    int tot_charge;
    int tot_spin;
    int species;
    int cell_index;
    vector<int > sc_index;
    vector<int > ele;
    vector<vector<double > > coord;
    vector<double > bgc_chg;
    vector<vector<double > > bgc_crd;
    Monomer(vector<int> e, vector<vector<double > > c, vector<double > bchg, vector<vector<double > > bcrd)
            :ele(e),coord(c),bgc_chg(bchg),bgc_crd(bcrd){}
    Monomer(vector<int> e, vector<vector<double > > c, vector<double > bchg, vector<vector<double > > bcrd,
            int tc, int ts)
            :ele(e),coord(c),bgc_chg(bchg),bgc_crd(bcrd),tot_charge(tc),tot_spin(ts){}
    Monomer(vector<int> e, vector<vector<double > > c, vector<double > bchg, vector<vector<double > > bcrd,
            int tc, int ts, int sp)
            :ele(e),coord(c),bgc_chg(bchg),bgc_crd(bcrd),tot_charge(tc),tot_spin(ts),species(sp){}
    Monomer(vector<SingleMol> smols, int index_, vector<double > dis_vec, double rcut){
        form_monomer_from_smols(smols,index_,dis_vec,rcut);
    }
    Monomer(vector<SingleMol> smols, int index_, int ic, vector<double > dis_vec, double rcut){
        form_monomer_from_smols(smols,index_,ic,dis_vec,rcut);
    }
    Monomer(Pbc_cell_smols & uc, Pbc_cell_smols & sc, vector<int > & num_uc, int stt_index,
            int index_, int ic, vector<double > & dis_vec, double rcut){
        form_monomer_from_pbccells(uc, sc, num_uc, stt_index, index_, ic, dis_vec, rcut);
    };
    void get_energy(const double ene){ energy = ene; };
    const double & get_energy(void) const { return energy; };

    void get_gradient(const vector<vector<double > > grad){ gradient = grad; };
    const vector<vector<double > > & get_gradient(void) const { return gradient; };

    void get_hessian(const vector<vector<double > > hess){ hessian = hess; };
    const vector<vector<double > > & get_hessian(void) const { return hessian; };

    void get_dipole_moment(const vector<double > dp) { dipole_moment = dp; };
    const vector<double > & get_dipole_moment(void) const { return dipole_moment; };

    void get_dipole_moment_der(const vector<vector<double > > dp_der) { dipole_moment_derivative = dp_der; };
    const vector<vector<double > > & get_dipole_moment_der(void) const { return dipole_moment_derivative; };

    void get_trans_dip_mnt(const double tdp) { trans_dip_mnt = tdp; };
    const double & get_trans_dip_mnt(void) const { return trans_dip_mnt; };

    bool is_frozen(void) const { return frozen; };
private:
    void
    form_monomer_from_smols(const vector<SingleMol > & smols,
                            const int & index_,
                            const vector<double > & dis_vec,
                            const double & rcut){
        SingleMol mol = smols[index_];
        frozen = mol.is_frozen();
        index = index_;
        species = mol.species;
        ele = mol.ele;
        coord = mol.coord;
        tot_charge = mol.tot_charge;
        tot_spin = mol.tot_spin;
        bgc_chg.clear();
        bgc_crd.clear();
        for (int i=0;i<smols.size();i++){
            if (i!=index and dis_vec[i] < rcut){
                bgc_chg.insert(bgc_chg.end(),smols[i].charge.begin(),smols[i].charge.end());
                bgc_crd.insert(bgc_crd.end(),smols[i].coord.begin(),smols[i].coord.end());
            }
        }
    }
    void
    form_monomer_from_smols(const vector<SingleMol > & smols,
                            const int & index_,
                            const int & ic,
                            const vector<double > & dis_vec,
                            const double & rcut){
        SingleMol mol = smols[index_];
        frozen = mol.is_frozen();
        index = index_;
        cell_index = ic;
        species = mol.species;
        ele = mol.ele;
        coord = mol.coord;
        tot_charge = mol.tot_charge;
        tot_spin = mol.tot_spin;
        bgc_chg.clear();
        bgc_crd.clear();
        for (int i=0;i<smols.size();i++){
            if (i!=index and dis_vec[i] < rcut){
                bgc_chg.insert(bgc_chg.end(),smols[i].charge.begin(),smols[i].charge.end());
                bgc_crd.insert(bgc_crd.end(),smols[i].coord.begin(),smols[i].coord.end());
            }
        }
    }
    void
    form_monomer_from_pbccells(Pbc_cell_smols & uc, Pbc_cell_smols & sc,
            vector<int > & num_uc, int stt_index,
            int index_, int ic, vector<double > & dis_vec, double rcut_sr){
        SingleMol & mol = sc.smols[index_+stt_index];
        frozen = mol.is_frozen();
        index = index_ + stt_index;
        cell_index = ic;
        sc_index = mol.sc_index;
        species = mol.species;
        ele = mol.ele;
        coord = mol.coord;
        tot_charge = mol.tot_charge;
        tot_spin = mol.tot_spin;
        bgc_chg.clear();
        bgc_crd.clear();
        for(int i=0; i<sc.smols.size();i++){
            if ((i != index) ){
                bgc_chg.insert(bgc_chg.end(),sc.smols[i].charge.begin(),sc.smols[i].charge.end());
                bgc_crd.insert(bgc_crd.end(),sc.smols[i].coord.begin(),sc.smols[i].coord.end());
            }
        }
    }
    double energy;
    vector<vector<double > > gradient;
    bool frozen;
    vector<vector<double > > hessian;
    vector<double > dipole_moment;
    vector<vector<double > > dipole_moment_derivative;
    double trans_dip_mnt;
};

class Dimer{
public:
    double rij;
    int index1;
    int index2;
    int tot_charge;
    int tot_spin;
    int tc1;
    int tc2;
    int ts1;
    int ts2;
    int species1;
    int species2;
    int mol_index1_uc;
    int mol_index2_uc;
    vector<int > sc_index1;
    vector<int > sc_index2;
    vector<int > ele;
    vector<vector<double > > coord;
    vector<double > charge;
    vector<double > bgc_chg;
    vector<vector<double > > bgc_crd;
    vector<double > com1;
    vector<double > com2;
    Dimer(const vector<SingleMol> smols, int i1, int i2, const vector<vector<double > > dis_mat,
          double rcut_qm, double rcut_sr){
        form_dimer_from_smols(smols,i1,i2,dis_mat,rcut_qm,rcut_sr);
    }
    Dimer(const Pbc_cell_smols & uc, const Pbc_cell_smols & sc, const vector<vector<double > > & dis_mat,
            int i1, int i2, int stt_index, double rcut_sr){
        form_dimer_from_pbccells(uc, sc, dis_mat, i1, i2, stt_index, rcut_sr);
    }

    Monomer get_m1(void){
        Monomer m1(ele1,coord1,bgc_chg1,bgc_crd1,tc1,ts1,species1);
        return m1;
    }
    Monomer get_m2(void){
        Monomer m2(ele2,coord2,bgc_chg2,bgc_crd2,tc2,ts2,species2);
        return m2;
    }

    void get_energy(const double ene){ energy = ene; };
    void get_ene_m1(const double ene){ ene_m1 = ene; };
    void get_ene_m2(const double ene){ ene_m2 = ene; };
    const double & get_energy(void) const { return energy; };
    const double & get_ene_m1(void) const { return ene_m1; };
    const double & get_ene_m2(void) const { return ene_m2; };

    void get_gradient(const vector<vector<double > > grad){ gradient = grad; };
    void get_grad_m1 (const vector<vector<double > > grad){ grad_m1  = grad; };
    void get_grad_m2 (const vector<vector<double > > grad){ grad_m2  = grad; };
    const vector<vector<double > > & get_gradient(void) const { return gradient; };
    const vector<vector<double > > & get_grad_m1(void)  const { return grad_m1; };
    const vector<vector<double > > & get_grad_m2(void)  const { return grad_m2; };

    void get_hessian(const vector<vector<double > > hess){ hessian = hess; };
    void get_hess_m1(const vector<vector<double > > hess){ hess_m1 = hess; };
    void get_hess_m2(const vector<vector<double > > hess){ hess_m2 = hess; };
    const vector<vector<double > > & get_hessian(void) const { return hessian; };
    const vector<vector<double > > & get_hess_m1(void) const { return hess_m1; };
    const vector<vector<double > > & get_hess_m2(void) const { return hess_m2; };

    void get_dipole_moment(const vector<double > dp) { dipole_moment = dp; };
    void get_dp_m1(const vector<double > dp ) { dp_m1 = dp; };
    void get_dp_m2(const vector<double > dp ) { dp_m2 = dp; };
    const vector<double > & get_dipole_moment(void) const { return dipole_moment; };
    const vector<double > & get_dp_m1(void) const { return dp_m1; };
    const vector<double > & get_dp_m2(void) const { return dp_m2; };

    void get_dipole_moment_der(const vector<vector<double > > dp_der) { dipole_moment_derivative = dp_der; };
    void get_dp_der_m1(const vector<vector<double > > dp_der) { dp_der_m1 = dp_der; };
    void get_dp_der_m2(const vector<vector<double > > dp_der) { dp_der_m2 = dp_der; };
    const vector<vector<double > > & get_dipole_moment_der(void) const { return dipole_moment_derivative; };
    const vector<vector<double > > & get_dp_der_m1(void) const { return dp_der_m1; };
    const vector<vector<double > > & get_dp_der_m2(void) const { return dp_der_m2; };

    void get_trans_dip_mnt(const double tdp) { trans_dip_mnt = tdp; };
    const double & get_trans_dip_mnt(void) const { return trans_dip_mnt; };
    void get_trans_dp_m1(const double tdp) { trans_dp_m1 = tdp; };
    const double & get_trans_dp_m1(void) const { return trans_dp_m1; };
    void get_trans_dp_m2(const double tdp) { trans_dp_m2 = tdp; };
    const double & get_trans_dp_m2(void) const { return trans_dp_m2; };

    int get_natom_m1(){ return ele1.size(); };
    int get_natom_m2(){ return ele2.size(); };

    vector<int > ele1;
    vector<int > ele2;
    vector<vector<double > > coord1;
    vector<vector<double > > coord2;
    vector<double > charge1;
    vector<double > charge2;
    vector<double > bgc_chg1;
    vector<double > bgc_chg2;
    vector<vector<double > > bgc_crd1;
    vector<vector<double > > bgc_crd2;

    void set_dd_ene_grad(void) { 
        set_ene_m1_zero();
        set_ene_m2_zero();
        set_grad_m1_zero();
        set_grad_m2_zero();
    }
    void set_dd_hess(void) {
        set_hess_m1_zero();
        set_hess_m2_zero();
    }
    void set_dd_dp_mnt(void) {
        set_dp_mnt_m1_zero();
        set_dp_mnt_m2_zero();
    }
    void set_dd_dp_der(void) {
        set_dp_der_m1_zero();
        set_dp_der_m2_zero();
    }
    void set_dd_tdp(void) {
        set_tdp_m1_zero();
        set_tdp_m2_zero();
    }
    void set_partial_grad(const vector<vector<double > > & temp_grad,
                          const int start_index,
                          const int final_index){ 
        if (temp_grad.size() == (final_index-start_index)){
            for (int ii=0;ii<ele.size();ii++){
                if ((ii >= start_index) and (ii<final_index)) {
                    gradient[ii][0] = temp_grad[ii-start_index][0];
                    gradient[ii][1] = temp_grad[ii-start_index][1];
                    gradient[ii][2] = temp_grad[ii-start_index][2];
                }
            }
        }
        else{
            cout << "temp_grad.size() != (final_index-start_index), please check your code." << endl;
        }
    }
    void set_partial_hess(const vector<vector<double > > & temp_hess,
                          const int start_index,
                          const int final_index){
        if (temp_hess.size() == (final_index-start_index)){
            for (int ii=0;ii<ele.size()*3;ii++){
                if ((ii >= start_index) and (ii<final_index)) {
                    for (int jj=start_index;jj<final_index;jj++){
                        hessian[ii][jj] = temp_hess[ii-start_index][jj-start_index];
                    }
                }
            }
        }
        else{
            cout << "temp_hess.size() != (final_index-start_index), please check your code." << endl;
        }
    }

    void set_partial_dp_der(const vector<vector<double > > & temp_dp_der,
                            const int start_index,
                            const int final_index){
        if (temp_dp_der.size() == (final_index - start_index)){
            for (int ii=0;ii<ele.size();ii++){
                if ((ii >= start_index) and (ii<final_index)) {
                    for (int jj=0;jj<9;jj++){
                        dipole_moment_derivative[ii][jj] = temp_dp_der[ii-start_index][jj];
                    }
                }
            }
        }
        else{
            cout << "temp_dp_der.size() != (final_index-start_index), please check your code." << endl;
        }
    }

private:
    void form_dimer_from_smols(const vector<SingleMol> & smols,
                               const int & i1, const int & i2,
                               const vector<vector<double > > & dis_mat,
                               const double & rcut_qm, const double & rcut_sr){
        index1 = i1;
        index2 = i2;
        ele.clear();
        coord.clear();
        charge.clear();
        bgc_crd.clear();
        bgc_chg.clear();
        SingleMol mol1 = smols[i1];
        SingleMol mol2 = smols[i2];
        rij = dis_mat[i1][i2];
        
        com1 = mol1.comass;
        com2 = mol2.comass;
        species1=mol1.species;
        species2=mol2.species;
        
        ele.insert(ele.end(),mol1.ele.begin(),mol1.ele.end());
        ele.insert(ele.end(),mol2.ele.begin(),mol2.ele.end());
        ele1 = mol1.ele;
        ele2 = mol2.ele;
        
        coord.insert(coord.end(),mol1.coord.begin(),mol1.coord.end());
        coord.insert(coord.end(),mol2.coord.begin(),mol2.coord.end());
        coord1 = mol1.coord;
        coord2 = mol2.coord;
        
        charge.insert(charge.end(),mol1.charge.begin(),mol1.charge.end());
        charge.insert(charge.end(),mol2.charge.begin(),mol2.charge.end());
        charge1 = mol1.charge;
        charge2 = mol2.charge;
        
        tc1 = mol1.tot_charge;
        tc2 = mol2.tot_charge;
        ts1 = mol1.tot_spin;
        ts2 = mol2.tot_spin;
        tot_charge = tc1 + tc2;
        tot_spin = 2 * (tot_charge) + 1;
        if (ts2 != 1 or ts1 != 1){
            tot_spin = max(ts1,ts2);
        }
        
        for (int i=0;i<smols.size();i++){
            double r1 = dis_mat[i1][i];
            double r2 = dis_mat[i2][i];
            if ( (i!=i1 and i!=i2) and (r1<rcut_sr or r2<rcut_sr)){
                bgc_chg.insert(bgc_chg.end(),smols[i].charge.begin(),smols[i].charge.end());
                bgc_crd.insert(bgc_crd.end(),smols[i].coord.begin(),smols[i].coord.end());
            }
        }
        bgc_chg1.assign(bgc_chg.begin(),bgc_chg.end());
        bgc_chg1.insert(bgc_chg1.end(),mol2.charge.begin(),mol2.charge.end());
        bgc_chg2.assign(bgc_chg.begin(),bgc_chg.end());
        bgc_chg2.insert(bgc_chg2.end(),mol1.charge.begin(),mol1.charge.end());
        bgc_crd1.assign(bgc_crd.begin(),bgc_crd.end());
        bgc_crd1.insert(bgc_crd1.end(),mol2.coord.begin(),mol2.coord.end());
        bgc_crd2.assign(bgc_crd.begin(),bgc_crd.end());
        bgc_crd2.insert(bgc_crd2.end(),mol1.coord.begin(),mol1.coord.end());
    }
    void form_dimer_from_pbccells(const Pbc_cell_smols & uc, const Pbc_cell_smols & sc,
                               const vector<vector<double > > & dis_mat,
                               const int & i1, const int & i2,
                               const int & stt_index,
                               const double & rcut_sr){
        index1 = i1;
        index2 = i2;
        ele.clear();
        coord.clear();
        charge.clear();
        bgc_crd.clear();
        bgc_chg.clear();
        SingleMol mol1 = sc.smols[i1];
        SingleMol mol2 = sc.smols[i2];
        rij = dis_mat[i1][i2];
        
        com1 = mol1.comass;
        com2 = mol2.comass;
        species1=mol1.species;
        species2=mol2.species;
        mol_index1_uc=mol1.mol_index;
        mol_index2_uc=mol2.mol_index;
        
        sc_index1 = mol1.sc_index;
        sc_index2 = mol2.sc_index;
        
        ele.insert(ele.end(),mol1.ele.begin(),mol1.ele.end());
        ele.insert(ele.end(),mol2.ele.begin(),mol2.ele.end());
        ele1 = mol1.ele;
        ele2 = mol2.ele;
        
        coord.insert(coord.end(),mol1.coord.begin(),mol1.coord.end());
        coord.insert(coord.end(),mol2.coord.begin(),mol2.coord.end());
        coord1 = mol1.coord;
        coord2 = mol2.coord;
        
        charge.insert(charge.end(),mol1.charge.begin(),mol1.charge.end());
        charge.insert(charge.end(),mol2.charge.begin(),mol2.charge.end());
        charge1 = mol1.charge;
        charge2 = mol2.charge;
        
        tc1 = mol1.tot_charge;
        tc2 = mol2.tot_charge;
        ts1 = mol1.tot_spin;
        ts2 = mol2.tot_spin;
        tot_charge = tc1 + tc2;
        tot_spin = 2 * (tot_charge) + 1;
        if (ts2 != 1 or ts1 != 1){
            tot_spin = max(ts1,ts2);
        }
        
        for (int i=0;i<sc.smols.size();i++){
            if ( (i!=i1 and i!=i2) ){
                
                bgc_chg.insert(bgc_chg.end(),sc.smols[i].charge.begin(),sc.smols[i].charge.end());
                bgc_crd.insert(bgc_crd.end(),sc.smols[i].coord.begin(),sc.smols[i].coord.end());
                
            }
        }
        
        bgc_chg1.assign(mol2.charge.begin(),mol2.charge.end());
        bgc_chg1.insert(bgc_chg1.end(),bgc_chg.begin(),bgc_chg.end());
        bgc_chg2.assign(mol1.charge.begin(),mol1.charge.end());
        bgc_chg2.insert(bgc_chg2.end(),bgc_chg.begin(),bgc_chg.end());
        bgc_crd1.assign(mol2.coord.begin(),mol2.coord.end());
        bgc_crd1.insert(bgc_crd1.end(),bgc_crd.begin(),bgc_crd.end());
        bgc_crd2.assign(mol1.coord.begin(),mol1.coord.end());
        bgc_crd2.insert(bgc_crd2.end(),bgc_crd.begin(),bgc_crd.end());
        
    }
    double energy;
    double ene_m1;
    double ene_m2;
    vector<vector<double > > gradient;
    vector<vector<double > > grad_m1;
    vector<vector<double > > grad_m2;
    void set_ene_m1_zero(void){ ene_m1=0.0; };
    void set_ene_m2_zero(void){ ene_m2=0.0; };
    void set_grad_m1_zero(void) { set_matrix_zero(grad_m1, ele1.size(), 3); };
    void set_grad_m2_zero(void) { set_matrix_zero(grad_m2, ele2.size(), 3); };
    vector<vector<double > > hessian;
    vector<vector<double > > hess_m1;
    vector<vector<double > > hess_m2;
    void set_hess_m1_zero(void) { set_matrix_zero(hess_m1, ele1.size()*3, ele1.size()*3); };
    void set_hess_m2_zero(void) { set_matrix_zero(hess_m2, ele2.size()*3, ele2.size()*3); };
    vector<double > dipole_moment;
    vector<double > dp_m1;
    vector<double > dp_m2;
    vector<vector<double > > dipole_moment_derivative;
    vector<vector<double > > dp_der_m1;
    vector<vector<double > > dp_der_m2;
    void set_dp_mnt_m1_zero(void) { set_vector_zero(dp_m1,3); };
    void set_dp_mnt_m2_zero(void) { set_vector_zero(dp_m2,3); };
    void set_dp_der_m1_zero(void) { set_matrix_zero(dp_der_m1, ele1.size(), 9); };
    void set_dp_der_m2_zero(void) { set_matrix_zero(dp_der_m2, ele2.size(), 9); };
    double trans_dip_mnt;
    double trans_dp_m1;
    double trans_dp_m2;
    void set_tdp_m1_zero(void) { trans_dp_m1=0.0; };
    void set_tdp_m2_zero(void) { trans_dp_m2=0.0; };
};

class Trimer{
public:
    double rij;
    double rik;
    double rjk;
    int index1;
    int index2;
    int index3;
    int tot_charge;
    int tot_spin;
    int tc1;
    int tc2;
    int tc3;
    int ts1;
    int ts2;
    int ts3;
    int species1;
    int species2;
    int species3;
    int mol_index1_uc;
    int mol_index2_uc;
    int mol_index3_uc;
    vector<int > sc_index1;
    vector<int > sc_index2;
    vector<int > sc_index3;
    vector<int > ele;
    vector<vector<double > > coord;
    vector<double > charge;
    vector<double > bgc_chg;
    vector<vector<double > > bgc_crd;
    vector<double > com1;
    vector<double > com2;
    vector<double > com3;
    Trimer(const vector<SingleMol> & smols, int i1, int i2, int i3,
           const vector<vector<double > > & dis_mat,
           double rcut_qm, double rcut_sr){
        form_trimer_from_smols(smols,i1,i2,i3,dis_mat,rcut_qm,rcut_sr);
    }
    Trimer(const Pbc_cell_smols & uc, const Pbc_cell_smols & sc, const int & i, const int & j, const int & k,
            const vector<vector<double > > & dis_mat, const int & stt_index){
        form_trimer_from_pbccells(uc, sc, i, j, k, dis_mat, stt_index);
    }
    Monomer get_m1(void){
        Monomer m1(ele1,coord1,bgc_chg1,bgc_crd1,tc1,ts1,species1);
        return m1;
    }
    Monomer get_m2(void){
        Monomer m2(ele2,coord2,bgc_chg2,bgc_crd2,tc2,ts2,species2);
        return m2;
    }
    Monomer get_m3(void){
        Monomer m3(ele3,coord3,bgc_chg3,bgc_crd3,tc3,ts3,species3);
    }
    

    void get_energy(const double ene){ energy = ene; };
    void get_ene_m1(const double ene){ ene_m1 = ene; };
    void get_ene_m2(const double ene){ ene_m2 = ene; };
    void get_ene_m3(const double ene){ ene_m3 = ene; };
    void get_ene_d12(const double ene){ ene_d12 = ene; };
    void get_ene_d13(const double ene){ ene_d13 = ene; };
    void get_ene_d23(const double ene){ ene_d23 = ene; };

    double get_energy(void) const { return energy; };
    double get_ene_m1(void) const { return ene_m1; };
    double get_ene_m2(void) const { return ene_m2; };
    double get_ene_m3(void) const { return ene_m3; };
    double get_ene_d12(void) const { return ene_d12; };
    double get_ene_d13(void) const { return ene_d13; };
    double get_ene_d23(void) const { return ene_d23; };

    void get_gradient(const vector<vector<double > > grad){ gradient = grad; };
    void get_grad_m1 (const vector<vector<double > > grad){ grad_m1  = grad; };
    void get_grad_m2 (const vector<vector<double > > grad){ grad_m2  = grad; };
    void get_grad_m3 (const vector<vector<double > > grad){ grad_m3  = grad; };
    void get_grad_d12 (const vector<vector<double > > grad) {grad_d12 = grad; };
    void get_grad_d13 (const vector<vector<double > > grad) {grad_d13 = grad; };
    void get_grad_d23 (const vector<vector<double > > grad) {grad_d23 = grad; };
    const vector<vector<double > > & get_gradient(void) const { return gradient; };
    const vector<vector<double > > & get_grad_m1(void)  const { return grad_m1; };
    const vector<vector<double > > & get_grad_m2(void)  const { return grad_m2; };
    const vector<vector<double > > & get_grad_m3(void)  const { return grad_m3; };
    const vector<vector<double > > & get_grad_d12(void) const { return grad_d12; };
    const vector<vector<double > > & get_grad_d13(void) const { return grad_d13; };
    const vector<vector<double > > & get_grad_d23(void) const { return grad_d23; };

    void get_hessian(const vector<vector<double > > hess){ hessian = hess; };
    void get_hess_m1(const vector<vector<double > > hess){ hess_m1 = hess; };
    void get_hess_m2(const vector<vector<double > > hess){ hess_m2 = hess; };
    void get_hess_m3(const vector<vector<double > > hess){ hess_m3 = hess; };
    void get_hess_d12(const vector<vector<double > > hess){ hess_d12 = hess; };
    void get_hess_d13(const vector<vector<double > > hess){ hess_d13 = hess; };
    void get_hess_d23(const vector<vector<double > > hess){ hess_d23 = hess; };
    const vector<vector<double > > & get_hessian(void) const { return hessian; };
    const vector<vector<double > > & get_hess_m1(void) const { return hess_m1; };
    const vector<vector<double > > & get_hess_m2(void) const { return hess_m2; };
    const vector<vector<double > > & get_hess_m3(void) const { return hess_m3; };
    const vector<vector<double > > & get_hess_d12(void) const { return hess_d12; };
    const vector<vector<double > > & get_hess_d13(void) const { return hess_d13; };
    const vector<vector<double > > & get_hess_d23(void) const { return hess_d23; };

    void get_dipole_moment(const vector<double > dp) { dipole_moment = dp; };
    void get_dp_m1(const vector<double > dp ) { dp_m1 = dp; };
    void get_dp_m2(const vector<double > dp ) { dp_m2 = dp; };
    void get_dp_m3(const vector<double > dp ) { dp_m3 = dp; };
    void get_dp_d12(const vector<double > dp ) { dp_d12 = dp; };
    void get_dp_d13(const vector<double > dp ) { dp_d13 = dp; };
    void get_dp_d23(const vector<double > dp ) { dp_d23 = dp; };
    const vector<double > & get_dipole_moment(void) const { return dipole_moment; };
    const vector<double > & get_dp_m1(void) const { return dp_m1; };
    const vector<double > & get_dp_m2(void) const { return dp_m2; };
    const vector<double > & get_dp_m3(void) const { return dp_m3; };
    const vector<double > & get_dp_d12(void) const { return dp_d12; };
    const vector<double > & get_dp_d13(void) const { return dp_d13; };
    const vector<double > & get_dp_d23(void) const { return dp_d23; };

    void get_dipole_moment_der(const vector<vector<double > > dp_der) { dipole_moment_derivative = dp_der; };
    void get_dp_der_m1(const vector<vector<double > > dp_der) { dp_der_m1 = dp_der; };
    void get_dp_der_m2(const vector<vector<double > > dp_der) { dp_der_m2 = dp_der; };
    void get_dp_der_m3(const vector<vector<double > > dp_der) { dp_der_m3 = dp_der; };
    void get_dp_der_d12(const vector<vector<double > > dp_der){ dp_der_d12= dp_der; };
    void get_dp_der_d13(const vector<vector<double > > dp_der){ dp_der_d13= dp_der; };
    void get_dp_der_d23(const vector<vector<double > > dp_der){ dp_der_d23= dp_der; };
    const vector<vector<double > > & get_dipole_moment_der(void) const { return dipole_moment_derivative; };
    const vector<vector<double > > & get_dp_der_m1(void) const { return dp_der_m1; };
    const vector<vector<double > > & get_dp_der_m2(void) const { return dp_der_m2; };
    const vector<vector<double > > & get_dp_der_m3(void) const { return dp_der_m3; };
    const vector<vector<double > > & get_dp_der_d12(void) const {return dp_der_d12;};
    const vector<vector<double > > & get_dp_der_d13(void) const {return dp_der_d13;};
    const vector<vector<double > > & get_dp_der_d23(void) const {return dp_der_d23;};

    void get_trans_dip_mnt(const double tdp) { trans_dip_mnt = tdp; };
    void get_trans_dp_m1(const double tdp) { trans_dp_m1 = tdp; };
    void get_trans_dp_m2(const double tdp) { trans_dp_m2 = tdp; };
    void get_trans_dp_m3(const double tdp) { trans_dp_m3 = tdp; };
    void get_trans_dp_d12(const double tdp) { trans_dp_d12 = tdp; };
    void get_trans_dp_d13(const double tdp) { trans_dp_d13 = tdp; };
    void get_trans_dp_d23(const double tdp) { trans_dp_d23 = tdp; };
    const double & get_trans_dip_mnt(void) const { return trans_dip_mnt; };
    const double & get_trans_dp_m1(void) const { return trans_dp_m1; };
    const double & get_trans_dp_m2(void) const { return trans_dp_m2; };
    const double & get_trans_dp_m3(void) const { return trans_dp_m3; };
    const double & get_trans_dp_d12(void) const { return trans_dp_d12; };
    const double & get_trans_dp_d13(void) const { return trans_dp_d13; };
    const double & get_trans_dp_d23(void) const { return trans_dp_d23; };

    int get_natom_m1(){ return ele1.size(); };
    int get_natom_m2(){ return ele2.size(); };
    int get_natom_m3(){ return ele3.size(); };

    vector<int > ele1;
    vector<int > ele2;
    vector<int > ele3;
    vector<vector<double > > coord1;
    vector<vector<double > > coord2;
    vector<vector<double > > coord3;
    vector<double > charge1;
    vector<double > charge2;
    vector<double > charge3;
    vector<double > bgc_chg1;
    vector<double > bgc_chg2;
    vector<double > bgc_chg3;
    vector<vector<double > > bgc_crd1;
    vector<vector<double > > bgc_crd2;
    vector<vector<double > > bgc_crd3;

    void set_dd_ene_grad(void) { 
        set_ene_m1_zero();
        set_ene_m2_zero();
        set_ene_m3_zero();
        set_grad_m1_zero();
        set_grad_m2_zero();
        set_grad_m3_zero();
        set_ene_d12_zero();
        set_ene_d13_zero();
        set_ene_d23_zero();
        set_grad_d12_zero();
        set_grad_d13_zero();
        set_grad_d23_zero();
    }

    void set_dd_hess(void) {
        set_hess_m1_zero();
        set_hess_m2_zero();
        set_hess_m3_zero();
        set_hess_d12_zero();
        set_hess_d13_zero();
        set_hess_d23_zero();
    }

    void set_dd_dp_mnt(void) {
        set_dp_mnt_m1_zero();
        set_dp_mnt_m2_zero();
        set_dp_mnt_m3_zero();
        set_dp_mnt_d12_zero();
        set_dp_mnt_d13_zero();
        set_dp_mnt_d23_zero();
    }
    void set_dd_dp_der(void) {
        set_dp_der_m1_zero();
        set_dp_der_m2_zero();
        set_dp_der_m3_zero();
        set_dp_der_d12_zero();
        set_dp_der_d13_zero();
        set_dp_der_d23_zero();
    }

    void set_dd_tdp(void) {
        set_tdp_m1_zero();
        set_tdp_m2_zero();
        set_tdp_m3_zero();
        set_tdp_d12_zero();
        set_tdp_d13_zero();
        set_tdp_d23_zero();
    }
private:
    void form_trimer_from_smols(const vector<SingleMol> & smols,
                                const int & i1, const int & i2, const int & i3,
                                const vector<vector<double > > & dis_mat,
                                const double & rcut_qm, const double & rcut_sr){
        index1 = i1;
        index2 = i2;
        index3 = i3;
        ele.clear();
        coord.clear();
        charge.clear();
        bgc_crd.clear();
        bgc_chg.clear();
        SingleMol mol1 = smols[i1];
        SingleMol mol2 = smols[i2];
        SingleMol mol3 = smols[i3];
        rij = dis_mat[i1][i2];
        rik = dis_mat[i1][i3];
        rij = dis_mat[i2][i3];
        
        com1 = mol1.comass;
        com2 = mol2.comass;
        com3 = mol3.comass;
        
        species1=mol1.species;
        species2=mol2.species;
        species3=mol3.species;
        
        ele.insert(ele.end(),mol1.ele.begin(),mol1.ele.end());
        ele.insert(ele.end(),mol2.ele.begin(),mol2.ele.end());
        ele.insert(ele.end(),mol3.ele.begin(),mol3.ele.end());
        ele1 = mol1.ele;
        ele2 = mol2.ele;
        ele3 = mol3.ele;
        
        coord.insert(coord.end(),mol1.coord.begin(),mol1.coord.end());
        coord.insert(coord.end(),mol2.coord.begin(),mol2.coord.end());
        coord.insert(coord.end(),mol3.coord.begin(),mol3.coord.end());
        coord1 = mol1.coord;
        coord2 = mol2.coord;
        coord3 = mol3.coord;
        
        charge.insert(charge.end(),mol1.charge.begin(),mol1.charge.end());
        charge.insert(charge.end(),mol2.charge.begin(),mol2.charge.end());
        charge.insert(charge.end(),mol3.charge.begin(),mol3.charge.end());
        charge1 = mol1.charge;
        charge2 = mol2.charge;
        charge3 = mol3.charge;
        
        tc1 = mol1.tot_charge;
        tc2 = mol2.tot_charge;
        tc3 = mol3.tot_charge;
        ts1 = mol1.tot_spin;
        ts2 = mol2.tot_spin;
        ts3 = mol3.tot_spin;
        tot_charge = tc1 + tc2 + tc3;
        tot_spin = 2 * (tot_charge) + 1;
        if (ts2 != 1 or ts1 != 1 or ts3 != 1){
            tot_spin = max(ts1,ts2);
        }
        
        for (int i=0;i<smols.size();i++){
            double r1 = dis_mat[i1][i];
            double r2 = dis_mat[i2][i];
            double r3 = dis_mat[i3][i];
            if ( (i!=i1 and i!=i2 and i!=i3) and (r1<rcut_sr or r2<rcut_sr or r3<rcut_sr)){
                bgc_chg.insert(bgc_chg.end(),smols[i].charge.begin(),smols[i].charge.end());
                bgc_crd.insert(bgc_crd.end(),smols[i].coord.begin(),smols[i].coord.end());
            }
        }
        bgc_chg1.assign(bgc_chg.begin(),bgc_chg.end());
        bgc_chg1.insert(bgc_chg1.end(),mol2.charge.begin(),mol2.charge.end());
        bgc_chg1.insert(bgc_chg1.end(),mol3.charge.begin(),mol3.charge.end());
        bgc_chg2.assign(bgc_chg.begin(),bgc_chg.end());
        bgc_chg2.insert(bgc_chg2.end(),mol1.charge.begin(),mol1.charge.end());
        bgc_chg2.insert(bgc_chg2.end(),mol3.charge.begin(),mol3.charge.end());
        bgc_chg3.assign(bgc_chg.begin(),bgc_chg.end());
        bgc_chg3.insert(bgc_chg3.end(),mol1.charge.begin(),mol1.charge.end());
        bgc_chg3.insert(bgc_chg3.end(),mol2.charge.begin(),mol2.charge.end());
        bgc_crd1.assign(bgc_crd.begin(),bgc_crd.end());
        bgc_crd1.insert(bgc_crd1.end(),mol2.coord.begin(),mol2.coord.end());
        bgc_crd1.insert(bgc_crd1.end(),mol3.coord.begin(),mol3.coord.end());
        bgc_crd2.assign(bgc_crd.begin(),bgc_crd.end());
        bgc_crd2.insert(bgc_crd2.end(),mol1.coord.begin(),mol1.coord.end());
        bgc_crd2.insert(bgc_crd2.end(),mol3.coord.begin(),mol3.coord.end());
        bgc_crd3.assign(bgc_crd.begin(),bgc_crd.end());
        bgc_crd3.insert(bgc_crd3.end(),mol1.coord.begin(),mol1.coord.end());
        bgc_crd3.insert(bgc_crd3.end(),mol2.coord.begin(),mol2.coord.end());
    }
    void form_trimer_from_pbccells(const Pbc_cell_smols & uc, const Pbc_cell_smols & sc,
                                   const int & i1, const int & i2, const int & i3,
                                   const vector<vector<double > > & dis_mat, const int & stt_index){
        index1 = i1;
        index2 = i2;
        index3 = i3;
        ele.clear();
        coord.clear();
        charge.clear();
        bgc_crd.clear();
        bgc_chg.clear();
        SingleMol mol1 = sc.smols[i1];
        SingleMol mol2 = sc.smols[i2];
        SingleMol mol3 = sc.smols[i3];
        rij = dis_mat[i1][i2];
        rik = dis_mat[i1][i3];
        rij = dis_mat[i2][i3];
        
        com1 = mol1.comass;
        com2 = mol2.comass;
        com3 = mol3.comass;
        species1=mol1.species;
        species2=mol2.species;
        species3=mol3.species;
        
        sc_index1 = mol1.sc_index;
        sc_index2 = mol2.sc_index;
        sc_index3 = mol3.sc_index;
        
        mol_index1_uc=mol1.mol_index;
        mol_index2_uc=mol2.mol_index;
        mol_index3_uc=mol3.mol_index;
        
        ele.insert(ele.end(),mol1.ele.begin(),mol1.ele.end());
        ele.insert(ele.end(),mol2.ele.begin(),mol2.ele.end());
        ele.insert(ele.end(),mol3.ele.begin(),mol3.ele.end());
        ele1 = mol1.ele;
        ele2 = mol2.ele;
        ele3 = mol3.ele;
        
        coord.insert(coord.end(),mol1.coord.begin(),mol1.coord.end());
        coord.insert(coord.end(),mol2.coord.begin(),mol2.coord.end());
        coord.insert(coord.end(),mol3.coord.begin(),mol3.coord.end());
        coord1 = mol1.coord;
        coord2 = mol2.coord;
        coord3 = mol3.coord;
        
        charge.insert(charge.end(),mol1.charge.begin(),mol1.charge.end());
        charge.insert(charge.end(),mol2.charge.begin(),mol2.charge.end());
        charge.insert(charge.end(),mol3.charge.begin(),mol3.charge.end());
        charge1 = mol1.charge;
        charge2 = mol2.charge;
        charge3 = mol3.charge;
        
        tc1 = mol1.tot_charge;
        tc2 = mol2.tot_charge;
        tc3 = mol3.tot_charge;
        ts1 = mol1.tot_spin;
        ts2 = mol2.tot_spin;
        ts3 = mol3.tot_spin;
        tot_charge = tc1 + tc2 + tc3;
        tot_spin = 2 * (tot_charge) + 1;
        if (ts2 != 1 or ts1 != 1 or ts3 != 1){
            tot_spin = max(ts1,ts2);
        }
        
        for (int i=0;i<sc.smols.size();i++){
            if ( (i!=i1 and i!=i2 and i!=i3) ){
                bgc_chg.insert(bgc_chg.end(),sc.smols[i].charge.begin(),sc.smols[i].charge.end());
                bgc_crd.insert(bgc_crd.end(),sc.smols[i].coord.begin(),sc.smols[i].coord.end());
            }
        }
        //bgc_chg.insert(bgc_chg.end(),sc.external_charges_chg.begin(),sc.external_charges_chg.end());
        //bgc_crd.insert(bgc_crd.end(),sc.external_charges_crd.begin(),sc.external_charges_crd.end());
        bgc_chg1.assign(bgc_chg.begin(),bgc_chg.end());
        bgc_chg1.insert(bgc_chg1.end(),mol2.charge.begin(),mol2.charge.end());
        bgc_chg1.insert(bgc_chg1.end(),mol3.charge.begin(),mol3.charge.end());
        bgc_chg2.assign(bgc_chg.begin(),bgc_chg.end());
        bgc_chg2.insert(bgc_chg2.end(),mol1.charge.begin(),mol1.charge.end());
        bgc_chg2.insert(bgc_chg2.end(),mol3.charge.begin(),mol3.charge.end());
        bgc_chg3.assign(bgc_chg.begin(),bgc_chg.end());
        bgc_chg3.insert(bgc_chg3.end(),mol1.charge.begin(),mol1.charge.end());
        bgc_chg3.insert(bgc_chg3.end(),mol2.charge.begin(),mol2.charge.end());
        bgc_crd1.assign(bgc_crd.begin(),bgc_crd.end());
        bgc_crd1.insert(bgc_crd1.end(),mol2.coord.begin(),mol2.coord.end());
        bgc_crd1.insert(bgc_crd1.end(),mol3.coord.begin(),mol3.coord.end());
        bgc_crd2.assign(bgc_crd.begin(),bgc_crd.end());
        bgc_crd2.insert(bgc_crd2.end(),mol1.coord.begin(),mol1.coord.end());
        bgc_crd2.insert(bgc_crd2.end(),mol3.coord.begin(),mol3.coord.end());
        bgc_crd3.assign(bgc_crd.begin(),bgc_crd.end());
        bgc_crd3.insert(bgc_crd3.end(),mol1.coord.begin(),mol1.coord.end());
        bgc_crd3.insert(bgc_crd3.end(),mol2.coord.begin(),mol2.coord.end());
    }
    double energy;
    double ene_m1;
    double ene_m2;
    double ene_m3;
    double ene_d12;
    double ene_d13;
    double ene_d23;
    vector<vector<double > > gradient;
    vector<vector<double > > grad_m1;
    vector<vector<double > > grad_m2;
    vector<vector<double > > grad_m3;
    vector<vector<double > > grad_d12;
    vector<vector<double > > grad_d13;
    vector<vector<double > > grad_d23;
    void set_ene_m1_zero(void){ ene_m1=0.0; };
    void set_ene_m2_zero(void){ ene_m2=0.0; };
    void set_ene_m3_zero(void){ ene_m3=0.0; };
    void set_grad_m1_zero(void) { set_matrix_zero(grad_m1, ele1.size(), 3); };
    void set_grad_m2_zero(void) { set_matrix_zero(grad_m2, ele2.size(), 3); };
    void set_grad_m3_zero(void) { set_matrix_zero(grad_m3, ele3.size(), 3); };
    void set_ene_d12_zero(void){ ene_d12=0.0; };
    void set_ene_d13_zero(void){ ene_d13=0.0; };
    void set_ene_d23_zero(void){ ene_d23=0.0; };
    void set_grad_d12_zero(void) { set_matrix_zero(grad_d12,(ele1.size()+ele2.size()),3); };
    void set_grad_d13_zero(void) { set_matrix_zero(grad_d13,(ele1.size()+ele3.size()),3); };
    void set_grad_d23_zero(void) { set_matrix_zero(grad_d23,(ele2.size()+ele3.size()),3); };

    vector<vector<double > > hessian;
    vector<vector<double > > hess_m1;
    vector<vector<double > > hess_m2;
    vector<vector<double > > hess_m3;
    vector<vector<double > > hess_d12;
    vector<vector<double > > hess_d13;
    vector<vector<double > > hess_d23;
    void set_hess_m1_zero(void){ set_matrix_zero(hess_m1, ele1.size()*3, ele1.size()*3); };
    void set_hess_m2_zero(void){ set_matrix_zero(hess_m2, ele2.size()*3, ele2.size()*3); };
    void set_hess_m3_zero(void){ set_matrix_zero(hess_m3, ele3.size()*3, ele3.size()*3); };
    void set_hess_d12_zero(void){ set_matrix_zero(hess_d12, (ele1.size()+ele2.size())*3, (ele1.size()+ele2.size())*3); };
    void set_hess_d13_zero(void){ set_matrix_zero(hess_d13, (ele1.size()+ele3.size())*3, (ele1.size()+ele3.size())*3); };
    void set_hess_d23_zero(void){ set_matrix_zero(hess_d23, (ele2.size()+ele3.size())*3, (ele2.size()+ele3.size())*3); };

    vector<double > dipole_moment;
    vector<double > dp_m1;
    vector<double > dp_m2;
    vector<double > dp_m3;
    vector<double > dp_d12;
    vector<double > dp_d13;
    vector<double > dp_d23;
    vector<vector<double > > dipole_moment_derivative;
    vector<vector<double > > dp_der_m1;
    vector<vector<double > > dp_der_m2;
    vector<vector<double > > dp_der_m3;
    vector<vector<double > > dp_der_d12;
    vector<vector<double > > dp_der_d13;
    vector<vector<double > > dp_der_d23;
    void set_dp_mnt_m1_zero(void) { set_vector_zero(dp_m1,3); };
    void set_dp_mnt_m2_zero(void) { set_vector_zero(dp_m2,3); };
    void set_dp_mnt_m3_zero(void) { set_vector_zero(dp_m3,3); };
    void set_dp_mnt_d12_zero(void){ set_vector_zero(dp_d12,3); };
    void set_dp_mnt_d13_zero(void){ set_vector_zero(dp_d13,3); };
    void set_dp_mnt_d23_zero(void){ set_vector_zero(dp_d23,3); };
    void set_dp_der_m1_zero(void) { set_matrix_zero(dp_der_m1, ele1.size(), 9); };
    void set_dp_der_m2_zero(void) { set_matrix_zero(dp_der_m2, ele2.size(), 9); };
    void set_dp_der_m3_zero(void) { set_matrix_zero(dp_der_m3, ele3.size(), 9); };
    void set_dp_der_d12_zero(void){ set_matrix_zero(dp_der_d12, (ele1.size()+ele2.size()),9); };
    void set_dp_der_d13_zero(void){ set_matrix_zero(dp_der_d13, (ele1.size()+ele3.size()),9); };
    void set_dp_der_d23_zero(void){ set_matrix_zero(dp_der_d23, (ele2.size()+ele3.size()),9); };

    double trans_dip_mnt;
    double trans_dp_d12;
    double trans_dp_d13;
    double trans_dp_d23;
    double trans_dp_m1;
    double trans_dp_m2;
    double trans_dp_m3;
    void set_tdp_m1_zero(void) { trans_dp_m1=0.0; };
    void set_tdp_m2_zero(void) { trans_dp_m2=0.0; };
    void set_tdp_m3_zero(void) { trans_dp_m3=0.0; };
    void set_tdp_d12_zero(void) { trans_dp_d12=0.0; };
    void set_tdp_d13_zero(void) { trans_dp_d13=0.0; };
    void set_tdp_d23_zero(void) { trans_dp_d23=0.0; };
};


#endif //MBE_DRIVER_STRUCTS_H

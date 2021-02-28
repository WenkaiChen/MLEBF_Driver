#include "../include/build_minimer.h"


void
build_monomers(vector<Monomer > & monomers,
        const vector<SingleMol > & smols,
        const vector<vector<double > > & dis_mat,
        const double & rcut,
        map<string, int> & map_monomer){
    int nmol = smols.size();
    monomers.clear();

    for (int i=0;i<nmol;i++){
        Monomer monomer(smols,i,dis_mat[i],rcut);
        monomers.push_back(monomer);
        stringstream mono_name;
        mono_name << monomer.index;
        map_monomer[mono_name.str()] = monomer.index;
    }
}

void
build_monomers(vector<Monomer > & monomers,
               Pbc_cell_smols & uc_smols,
               Pbc_cell_smols & sc_smols,
               const int & stt_index,
               vector<int > & num_uc,
               vector<vector<double > > & dis_mat,
               const double & rcut,
               map<string, int> & map_monomer){
    
    int nmol = uc_smols.smols.size();
    monomers.clear();
    for (int i=0; i<nmol; i++){
        Monomer monomer(uc_smols,sc_smols,num_uc,stt_index,i,0,dis_mat[stt_index+i],rcut);
        monomers.push_back(monomer);
        stringstream mono_name;
        mono_name << monomer.index;
        map_monomer[mono_name.str()] = monomer.index;
    }
}

void
build_dimers(vector<Dimer > & dimers,
        const vector<SingleMol > & smols,
        const vector<vector< double> > & dis_mat,
        const double & rcut_qm,
        const double & rcut_sr,
        map<string, int> & map_dimer){
    int nmol = smols.size();
    dimers.clear();
    int count_dimer = 0;
    for (int i=0;i<nmol;i++){
        for (int j=i+1;j<nmol;j++){
            double rij = dis_mat[i][j];
            
            if (rij < rcut_qm){
                Dimer dimer(smols,i,j,dis_mat,rcut_qm,rcut_sr);
                dimers.push_back(dimer);
                stringstream di_name;
                di_name << i << "_" << j;
                map_dimer[di_name.str()] = count_dimer;
                count_dimer += 1;
            }
        }
    }
}

void
build_dimers(vector<Dimer > & dimers,
             Pbc_cell_smols & uc_smols,
             Pbc_cell_smols & sc_smols,
             const int & stt_index,
             vector<int > & num_uc,
             const vector<vector<double > > & dis_mat,
             const double & rcut_qm,
             const double & rcut_sr,
             map<string, int> & map_dimer){
    
    int nmol = uc_smols.smols.size();
    int sc_nmol = sc_smols.smols.size();
    dimers.clear();
    int count_dimer = 0;
    for (int i=stt_index;i<stt_index+nmol;i++){
        for (int j=0;j<stt_index-1;j++){
            double rij = dis_mat[i][j];
            if (rij < rcut_sr){
                Dimer dimer(uc_smols,sc_smols,dis_mat,i,j,stt_index,rcut_sr);
                dimers.push_back(dimer);
                stringstream di_name;
                di_name << i << "_" << j;
                map_dimer[di_name.str()] = count_dimer;
                count_dimer += 1;
            }
        }
        for (int j=i+1; j<stt_index+nmol;j++){
            double rij = dis_mat[i][j];
            if (rij < rcut_sr){
                Dimer dimer(uc_smols,sc_smols,dis_mat,i,j,stt_index,rcut_sr);
                dimers.push_back(dimer);
                stringstream di_name;
                di_name << i << "_" << j;
                map_dimer[di_name.str()] = count_dimer;
                count_dimer += 1;
            }
        }
        for (int j=stt_index+nmol;j<sc_nmol;j++){
            double rij = dis_mat[i][j];
            if (rij < rcut_sr){
                Dimer dimer(uc_smols,sc_smols,dis_mat,i,j,stt_index,rcut_sr);
                dimers.push_back(dimer);
                stringstream di_name;
                di_name << i << "_" << j;
                map_dimer[di_name.str()] = count_dimer;
                count_dimer += 1;
            }
        }
    }
}

void
build_trimers(vector<Trimer > & trimers,
              const vector<SingleMol > & smols,
              const vector<vector<double > > & dis_mat,
              const double & rcut_qm,
              const double & rcut_sr,
              map<string, int> & map_trimer,
              bool reduce_order,
              vector<int> ro_index){
    int nmol = smols.size();
    trimers.clear();

    int count_trimer = 0;
    for (int i=0;i<nmol;i++){
        for (int j=i+1;j<nmol;j++){
            for (int k=j+1;k<nmol;k++){
                double rij = dis_mat[i][j];
                double rik = dis_mat[i][k];
                double rjk = dis_mat[j][k];
                
                if ( (rij < rcut_qm) && (rik < rcut_qm) && (rjk < rcut_qm) ){
                    if (not reduce_order) {
                        Trimer trimer(smols, i, j, k, dis_mat, rcut_qm, rcut_sr);
                        trimers.push_back(trimer);
                        stringstream tri_name;
                        tri_name << i << "_" << j << "_" << k;
                        map_trimer[tri_name.str()] = count_trimer;
                        count_trimer += 1;
                    }
                    else {
                        bool bool_ro;
                        bool_ro =( is_element_in_vector(ro_index, smols[i].species)
                                   or is_element_in_vector(ro_index, smols[j].species)
                                   or is_element_in_vector(ro_index, smols[k].species));
                        if (not bool_ro){
                            Trimer trimer(smols, i, j, k, dis_mat, rcut_qm, rcut_sr);
                            trimers.push_back(trimer);
                            stringstream tri_name;
                            tri_name << i << "_" << j << "_" << k;
                            map_trimer[tri_name.str()] = count_trimer;
                            count_trimer += 1;
                        }
                    }
                }
            }
        }
    }
}

void
build_trimers(vector<Trimer > & trimers,
              Pbc_cell_smols & uc_smols,
              Pbc_cell_smols & sc_smols,
              const int & stt_index,
              vector<int > & num_uc,
              const vector<vector<double > > & dis_mat,
              const double & rcut_qm,
              const double & rcut_sr,
              map<string, int> & map_trimer,
              bool reduce_order,
              vector<int > ro_index) {
    
    int nmol = uc_smols.smols.size();
    int sc_nmol = sc_smols.smols.size();
    trimers.clear();
    int count_trimer = 0;
    for (int i=stt_index; i<stt_index+nmol; i++){
        for (int j=0; j<stt_index-1; j++){
            for (int k=j+1; k<stt_index-1; k++){
                double rij = dis_mat[i][j];
                double rik = dis_mat[i][k];
                double rjk = dis_mat[j][k];
                if ( (rij < rcut_qm) && (rik < rcut_qm) && (rjk < rcut_qm) ){
                    if (not reduce_order) {
                        Trimer trimer(uc_smols, sc_smols, i, j, k, dis_mat, stt_index);
                        trimers.push_back(trimer);
                        stringstream tri_name;
                        tri_name << i << "_" << j << "_" << k;
                        map_trimer[tri_name.str()] = count_trimer;
                        count_trimer += 1;
                    }
                    else {
                        cout << "reduce order method is not implemented for PBC condition.\nProgram will exit with code 1234\n";
                        exit(1234);
                    }
                }
            }
            for (int k=i+1; k<stt_index+nmol; k++){
                double rij = dis_mat[i][j];
                double rik = dis_mat[i][k];
                double rjk = dis_mat[j][k];
                if ( (rij < rcut_qm) && (rik < rcut_qm) && (rjk < rcut_qm) ){
                    if (not reduce_order) {
                        Trimer trimer(uc_smols, sc_smols, i, j, k, dis_mat, stt_index);
                        trimers.push_back(trimer);
                        stringstream tri_name;
                        tri_name << i << "_" << j << "_" << k;
                        map_trimer[tri_name.str()] = count_trimer;
                        count_trimer += 1;
                    }
                    else {
                        cout << "reduce order method is not implemented for PBC condition.\nProgram will exit with code 1234\n";
                        exit(1234);
                    }
                }
            }
            for (int k=stt_index+nmol; k<sc_nmol; k++){
                double rij = dis_mat[i][j];
                double rik = dis_mat[i][k];
                double rjk = dis_mat[j][k];
                if ( (rij < rcut_qm) && (rik < rcut_qm) && (rjk < rcut_qm) ){
                    if (not reduce_order) {
                        Trimer trimer(uc_smols, sc_smols, i, j, k, dis_mat, stt_index);
                        trimers.push_back(trimer);
                        stringstream tri_name;
                        tri_name << i << "_" << j << "_" << k;
                        map_trimer[tri_name.str()] = count_trimer;
                        count_trimer += 1;
                    }
                    else {
                        cout << "reduce order method is not implemented for PBC condition.\nProgram will exit with code 1234\n";
                        exit(1234);
                    }
                }
            }
        }
        for (int j=i+1; j<stt_index+nmol; j++){
            for (int k=j+1; k<stt_index+nmol; k++){
                double rij = dis_mat[i][j];
                double rik = dis_mat[i][k];
                double rjk = dis_mat[j][k];
                if ( (rij < rcut_qm) && (rik < rcut_qm) && (rjk < rcut_qm) ){
                    if (not reduce_order) {
                        Trimer trimer(uc_smols, sc_smols, i, j, k, dis_mat, stt_index);
                        trimers.push_back(trimer);
                        stringstream tri_name;
                        tri_name << i << "_" << j << "_" << k;
                        map_trimer[tri_name.str()] = count_trimer;
                        count_trimer += 1;
                    }
                    else {
                        cout << "reduce order method is not implemented for PBC condition.\nProgram will exit with code 1234\n";
                        exit(1234);
                    }
                }
            }
            for (int k=stt_index+nmol; k<sc_nmol; k++){
                double rij = dis_mat[i][j];
                double rik = dis_mat[i][k];
                double rjk = dis_mat[j][k];
                if ( (rij < rcut_qm) && (rik < rcut_qm) && (rjk < rcut_qm) ){
                    if (not reduce_order) {
                        Trimer trimer(uc_smols, sc_smols, i, j, k, dis_mat, stt_index);
                        trimers.push_back(trimer);
                        stringstream tri_name;
                        tri_name << i << "_" << j << "_" << k;
                        map_trimer[tri_name.str()] = count_trimer;
                        count_trimer += 1;
                    }
                    else {
                        cout << "reduce order method is not implemented for PBC condition.\nProgram will exit with code 1234\n";
                        exit(1234);
                    }
                }
            }
        }
        for (int j=stt_index+nmol;j<sc_nmol;j++){
            for (int k=j+1; k<sc_nmol; k++){
                double rij = dis_mat[i][j];
                double rik = dis_mat[i][k];
                double rjk = dis_mat[j][k];
                if ( (rij < rcut_qm) && (rik < rcut_qm) && (rjk < rcut_qm) ){
                    if (not reduce_order) {
                        Trimer trimer(uc_smols, sc_smols, i, j, k, dis_mat, stt_index);
                        trimers.push_back(trimer);
                        stringstream tri_name;
                        tri_name << i << "_" << j << "_" << k;
                        map_trimer[tri_name.str()] = count_trimer;
                        count_trimer += 1;
                    }
                    else {
                        cout << "reduce order method is not implemented for PBC condition.\nProgram will exit with code 1234\n";
                        exit(1234);
                    }
                }
            }
        }
    }
}

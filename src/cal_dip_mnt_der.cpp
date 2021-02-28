#include "../include/cal_dip_mnt_der.h"
void
cal_dipole_moment_der(
        vector<vector<double > > & dip_mnt_der,
        vector<Monomer > & monomers,
        vector<Dimer > & dimers,
        vector<Trimer > & trimers,
        const double & rcut_qm,
        const double & rcut_sr,
        const double & rcut_lr,
        const map<string, int> & map_monomer,
        const map<string, int> & map_dimer,
        const map<string, int> & map_trimer,
        const bool & directlyDelta,
        const vector<int> & dd_index,
        const bool & reduce_order,
        const vector<int> & ro_index,
        const bool & partial_nobgc,
        const bool & mono_bgc,
        const bool & cal_grad,
        const int & order,
        const bool & use_state,
        const int & tot_state,
        const int & current_state,
        const bool & pbc,
        const int & stt_index){
    
    if (order == 2) {
        read_dp_der_mono(monomers);
        read_dp_der_dimer(dimers,monomers,rcut_qm,rcut_sr,rcut_lr,map_monomer,map_dimer,
                          directlyDelta,dd_index,reduce_order,ro_index,partial_nobgc,mono_bgc,pbc);
        dip_mnt_der = mbe_dipole_moment_der(monomers,dimers,rcut_qm,rcut_sr,rcut_lr,stt_index);
    }
    else if (order == 3) {
        read_dp_der_mono(monomers);
        read_dp_der_dimer(dimers,monomers,rcut_qm,rcut_sr,rcut_lr,map_monomer,map_dimer,
                          directlyDelta,dd_index,reduce_order,ro_index,partial_nobgc,mono_bgc,pbc);
        read_dp_der_trimer(trimers,dimers,monomers,rcut_qm,rcut_sr,rcut_lr,map_monomer,map_dimer,
                       map_trimer,directlyDelta,dd_index,pbc);
        dip_mnt_der = mbe_dipole_moment_der(monomers,dimers,trimers,rcut_qm,rcut_sr,rcut_lr,directlyDelta,dd_index,stt_index);
    }
    write_output_dp_der(dip_mnt_der);
}

void
read_dp_der_mono(vector<Monomer > & monomers){
    int n_mono = monomers.size();
    cout << "\tReading dipole moment derivatives file for monomers...\n";
    for (int i=0; i<n_mono;i++) {
        Monomer &monomer = monomers[i];
        stringstream dp_der_name;

        dp_der_name << "monomer_" << monomer.index << ".dp_der";
        {
            ifstream dp_der_file;
            dp_der_file.open(dp_der_name.str());
            if (!dp_der_file.is_open()) {
                cout << "Could not open the file " << dp_der_name.str() << endl;
                cout << "Program terminating.\n";
                exit(EXIT_FAILURE);
            }
            int natom = monomer.ele.size();
            string line; 
            vector<vector<double > > dipole_moment_der;
            dipole_moment_der.clear();
            dipole_moment_der.reserve(natom);
            dipole_moment_der.resize(natom);
            for (int ii=0;ii<natom;ii++){
                getline(dp_der_file,line);
                vector<string > data = split(line," ");
                dipole_moment_der[ii].reserve(9);
                dipole_moment_der[ii].resize(9);
                for (int jj=0;jj<9;jj++) {
                    dipole_moment_der[ii][jj] = atof(data[jj].c_str());
                }
            }
            monomer.get_dipole_moment_der(dipole_moment_der);
            dp_der_file.close();
        }
    }
}

void
read_dp_der_dimer(vector<Dimer > & dimers,
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
                  const bool & pbc){
    int n_dimer = dimers.size();
    cout << "\tReading dipole moment derivatives files for dimers...\n";
    for (int i=0;i<n_dimer;i++) {
        Dimer &dimer = dimers[i];
        if (dimer.rij < rcut_qm) {
            stringstream dp_der_name;
            dp_der_name << "dimer_" << i << ".dp_der";
            {
                ifstream dp_der_file;
                dp_der_file.open(dp_der_name.str());
                if (!dp_der_file.is_open()) {
                    cout << "Could not open the file " << dp_der_name.str() << endl;
                    cout << "Program terminating.\n";
                    exit(EXIT_FAILURE);
                }
                int natom = dimer.ele.size();
                string line; 
                vector<vector<double> > dipole_moment_der;
                dipole_moment_der.clear();
                dipole_moment_der.reserve(natom);
                dipole_moment_der.resize(natom);
                for (int ii=0;ii<natom;ii++){
                    getline(dp_der_file,line);
                    vector<string > data = split(line," ");
                    dipole_moment_der[ii].reserve(9);
                    dipole_moment_der[ii].resize(9);
                    for (int jj=0;jj<9;jj++) {
                        dipole_moment_der[ii][jj] = atof(data[jj].c_str());
                    }
                }
                dimer.get_dipole_moment_der(dipole_moment_der);
                dp_der_file.close();
            }
            if (not pbc){
                {
                    int i1 = dimer.index1;
                    int i2 = dimer.index2;
                    int mono_i1 = i1;
                    int mono_i2 = i2;
                    bool already_done_mono = false;
                    if (directlyDelta) {
                        bool bool_dd;
                        bool_dd = is_element_in_vector(dd_index, monomers[mono_i1].species) and
                                  is_element_in_vector(dd_index, monomers[mono_i2].species);
                        if (bool_dd) {
                            dimer.set_dd_dp_der();
                            already_done_mono = true;
                        }
                    }
                    if (not already_done_mono) {
                        if (partial_nobgc) {
                            
                            dimer.get_dp_der_m1(monomers[mono_i1].get_dipole_moment_der());
                            if (not is_element_in_vector(ro_index, monomers[mono_i1].species)) {
                                stringstream dp_der_name_nobgc;
                                dp_der_name_nobgc << "dimer_" << i << "_nobgc.dp_der";
                                vector<vector<double > > tmp_dp_mnt_der;
                                cout << "This part is not finished yet, cause partial_nobgc is a bad method." << endl;
                                exit(1);
                            }
                            
                            dimer.get_dp_der_m2(monomers[mono_i2].get_dipole_moment_der());
                            if (not is_element_in_vector(ro_index, monomers[mono_i2].species)) {
                                stringstream dp_der_name_nobgc;
                                dp_der_name_nobgc << "dimer_" << i << "_nobgc.dp_der";
                                vector<vector<double > > tmp_dp_mnt_der;
                                cout << "This part is not finished yet, cause partial_nobgc is a bad method." << endl;
                                exit(1);
                            }
                        } else if (mono_bgc) {
                            
                            if (not is_element_in_vector(ro_index, monomers[mono_i1].species)) {
                                stringstream dp_der_name_bgc;
                                dp_der_name_bgc << "monomer_" << mono_i1 << "_bgc.dp_der";
                                dimer.get_dp_der_m1(read_mono_dp_mnt_der_tmp(dp_der_name_bgc,monomers[mono_i1].ele.size()));
                            } else {
                                dimer.get_dp_der_m1(monomers[mono_i1].get_dipole_moment_der());
                            }
                      
                            if (not is_element_in_vector(ro_index, monomers[mono_i2].species)) {
                                stringstream dp_der_name_bgc;
                                dp_der_name_bgc << "monomer_" << mono_i2 << "_bgc.dp_der";
                                dimer.get_dp_der_m2(read_mono_dp_mnt_der_tmp(dp_der_name_bgc,monomers[mono_i2].ele.size()));
                            } else {
                                dimer.get_dp_der_m2(monomers[mono_i2].get_dipole_moment_der());
                            }
                        } else {
                           
                            dimer.get_dp_der_m1(monomers[mono_i1].get_dipole_moment_der());
                            dimer.get_dp_der_m2(monomers[mono_i2].get_dipole_moment_der());
                        }
                    }
                }
            }
            else{ 
                cout << "Dipole moment calculation for PBC condition is not implemented yet." << endl;
                exit(1);
            }
        }
    }
}

void
read_dp_der_trimer(vector<Trimer > & trimers,
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
                   const bool & pbc){
    int n_trimer = trimers.size();
    cout << "\tReading dipole moment derivatives files for trimers...\n";
    for (int i=0;i<n_trimer;i++){
        Trimer & trimer = trimers[i];
        if ( (trimer.rij < rcut_qm) && (trimer.rik < rcut_qm) && (trimer.rjk < rcut_qm) ){
            stringstream dp_der_name;
            dp_der_name << "trimer_" << i << ".dp_der";
            {
                ifstream dp_der_file;
                dp_der_file.open(dp_der_name.str());
                if (!dp_der_file.is_open()) {
                    cout << "Could not open the file " << dp_der_name.str() << endl;
                    cout << "Program terminating.\n";
                    exit(EXIT_FAILURE);
                }
                int natom = trimer.ele.size();
                string line; 
                vector<vector<double> >  dipole_moment_der;
                set_matrix_zero(dipole_moment_der, natom, 9);
                for (int ii = 0; ii < natom; ii++) {
                    getline(dp_der_file, line);
                    vector<string> data = split(line, " ");
                    for (int jj = 0; jj < 9; jj++) {
                        dipole_moment_der[ii][jj] = atof(data[jj].c_str());
                    }
                }
                trimer.get_dipole_moment_der(dipole_moment_der);
                dp_der_file.close();
            }
            if (not pbc) {
                {
                    int i1 = trimer.index1;
                    int i2 = trimer.index2;
                    int i3 = trimer.index3;
                    int mono_i1 = i1;
                    int mono_i2 = i2;
                    int mono_i3 = i3;
                    stringstream di_ij_name;
                    stringstream di_ik_name;
                    stringstream di_jk_name;
                    di_ij_name << i1 << "_" << i2;
                    di_ik_name << i1 << "_" << i3;
                    di_jk_name << i2 << "_" << i3;
                    int di_12 = map_dimer.at(di_ij_name.str());
                    int di_13 = map_dimer.at(di_ik_name.str());
                    int di_23 = map_dimer.at(di_jk_name.str());
                    if (directlyDelta) {
                        bool bool_dd;
                        bool_dd = (is_element_in_vector(dd_index, monomers[mono_i1].species)
                                   and is_element_in_vector(dd_index, monomers[mono_i2].species)
                                   and is_element_in_vector(dd_index, monomers[mono_i3].species));
                        if (bool_dd) {
                            trimer.set_dd_dp_der();
                        } else {
                            
                            trimer.get_dp_der_m1(monomers[mono_i1].get_dipole_moment_der());
                            trimer.get_dp_der_m2(monomers[mono_i2].get_dipole_moment_der());
                            trimer.get_dp_der_m3(monomers[mono_i3].get_dipole_moment_der());
                            trimer.get_dp_der_d12(dimers[di_12].get_dipole_moment_der());
                            trimer.get_dp_der_d13(dimers[di_13].get_dipole_moment_der());
                            trimer.get_dp_der_d23(dimers[di_23].get_dipole_moment_der());
                        }
                    } else {
                        
                        trimer.get_dp_der_m1(monomers[mono_i1].get_dipole_moment_der());
                        trimer.get_dp_der_m2(monomers[mono_i2].get_dipole_moment_der());
                        trimer.get_dp_der_m3(monomers[mono_i3].get_dipole_moment_der());
                        trimer.get_dp_der_d12(dimers[di_12].get_dipole_moment_der());
                        trimer.get_dp_der_d13(dimers[di_13].get_dipole_moment_der());
                        trimer.get_dp_der_d23(dimers[di_23].get_dipole_moment_der());
                    }
                }
            }
            else { 
                cout << "dipole moment calculation for PBC condition is not implemented yet." << endl;
                exit(1);
            }
        }
    }
}

vector<vector<double > > mbe_dipole_moment_der(
        const vector<Monomer > & monomers,
        vector<Dimer > & dimers,
        const double & rcut_qm,
        const double & rcut_sr,
        const double & rcut_lr,
        const int & stt_index){
    map<int, int> natom_set;
    map<int, int> natom_cum; 
    int ntot_atom = 0;
    for (int i=0;i<monomers.size();i++){
        natom_set[monomers[i].index]=monomers[i].ele.size();
        ntot_atom += monomers[i].ele.size();
    }
    natom_cum[monomers[0].index] = 0;
    for (int i=1;i<monomers.size();i++){
        natom_cum[monomers[i].index] = natom_cum[monomers[i-1].index]+natom_set.at(monomers[i-1].index);
    }

    vector<vector<double > > dp_mnt_der(ntot_atom, vector<double > (9,0));
    
    cal_mono_dp_der(dp_mnt_der,monomers,natom_set,natom_cum);
    
    cal_dimer_dp_der(dp_mnt_der,dimers,rcut_qm,rcut_sr,rcut_lr,natom_set,natom_cum,stt_index);

    return dp_mnt_der;
}

vector<vector<double > > mbe_dipole_moment_der(
        const vector<Monomer > & monomers,
        vector<Dimer > & dimers,
        const vector<Trimer > & trimers,
        const double & rcut_qm,
        const double & rcut_sr,
        const double & rcut_lr,
        const bool & directlyDelta,
        const vector<int> & dd_index,
        const int & stt_index){
    
    map<int, int> natom_set;
    map<int, int> natom_cum; 
    int ntot_atom = 0;
    for (int i=0;i<monomers.size();i++){
        natom_set[monomers[i].index]=monomers[i].ele.size();
        ntot_atom += monomers[i].ele.size();
    }
    natom_cum[monomers[0].index] = 0;
    for (int i=1;i<monomers.size();i++){
        natom_cum[monomers[i].index] = natom_cum[monomers[i-1].index]+natom_set.at(monomers[i-1].index);
    }
    vector<vector<double > > dp_mnt_der(ntot_atom, vector<double > (9,0));
    
    cal_mono_dp_der(dp_mnt_der,monomers,natom_set,natom_cum);
    
    cal_dimer_dp_der(dp_mnt_der,dimers,rcut_qm,rcut_sr,rcut_lr,natom_set,natom_cum,stt_index);
    
    cal_trimer_dp_der(dp_mnt_der,monomers,trimers,rcut_qm,rcut_sr,rcut_lr,natom_set,natom_cum,directlyDelta,dd_index,stt_index);

    return dp_mnt_der;
}

void
cal_mono_dp_der(vector<vector<double > > & dip_mnt_der,
                const vector<Monomer > & monomers,
                const map<int, int> & natom_set,
                const map<int, int> & natom_cum){
    int nmol = monomers.size();
    for (int i=0;i<nmol;i++){
        const Monomer & monomer = monomers[i];
        int index = monomer.index;
        int nstart = natom_cum.at(index);
        int natom = natom_set.at(index);
        vector<vector<double > > dp_mnt_der_mono = monomer.get_dipole_moment_der();
        
        for (int j=0;j<natom;j++){
            for (int k=0;k<9;k++) {
                dip_mnt_der[j + nstart][k] += dp_mnt_der_mono[j][k];
            }
        }
    }
}

void
cal_dimer_dp_der(vector<vector<double > > & dip_mnt_der,
                 vector<Dimer > & dimers,
                 const double & rcut_qm,
                 const double & rcut_sr,
                 const double & rcut_lr,
                 const map<int, int> & natom_set,
                 const map<int, int> & natom_cum,
                 const int & stt_index){
    int nmol = natom_set.size();
    int ndimer = dimers.size();
    for (int i=0;i<ndimer;i++){
        Dimer & dimer = dimers[i];
        int i1 = dimer.index1;
        int i2 = dimer.index2;
        if (stt_index >=0) {
            i1=dimer.mol_index1_uc+stt_index;
            i2=dimer.mol_index2_uc+stt_index;
        }
        bool monomer1_in_uc = false; 
        bool monomer2_in_uc = false;
        if (stt_index >= 0) {
            if ( (dimer.index1 >= stt_index) and (dimer.index1 < (stt_index+nmol)) ){
                monomer1_in_uc = true;
            }
            if ( (dimer.index2 >= stt_index) and (dimer.index2 < (stt_index+nmol)) ){
                monomer2_in_uc = true;
            }
        }
        int ns1 = natom_cum.at(i1);
        int ns2 = natom_cum.at(i2);
        int na1 = natom_set.at(i1);
        int na2 = natom_set.at(i2);
        if (dimer.rij < rcut_qm) {
            vector<vector<double> > dp_mnt_der_dim = dimer.get_dipole_moment_der();
            vector<vector<double> > dp_der_m1 = dimer.get_dp_der_m1();
            vector<vector<double> > dp_der_m2 = dimer.get_dp_der_m2();
            if (monomer1_in_uc or (stt_index < 0)) { // non-pbc or pbc-but-monomer1-in-uc
                for (int j = 0; j < na1; j++) {
                    for (int k = 0; k < 9; k++){
                        dip_mnt_der[j + ns1][k] += (dp_mnt_der_dim[j][k] - dp_der_m1[j][k]);
                    }
                }
            } else { // pbc-but-monomer1-not-in-uc
                // do nothing.
                cout << "PBC condition dipole moment der is not implemented." << endl;
                exit(1);
            }
            if (monomer2_in_uc or (stt_index < 0)) { // non-pbc or pbc-but-monomer2-in-uc
                for (int j = 0; j < na2; j++) {
                    for (int k = 0; k < 9; k++) {
                        dip_mnt_der[j + ns2][k] += (dp_mnt_der_dim[j + na1][k] - dp_der_m2[j][k]);
                    }
                }
            } else { // pbc-but-monomer1-not-in-uc
                // do nothing.
                cout << "PBC condition dipole moment der is not implemented." << endl;
                exit(1);
            }
        }
    }
}

void
cal_trimer_dp_der(vector<vector<double > > & dip_mnt_der,
                  const vector<Monomer > & monomers,
                  const vector<Trimer > & trimers,
                  const double & rcut_qm,
                  const double & rcut_sr,
                  const double & rcut_lr,
                  const map<int, int> & natom_set,
                  const map<int, int> & natom_cum,
                  const bool & directlyDelta,
                  const vector<int> & dd_index,
                  const int & stt_index){
    int nmol = natom_set.size();
    int ntrimer = trimers.size();
    for (int i=0;i<ntrimer;i++){
        const Trimer & trimer = trimers[i];
        int i1 = trimer.index1;
        int i2 = trimer.index2;
        int i3 = trimer.index3;
        if (stt_index >=0) {
            i1=trimer.mol_index1_uc+stt_index;
            i2=trimer.mol_index2_uc+stt_index;
            i3=trimer.mol_index3_uc+stt_index;
        }
        bool monomer1_in_uc = false; 
        bool monomer2_in_uc = false;
        bool monomer3_in_uc = false;
        if (stt_index >= 0) {
            if ( (trimer.index1 >= stt_index) and (trimer.index1 < (stt_index+nmol)) ){
                monomer1_in_uc = true;
            }
            if ( (trimer.index2 >= stt_index) and (trimer.index2 < (stt_index+nmol)) ){
                monomer2_in_uc = true;
            }
            if ( (trimer.index3 >= stt_index) and (trimer.index3 < (stt_index+nmol)) ){
                monomer3_in_uc = true;
            }
        }
        int ns1 = natom_cum.at(i1);
        int ns2 = natom_cum.at(i2);
        int ns3 = natom_cum.at(i3);
        int na1 = natom_set.at(i1);
        int na2 = natom_set.at(i2);
        int na3 = natom_set.at(i3);
        if ( (trimer.rij < rcut_qm) && (trimer.rik < rcut_qm) && (trimer.rjk < rcut_qm) ){
            vector<vector<double > > dp_der_tri = trimer.get_dipole_moment_der();
            vector<vector<double > > dp_der_d12 = trimer.get_dp_der_d12();
            vector<vector<double > > dp_der_d13 = trimer.get_dp_der_d13();
            vector<vector<double > > dp_der_d23 = trimer.get_dp_der_d23();
            vector<vector<double > > dp_der_m1  = trimer.get_dp_der_m1();
            vector<vector<double > > dp_der_m2  = trimer.get_dp_der_m2();
            vector<vector<double > > dp_der_m3  = trimer.get_dp_der_m3();
            if (not directlyDelta) {
                if (monomer1_in_uc or (stt_index < 0) ) { // non-pbc or pbc-but-monomer1-in-uc
                    for (int j = 0; j < na1; j++) {
                        for (int k = 0; k < 9; k++) {
                            dip_mnt_der[j + ns1][k] += (dp_der_tri[j][k] - dp_der_d12[j][k] - dp_der_d13[j][k] + dp_der_m1[j][k]);
                        }
                    }
                }
                if (monomer2_in_uc or (stt_index < 0) ) { // non-pbc or pbc-but-monomer2-in-uc
                    for (int j = 0; j < na2; j++) {
                        for (int k = 0; k < 9; k++) {
                            dip_mnt_der[j + ns2][k] += (dp_der_tri[j + na1][k] - dp_der_d12[j + na1][k] -
                                    dp_der_d23[j][k] + dp_der_m2[j][k]);
                        }
                    }
                }
                if (monomer3_in_uc or (stt_index < 0) ) { // non-pbc or pbc-but-monomer3-in-uc
                    for (int j = 0; j < na3; j++) {
                        for (int k = 0; k < 9; k++) {
                            dip_mnt_der[j + ns3][k] += (dp_der_tri[j + na1 + na2][k] - dp_der_d13[j + na1][k] -
                                    dp_der_d23[j + na2][k] + dp_der_m3[j][k]);
                        }
                    }
                }
            }
            else{ //directly delta
                // indexes
                int i1 = trimer.index1;
                int i2 = trimer.index2;
                int i3 = trimer.index3;
                int mono_i1 = i1;
                int mono_i2 = i2;
                int mono_i3 = i3;
                bool bool_dd;
                bool_dd = ( is_element_in_vector(dd_index, monomers[mono_i1].species)
                            and is_element_in_vector(dd_index, monomers[mono_i2].species)
                            and is_element_in_vector(dd_index, monomers[mono_i3].species));
                if (bool_dd){
                    if (monomer1_in_uc or (stt_index < 0) ) { // non-pbc or pbc-but-monomer1-in-uc
                        for (int j = 0; j < na1; j++) {
                            for (int k = 0; k < 9; k++) {
                                dip_mnt_der[j + ns1][k] += (dp_der_tri[j][k]);
                            }
                        }
                    }
                    if (monomer2_in_uc or (stt_index < 0) ) { // non-pbc or pbc-but-monomer2-in-uc
                        for (int j = 0; j < na2; j++) {
                            for (int k = 0; k < 9; k++) {
                                dip_mnt_der[j + ns2][k] += (dp_der_tri[j + na1][k]);
                            }
                        }
                    }
                    if (monomer3_in_uc or (stt_index < 0) ) { // non-pbc or pbc-but-monomer3-in-uc
                        for (int j = 0; j < na3; j++) {
                            for (int k = 0; k < 9; k++) {
                                dip_mnt_der[j + ns3][k] += (dp_der_tri[j + na1 + na2][k]);
                            }
                        }
                    }
                }
                else{
                    // indexes
                    int i1 = trimer.index1;
                    int i2 = trimer.index2;
                    int i3 = trimer.index3;
                    int mono_i1 = i1;
                    int mono_i2 = i2;
                    int mono_i3 = i3;
                    bool bool_dd;
                    bool_dd = ( is_element_in_vector(dd_index, monomers[mono_i1].species)
                                and is_element_in_vector(dd_index, monomers[mono_i2].species)
                                and is_element_in_vector(dd_index, monomers[mono_i3].species));
                    if (not bool_dd) { // if bool_dd == true, then the dimer and monomer energies are zero. no futher action is needed.
                        // trimer is not directly delta, but part of dimers may be directly delta
                        vector<vector<double > > dd_dp_der_d12(na1+na2, vector<double> (9,0.0));
                        vector<vector<double > > dd_dp_der_d13(na1+na3, vector<double> (9,0.0));
                        vector<vector<double > > dd_dp_der_d23(na2+na3, vector<double> (9,0.0));
                        bool bool_dd_dij = ( is_element_in_vector(dd_index, monomers[mono_i1].species)
                                             and is_element_in_vector(dd_index, monomers[mono_i2].species) );
                        bool bool_dd_dik = ( is_element_in_vector(dd_index, monomers[mono_i1].species)
                                             and is_element_in_vector(dd_index, monomers[mono_i3].species) );
                        bool bool_dd_djk = ( is_element_in_vector(dd_index, monomers[mono_i2].species)
                                             and is_element_in_vector(dd_index, monomers[mono_i3].species) );
                        if (bool_dd_dij) {
                            for (int j = 0; j < na1+na2; j++) {
                                for (int k = 0; k < 9; k++) {
                                    dd_dp_der_d12[j][k] = dp_der_d12[j][k];
                                }
                            }
                        } else {
                            for (int j = 0; j < na1; j++) {
                                for (int k = 0; k < 9; k++) {
                                    dd_dp_der_d12[j][k] = dp_der_d12[j][k]-dp_der_m1[j][k];
                                }
                            }
                            for (int j = 0; j < na2; j++) {
                                for (int k = 0; k < 9; k++) {
                                    dd_dp_der_d12[j + na1][k] = dp_der_d12[j + na1][k]-dp_der_m2[j][k];
                                }
                            }
                        }
                        if (bool_dd_dik) {
                            for (int j = 0; j < na1+na3; j++){
                                for (int k = 0; k < 9; k++) {
                                    dd_dp_der_d13[j][k] = dp_der_d13[j][k];
                                }
                            }
                        } else {
                            for (int j = 0; j < na1; j++) {
                                for (int k = 0; k < 9; k++) {
                                    dd_dp_der_d13[j][k] = dp_der_d13[j][k]-dp_der_m1[j][k];
                                }
                            }
                            for (int j = 0; j < na3; j++) {
                                for (int k = 0; k < 9; k++) {
                                    dd_dp_der_d13[j + na1][k] = dp_der_d13[j + na1][k]-dp_der_m3[j][k];
                                }
                            }
                        }
                        if (bool_dd_djk) {
                            for (int j = 0; j < na2+na3; j++){
                                for (int k = 0; k < 9; k++) {
                                    dd_dp_der_d23[j][k] = dp_der_d23[j][k];
                                }
                            }
                        } else {
                            for (int j = 0; j < na2; j++) {
                                for (int k = 0; k < 9; k++) {
                                    dd_dp_der_d23[j][k] = dp_der_d23[j][k]-dp_der_m2[j][k];
                                }
                            }
                            for (int j = 0; j < na3; j++) {
                                for (int k = 0; k < 9; k++) {
                                    dd_dp_der_d23[j + na2][k] = dp_der_d23[j + na2][k]-dp_der_m3[j][k];
                                }
                            }
                        }
                        // calculate total gradient
                        if (monomer1_in_uc or (stt_index < 0) ) { // non-pbc or pbc-but-monomer1-in-uc
                            for (int j = 0; j < na1; j++) {
                                for (int k = 0; k < 9; k++) {
                                    dip_mnt_der[j + ns1][k] += (dp_der_tri[j][k] - dd_dp_der_d12[j][k] - dd_dp_der_d13[j][k] -
                                            dp_der_m1[j][k]);
                                }
                            }
                        }
                        if (monomer2_in_uc or (stt_index < 0) ) { // non-pbc or pbc-but-monomer2-in-uc
                            for (int j = 0; j < na2; j++) {
                                for (int k = 0; k < 9; k++) {
                                    dip_mnt_der[j + ns2][k] += (dp_der_tri[j + na1][k] - dd_dp_der_d12[j + na1][k] -
                                                         dd_dp_der_d23[j][k] - dp_der_m2[j][k]);
                                }
                            }
                        }
                        if (monomer3_in_uc or (stt_index < 0) ) { // non-pbc or pbc-but-monomer3-in-uc
                            for (int j = 0; j < na3; j++) {
                                for (int k = 0; k < 9; k++) {
                                    dip_mnt_der[j + ns3][k] += (dp_der_tri[j + na1 + na2][k] - dd_dp_der_d13[j + na1][k] -
                                                         dd_dp_der_d23[j + na2][k] -
                                            dp_der_m3[j][k]);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

vector<vector<double > >
read_mono_dp_mnt_der_tmp(const stringstream & dp_der_name, const int & natom){
    ifstream dp_der_file;
    dp_der_file.open(dp_der_name.str());
    if (!dp_der_file.is_open()){
        cout << "Could not open the file " << dp_der_name.str() << endl;
        cout << "Program terminating.\n";
        exit(EXIT_FAILURE);
    }
    string line; // tmp string for storing the data read from file.
    vector<vector<double > > dp_der;
    dp_der.clear();
    dp_der.reserve(natom);
    dp_der.resize(natom);
    for (int ii=0;ii<natom;ii++){
        getline(dp_der_file,line);
        vector<string > data = split(line," ");
        dp_der[ii].reserve(9);
        dp_der[ii].resize(9);
        for (int jj=0;jj<9;jj++){
            dp_der[ii][jj] = atof(data[jj].c_str());
        }
    }
    dp_der_file.close();
    return dp_der;
}

void
write_output_dp_der(const vector<vector<double > > & dip_mnt_der){
    stringstream dp_der_name;
    dp_der_name << "tot_dp_der.dp_der";
    ofstream dp_der_file;
    dp_der_file.open(dp_der_name.str());
    dp_der_file << fixed << right;
    for (int ii = 0; ii < dip_mnt_der.size(); ii++) {
        for (int jj = 0; jj < 9; jj++) {
            dp_der_file << setw(18) << setprecision(9) << dip_mnt_der[ii][jj];
        }
        dp_der_file << "\n";
    }
    dp_der_file.close();
}

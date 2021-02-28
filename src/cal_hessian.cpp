//
// Created by Wen-Kai on 2020/3/3.
//
#include "../include/cal_hessian.h"

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
            const int & order,
            const bool & use_state,
            const int & tot_state,
            const int & current_state,
            const bool & pbc,
            const int & stt_index){
    if (order == 2) {
        read_hessian_mono(monomers);
        read_hessian_dimer(dimers,monomers,rcut_qm,rcut_sr,rcut_lr,map_monomer,map_dimer,
                           directlyDelta,dd_index,reduce_order,ro_index,partial_nobgc,mono_bgc,pbc);
        hess = mbe_hessian(monomers,dimers,rcut_qm,rcut_sr,rcut_lr,stt_index);
    }
    else if (order == 3) { // no short range and long range.
        read_hessian_mono(monomers);
        read_hessian_dimer(dimers,monomers,rcut_qm,rcut_sr,rcut_lr,map_monomer,map_dimer,
                           directlyDelta,dd_index,reduce_order,ro_index,partial_nobgc,mono_bgc,pbc);
        read_hessian_trimer(trimers,dimers,monomers,rcut_qm,rcut_sr,rcut_lr,map_monomer,map_dimer,
                map_trimer,directlyDelta,dd_index,pbc);
        hess = mbe_hessian(monomers,dimers,trimers,rcut_qm,rcut_sr,rcut_lr,directlyDelta,dd_index,stt_index);
    }
    if (has_freeze_mol) {
        set_freeze_hess(monomers,hess);
    }
    write_output_hessian(hess);
}

void
read_hessian_mono(vector<Monomer > & monomers){
    int n_mono = monomers.size();
    cout << "\tReading hessian files for monomers...\n";
    for (int i=0; i<n_mono;i++) {
        Monomer &monomer = monomers[i];
        stringstream hessname;

        hessname << "monomer_" << monomer.index << ".hess";
        {
            ifstream hessfile;
            hessfile.open(hessname.str());
            if (!hessfile.is_open()){
                cout << "Could not open the file " << hessname.str() << endl;
                cout << "Program terminating.\n";
                exit(EXIT_FAILURE);
            }
            int natom = monomer.ele.size();
            string line; 
            vector<vector<double > > hess;
            hess.clear();
            hess.reserve(natom*3);
            hess.resize(natom*3);
            for (int ii=0;ii<natom*3;ii++){
                getline(hessfile,line);
                vector<string > data = split(line," ");
                hess[ii].reserve(natom*3);
                hess[ii].resize(natom*3);
                for (int jj=0;jj<natom*3;jj++){
                    hess[ii][jj] = atof(data[jj].c_str());
                }
            }
            monomer.get_hessian(hess);
            hessfile.close();
        }
    }
}

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
                   const bool & pbc){
    int n_dimer = dimers.size();
    cout << "\tReading hessian files for dimers...\n";
        for (int i=0;i<n_dimer;i++) {
            Dimer &dimer = dimers[i];
        if (dimer.rij < rcut_qm) {
            stringstream hessname;
            hessname << "dimer_" << i << ".hess";
            {
                ifstream hessfile;
                hessfile.open(hessname.str());
                if (!hessfile.is_open()) {
                    cout << "Could not open the file " << hessname.str() << endl;
                    cout << "Program terminating.\n";
                    exit(EXIT_FAILURE);
                }
                int natom = dimer.ele.size();
                string line; 
                vector<vector<double> > hess;
                hess.clear();
                hess.reserve(natom * 3);
                hess.resize(natom * 3);
                for (int ii = 0; ii < natom * 3; ii++) {
                    getline(hessfile, line);
                    vector<string> data = split(line, " ");
                    hess[ii].reserve(natom * 3);
                    hess[ii].resize(natom * 3);
                    for (int jj = 0; jj < natom * 3; jj++) {
                        hess[ii][jj] = atof(data[jj].c_str());
                    }
                }
                dimer.get_hessian(hess);
                hessfile.close();
            }
            if (not pbc){
                {// deal with monomer and dimer hessian.
                    // indexes
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
                            dimer.set_dd_hess();
                            already_done_mono = true;
                        }
                    }
                    if (not already_done_mono) {
                        if (partial_nobgc) {
                            // monomer_1
                            dimer.get_hess_m1(monomers[mono_i1].get_hessian());
                            if (not is_element_in_vector(ro_index, monomers[mono_i1].species)) {
                                stringstream hessname_nobgc;
                                hessname_nobgc << "dimer_" << i << "_nobgc.hess";
                                int natom = dimer.ele.size();
                                vector<vector<double > > temp_hess = read_mono_hess_from_dimer(hessname_nobgc,
                                                                            natom*3, 0,
                                                                            monomers[mono_i1].ele.size()*3);
                                dimer.set_partial_hess(temp_hess, 0, monomers[mono_i1].ele.size()*3);
                            }
                            // monomer_2
                            dimer.get_hess_m2(monomers[mono_i2].get_hessian());
                            if (not is_element_in_vector(ro_index, monomers[mono_i2].species)) {
                                stringstream hessname_nobgc;
                                hessname_nobgc << "dimer_" << i << "_nobgc.hess";
                                int natom = dimer.ele.size();
                                vector<vector<double > > temp_hess = read_mono_hess_from_dimer(hessname_nobgc,
                                                                                               natom*3,
                                                                                               monomers[mono_i1].ele.size()*3,
                                                                                               natom*3);
                                dimer.set_partial_hess(temp_hess, monomers[mono_i1].ele.size()*3, natom*3);
                            }
                        } else if (mono_bgc) {
                            // monomer_1
                            if (not is_element_in_vector(ro_index, monomers[mono_i1].species)) {
                                stringstream hessname_bgc;
                                hessname_bgc << "monomer_" << mono_i1 << "_bgc.hess";
                                int natom = monomers[mono_i1].ele.size();
                                int nh = natom * 3;
                                dimer.get_hess_m1(read_mono_hess_tmp(hessname_bgc, nh));
                            } else {
                                dimer.get_hess_m1(monomers[mono_i1].get_hessian());
                            }
                            // monomer_2
                            if (not is_element_in_vector(ro_index, monomers[mono_i2].species)) {
                                stringstream hessname_bgc;
                                hessname_bgc << "monomer_" << mono_i2 << "_bgc.hess";
                                int natom = monomers[mono_i2].ele.size();
                                int nh = natom * 3;
                                dimer.get_hess_m2(read_mono_hess_tmp(hessname_bgc, nh));
                            } else {
                                dimer.get_hess_m2(monomers[mono_i2].get_hessian());
                            }
                        } else {
                            // hessian
                            dimer.get_hess_m1(monomers[mono_i1].get_hessian());
                            dimer.get_hess_m2(monomers[mono_i2].get_hessian());
                        }
                    }
                }
            }
            else{ 
                cout << "hessian calculation for PBC condition is not implemented yet." << endl;
                exit(1);
            }
        }
    }
}

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
                    const bool & pbc){
    int n_trimer = trimers.size();
    cout << "\tReading hessian files for trimers...\n";
    for (int i=0;i<n_trimer;i++){
        Trimer & trimer = trimers[i];
        if ( (trimer.rij < rcut_qm) && (trimer.rik < rcut_qm) && (trimer.rjk < rcut_qm) ){
            stringstream hessname;

            hessname << "trimer_" << i << ".hess";
            {
                ifstream hessfile;
                hessfile.open(hessname.str());
                if (!hessfile.is_open()) {
                    cout << "Could not open the file " << hessname.str() << endl;
                    cout << "Program terminating.\n";
                    exit(EXIT_FAILURE);
                }
                int natom = trimer.ele.size();
                string line; 
                vector<vector<double> > hess;
                hess.clear();
                hess.reserve(natom*3);
                hess.resize(natom*3);
                for (int ii = 0; ii < natom*3; ii++) {
                    getline(hessfile, line);
                    vector<string> data = split(line, " ");
                    hess[ii].reserve(natom*3);
                    hess[ii].resize(natom*3);
                    for (int jj = 0; jj < natom * 3; jj++) {
                        hess[ii][jj] = atof(data[jj].c_str());
                    }
                }
                trimer.get_hessian(hess);
                hessfile.close();
            }
            if (not pbc) {
                {// deal with monomer and dimer energies and gradients.
                    // indexes
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
                            trimer.set_dd_hess();
                        } else {
                            // hessian
                            trimer.get_hess_m1(monomers[mono_i1].get_hessian());
                            trimer.get_hess_m2(monomers[mono_i2].get_hessian());
                            trimer.get_hess_m3(monomers[mono_i3].get_hessian());
                            trimer.get_hess_d12(dimers[di_12].get_hessian());
                            trimer.get_hess_d13(dimers[di_13].get_hessian());
                            trimer.get_hess_d23(dimers[di_23].get_hessian());
                        }
                    } else {
                        // hessian
                        trimer.get_hess_m1(monomers[mono_i1].get_hessian());
                        trimer.get_hess_m2(monomers[mono_i2].get_hessian());
                        trimer.get_hess_m3(monomers[mono_i3].get_hessian());
                        trimer.get_hess_d12(dimers[di_12].get_hessian());
                        trimer.get_hess_d13(dimers[di_13].get_hessian());
                        trimer.get_hess_d23(dimers[di_23].get_hessian());
                    }
                }
            }
            else { // pbc condition
                cout << "hessian calculation for PBC condition is not implemented yet." << endl;
                exit(1);
            }
        }
    }
}

vector<vector<double > > mbe_hessian(
        const vector<Monomer > & monomers,
        vector<Dimer > & dimers,
        const double & rcut_qm,
        const double & rcut_sr,
        const double & rcut_lr,
        const int & stt_index){
    map<int, int> natom_set;
    map<int, int> natom_cum; // first number of the x^(nd) molecular.
    int ntot_atom = 0;
    for (int i=0;i<monomers.size();i++){
        natom_set[monomers[i].index]=monomers[i].ele.size();
        ntot_atom += monomers[i].ele.size();
    }
    natom_cum[monomers[0].index] = 0;
    for (int i=1;i<monomers.size();i++){
        natom_cum[monomers[i].index] = natom_cum[monomers[i-1].index]+natom_set.at(monomers[i-1].index);
    }

    vector<vector<double > > hess(ntot_atom*3, vector<double > (ntot_atom*3,0));
    cal_mono_hess(hess,monomers,natom_set,natom_cum);
    cal_dimer_hess(hess,dimers,rcut_qm,rcut_sr,rcut_lr,natom_set,natom_cum,stt_index);
    return hess;
}

vector<vector<double > > mbe_hessian(
        const vector<Monomer > & monomers,
        vector<Dimer > & dimers,
        const vector<Trimer > & trimers,
        const double & rcut_qm,
        const double & rcut_sr,
        const double & rcut_lr,
        const bool & directlyDelta,
        const vector<int> & dd_index,
        const int & stt_index){
    // for trimers.
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
    vector<vector<double > > hess(ntot_atom*3, vector<double > (ntot_atom*3,0));
    cal_mono_hess(hess,monomers,natom_set,natom_cum);
    cal_dimer_hess(hess,dimers,rcut_qm,rcut_sr,rcut_lr,natom_set,natom_cum,stt_index);
    cal_trimer_hess(hess,monomers,trimers,rcut_qm,rcut_sr,rcut_lr,natom_set,natom_cum,directlyDelta,dd_index,stt_index);

    return hess;
}

void
cal_mono_hess(vector<vector<double > > & hess,
              const vector<Monomer > & monomers,
              const map<int, int> & natom_set,
              const map<int, int> & natom_cum){
    int nmol = monomers.size();
    for (int i=0;i<nmol;i++){
        const Monomer & monomer = monomers[i];
        int index = monomer.index;
        int nstart = natom_cum.at(index)*3; 
        int natom = natom_set.at(index);
        vector<vector<double > > hess_mono = monomer.get_hessian();
        // record into gradient.
        for (int j=0;j<natom*3;j++){
            for (int k=0;k<natom*3;k++) {
                hess[j + nstart][k + nstart] += hess_mono[j][k];
            }
        }
    }
}

void
cal_dimer_hess(vector<vector<double > > & hess,
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
        int ns1 = natom_cum.at(i1)*3; 
        int ns2 = natom_cum.at(i2)*3; 
        int na1 = natom_set.at(i1);
        int na2 = natom_set.at(i2);
        int nh1 = na1*3; 
        int nh2 = na2*3; 
        if (dimer.rij < rcut_qm) {
            vector<vector<double> > hessian = dimer.get_hessian();
            vector<vector<double> > hess_m1 = dimer.get_hess_m1();
            vector<vector<double> > hess_m2 = dimer.get_hess_m2();
            if (monomer1_in_uc or (stt_index < 0) ) { 
                for (int j = 0; j < nh1; j++) {
                    for (int k = 0; k < nh1; k++) {
                        hess[j + ns1][k + ns1] += (hessian[j][k] - hess_m1[j][k]);
                    }
                    for (int k = nh1; k < nh1+nh2; k++) {
                        hess[j + ns1][k + ns2 - nh1] += hessian[j][k];
                    }
                }
            } else { 
                // do nothing.
            }
            if (monomer2_in_uc or (stt_index < 0) ) { 
                for (int j = 0; j < nh2; j++) {
                    for (int k = 0; k < nh1; k++){
                        hess[j + ns2][k + ns1] += hessian[j + nh1][k];
                    }
                    for (int k = nh1; k < nh1+nh2; k++){
                        hess[j + ns2][k + ns2 - nh1] += (hessian[j + nh1][k] - hess_m2[j][k - nh1]);
                    }
                }
            } else{ 
                // do nothing.
            }
        }
            
        else if (dimer.rij < rcut_sr){ 
            cout << "The calculations for short range hessian is not implemented. sorry." << endl;
            exit(1);
        }
    }
}

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
        int ns1 = natom_cum.at(i1)*3; 
        int ns2 = natom_cum.at(i2)*3; 
        int ns3 = natom_cum.at(i3)*3; 
        int na1 = natom_set.at(i1);
        int na2 = natom_set.at(i2);
        int na3 = natom_set.at(i3);
        int nh1 = na1*3; 
        int nh2 = na2*3; 
        int nh3 = na3*3; 
        if ( (trimer.rij < rcut_qm) && (trimer.rik < rcut_qm) && (trimer.rjk < rcut_qm) ){
            vector<vector<double > > hessian = trimer.get_hessian();
            vector<vector<double > > hess_d12 = trimer.get_hess_d12();
            vector<vector<double > > hess_d13 = trimer.get_hess_d13();
            vector<vector<double > > hess_d23 = trimer.get_hess_d23();
            vector<vector<double > > hess_m1  = trimer.get_hess_m1();
            vector<vector<double > > hess_m2  = trimer.get_hess_m2();
            vector<vector<double > > hess_m3  = trimer.get_hess_m3();
            if (not directlyDelta) {
                if (monomer1_in_uc or (stt_index < 0) ) { 
                    for (int j = 0; j < nh1; j++) {
                        for (int k = 0; k < nh1; k++){
                            hess[j + ns1][k + ns1] += (hessian[j][k] - hess_d12[j][k] - hess_d13[j][k] + hess_m1[j][k]);
                            
                        }
                        for (int k = nh1; k < (nh1+nh2); k++){
                            hess[j + ns1][k + ns2 - nh1] += (hessian[j][k] - hess_d12[j][k]);
                            
                        }
                        for (int k = (nh1+nh2); k < (nh1+nh2+nh3); k++){
                            hess[j + ns1][k + ns3 - nh1 -nh2] += (hessian[j][k] - hess_d13[j][k - nh2]);
                            
                        }
                    }
                }
                if (monomer2_in_uc or (stt_index < 0) ) { 
                    for (int j = 0; j < nh2; j++) {
                        for (int k = 0; k < nh1; k++){
                            hess[j + ns2][k + ns1] += (hessian[j+nh1][k] - hess_d12[j+nh1][k]);
                            
                        }
                        for (int k = nh1; k < (nh1+nh2); k++){
                            hess[j + ns2][k + ns2 - nh1] += (hessian[j+nh1][k] - hess_d12[j+nh1][k] - hess_d23[j][k-nh1] + hess_m2[j][k-nh1]);
                            
                        }
                        for (int k = (nh1+nh2); k < (nh1+nh2+nh3); k++){
                            hess[j + ns2][k + ns3 - nh1 - nh2] += (hessian[j+nh1][k] - hess_d23[j][k-nh1]);
                            
                        }
                    }
                }
                if (monomer3_in_uc or (stt_index < 0) ) { 
                    for (int j = 0; j < nh3; j++){
                        for (int k = 0; k < nh1; k++){
                            hess[j + ns3][k + ns1] += (hessian[j+nh1+nh2][k] - hess_d13[j+nh1][k]);
                            
                        }
                        for (int k = nh1; k < (nh1+nh2); k++){
                            hess[j + ns3][k + ns2 - nh1] += (hessian[j+nh1+nh2][k] - hess_d23[j+nh2][k-nh1]);
                            
                        }
                        for (int k = (nh1+nh2); k < (nh1+nh2+nh3); k++){
                            hess[j + ns3][k + ns3 - nh1 - nh2] += (hessian[j+nh1+nh2][k] - hess_d13[j+nh1][k-nh2] - hess_d23[j+nh2][k-nh1] + hess_m3[j][k-nh1-nh2]);
                            
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
                    if (monomer1_in_uc or (stt_index < 0) ) { 
                        for (int j = 0; j < nh1; j++) {
                            for (int k = 0; k < nh1; k++){
                                hess[j+ns1][k+ns1] += hessian[j][k];
                            }
                            for (int k = nh1; k < (nh1+nh2); k++){
                                hess[j+ns1][k+ns2-nh1] += hessian[j][k];
                            }
                            for (int k = (nh1+nh2); k < (nh1+nh2+nh3); k++){
                                hess[j+ns1][k+ns3-nh1-nh2] += hessian[j][k];
                            }
                        };
                    }
                    if (monomer2_in_uc or (stt_index < 0) ) { 
                        for (int j = 0; j < nh2; j++) {
                            for (int k = 0; k < nh1; k++){
                                hess[j+ns2][k+ns1] += hessian[j+nh1][k];
                            }
                            for (int k = nh1; k < (nh1+nh2); k++){
                                hess[j+ns2][k+ns2-nh1] += hessian[j+nh1][k];
                            }
                            for (int k = (nh1+nh2); k < (nh1+nh2+nh3); k++){
                                hess[j+ns2][k+ns3-nh1-nh2] += hessian[j+nh1][k];
                            }
                        }
                    }
                    if (monomer3_in_uc or (stt_index < 0) ) { 
                        for (int j = 0; j < nh3; j++) {
                            for (int k = 0; k < nh1; k++){
                                hess[j+ns3][k+ns1] += hessian[j+nh1+nh2][k];
                            }
                            for (int k = nh1; k < (nh1+nh2); k++){
                                hess[j+ns3][k+ns2-nh1] += hessian[j+nh1+nh2][k];
                            }
                            for (int k = (nh1+nh2); k < (nh1+nh2+nh3); k++){
                                hess[j+ns3][k+ns3-nh1-nh2] += hessian[j+nh1+nh2][k];
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
                    if (not bool_dd) { 
                        vector<vector<double > > dd_hess_d12(nh1+nh2, vector<double> (nh1+nh2,0.0));
                        vector<vector<double > > dd_hess_d13(nh1+nh3, vector<double> (nh1+nh2,0.0));
                        vector<vector<double > > dd_hess_d23(nh2+nh3, vector<double> (nh1+nh2,0.0));
                        bool bool_dd_dij = ( is_element_in_vector(dd_index, monomers[mono_i1].species)
                                             and is_element_in_vector(dd_index, monomers[mono_i2].species) );
                        bool bool_dd_dik = ( is_element_in_vector(dd_index, monomers[mono_i1].species)
                                             and is_element_in_vector(dd_index, monomers[mono_i3].species) );
                        bool bool_dd_djk = ( is_element_in_vector(dd_index, monomers[mono_i2].species)
                                             and is_element_in_vector(dd_index, monomers[mono_i3].species) );
                        if (bool_dd_dij) {
                            for (int j = 0; j < (nh1+nh2); j++) {
                                for (int k = 0; k < (nh1+nh2); k++){
                                    dd_hess_d12[j][k] = hess_d12[j][k];
                                }
                            }
                        } else {
                            for (int j = 0; j < nh1; j++) {
                                for (int k = 0; k < nh1; k++){
                                    dd_hess_d12[j][k] = hess_d12[j][k] - hess_m1[j][k];
                                }
                                for (int k = nh1; k < (nh1+nh2); k++){
                                    dd_hess_d12[j][k] = hess_d12[j][k];
                                }
                            }
                            for (int j = 0; j < nh2; j++) {
                                for (int k = 0; k < nh1; k++){
                                    dd_hess_d12[j+nh1][k] = hess_d12[j+nh1][k];
                                }
                                for (int k = nh1; k < (nh1+nh2); k++){
                                    dd_hess_d12[j+nh1][k] = hess_d12[j+nh1][k] - hess_m2[j][k-nh1];
                                }
                            }
                        }
                        if (bool_dd_dik) {
                            for (int j = 0; j < (nh1+nh3); j++) {
                                for (int k = 0; k < (nh1+nh3); k++){
                                    dd_hess_d13[j][k] = hess_d13[j][k];
                                }
                            }
                        } else {
                            for (int j = 0; j < nh1; j++) {
                                for (int k = 0; k < nh1; k++){
                                    dd_hess_d13[j][k] = hess_d13[j][k] - hess_m1[j][k];
                                }
                                for (int k = nh1; k < (nh1+nh3); k++){
                                    dd_hess_d13[j][k] = hess_d13[j][k];
                                }
                            }
                            for (int j = 0; j < nh3; j++) {
                                for (int k = 0; k < nh1; k++){
                                    dd_hess_d13[j+nh1][k] = hess_d13[j+nh1][k];
                                }
                                for (int k = nh1; k < (nh1+nh3); k++){
                                    dd_hess_d13[j+nh1][k] = hess_d13[j+nh1][k] - hess_m3[j][k-nh1];
                                }
                            }
                        }
                        if (bool_dd_djk) {
                            for (int j = 0; j < (nh2+nh3); j++) {
                                for (int k = 0; k < (nh2+nh3); k++){
                                    dd_hess_d23[j][k] = hess_d23[j][k];
                                }
                            }
                        } else {
                            for (int j = 0; j < nh2; j++) {
                                for (int k = 0; k < nh2; k++){
                                    dd_hess_d23[j][k] = hess_d23[j][k] - hess_m2[j][k];
                                }
                                for (int k = nh2; k < (nh2+nh3); k++){
                                    dd_hess_d23[j][k] = hess_d23[j][k];
                                }
                            }
                            for (int j = 0; j < nh3; j++) {
                                for (int k = 0; k < nh2; k++){
                                    dd_hess_d23[j+nh2][k] = hess_d23[j+nh2][k];
                                }
                                for (int k = nh2; k < (nh2+nh3); k++){
                                    dd_hess_d23[j+nh2][k] = hess_d13[j+nh2][k] - hess_m3[j][k-nh2];
                                }
                            }
                        }
                        // calculate total hessian
                        if (monomer1_in_uc or (stt_index < 0) ) { 
                            for (int j = 0; j < nh1; j++) {
                                for (int k = 0; k < nh1; k++){
                                    hess[j+ns1][k+ns1] += (hessian[j][k] - dd_hess_d12[j][k] -
                                                           dd_hess_d13[j][k] - hess_m1[j][k]);
                                }
                                for (int k = nh1; k < (nh1+nh2); k++){
                                    hess[j+ns1][k+ns2-nh1] += (hessian[j][k] - dd_hess_d12[j][k]);
                                }
                                for (int k = (nh1+nh2); k < (nh1+nh2+nh3); k++){
                                    hess[j+ns1][k+ns3-nh1-nh2] += (hessian[j][k] - dd_hess_d13[j][k-nh2]);
                                }
                            }
                        }
                        if (monomer2_in_uc or (stt_index < 0) ) { 
                            for (int j = 0; j < nh2; j++) {
                                for (int k = 0; k < nh1; k++){
                                    hess[j+ns2][k+ns1] += (hessian[j+nh1][k] - dd_hess_d12[j+nh1][k]);
                                }
                                for (int k = nh1; k < (nh1+nh2); k++){
                                    hess[j+ns2][k+ns2-nh1] += (hessian[j+nh1][k] - dd_hess_d12[j+nh1][k] -
                                                               dd_hess_d23[j][k-nh1]-hess_m2[j][k-nh1]);
                                }
                                for (int k = (nh1+nh2); k < (nh1+nh2+nh3); k++){
                                    hess[j+ns2][k+ns3-nh1-nh2] += (hessian[j+nh1][k] - dd_hess_d23[j][k-nh1]);
                                }
                            }
                        }
                        if (monomer3_in_uc or (stt_index < 0) ) { 
                            for (int j = 0; j < nh3; j++) {
                                for (int k = 0; k < nh1; k++){
                                    hess[j+ns3][k+ns1] += (hessian[j+nh1+nh2][k] - dd_hess_d13[j+nh1][k]);
                                }
                                for (int k = nh1; k < (nh1+nh2); k++){
                                    hess[j+ns3][k+ns2-nh1] += (hessian[j+nh1+nh2][k] - dd_hess_d23[j+nh2][k-nh1]);
                                }
                                for (int k = (nh1+nh2); k < (nh1+nh2+nh3); k++){
                                    hess[j+ns3][k+ns3-nh1-nh2] += (hessian[j+nh1+nh2][k] - dd_hess_d13[j+nh1][k-nh2] -
                                                                   dd_hess_d23[j+nh2][k-nh1] - hess_m3[j][k-nh1-nh2]);
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
read_mono_hess_from_dimer(const stringstream & hessname,
                          const int & nh, 
                          const int & start_index,
                          const int & final_index){
    ifstream hessfile;
    hessfile.open(hessname.str());
    if (!hessfile.is_open()){
        cout << "Could not open the file " << hessname.str() << endl;
        cout << "Program terminating.\n";
        exit(EXIT_FAILURE);
    }
    string line; 
    vector<vector<double > > hess;
    int nh_mono = final_index - start_index;
    hess.clear();
    hess.resize(nh_mono);
    hess.resize(nh_mono);
    for (int ii=0;ii<nh;ii++){
        getline(hessfile,line);
        hess[ii].reserve(nh_mono);
        hess[ii].resize(nh_mono);
        if ((ii >= start_index) and (ii<final_index)) {
            vector<string> data = split(line, " ");
            for (int jj=start_index; jj<final_index; jj++){
                hess[ii-start_index][jj-start_index] = atof(data[jj].c_str());
            }
        }
    }
    hessfile.close();
    return hess;
}

vector<vector<double > >
read_mono_hess_tmp(const stringstream & hessname, const int & nh){
    ifstream hessfile;
    hessfile.open(hessname.str());
    if (!hessfile.is_open()){
        cout << "Could not open the file " << hessname.str() << endl;
        cout << "Program terminating.\n";
        exit(EXIT_FAILURE);
    }
    string line; 
    vector<vector<double > > hess;
    hess.clear();
    hess.reserve(nh);
    hess.resize(nh);
    for (int ii=0;ii<nh;ii++){
        getline(hessfile,line);
        vector<string > data = split(line," ");
        hess[ii].reserve(nh);
        hess[ii].resize(nh);
        for (int jj=0;jj<nh;jj++){
            hess[ii][jj] = atof(data[jj].c_str());
        }
    }
    hessfile.close();
    return hess;
}

void
write_output_hessian(const vector<vector<double > > & tot_hess){
    stringstream hessname;
    hessname << "tot_hess.hess";
    ofstream hessfile;
    hessfile.open(hessname.str());
    hessfile << fixed << right;
    for (int ii = 0; ii < tot_hess.size(); ii++) {
        for (int jj = 0; jj < tot_hess.size(); jj++) {
            hessfile << setw(18) << setprecision(9) << tot_hess[ii][jj];
        }
        hessfile << "\n";
    }
    hessfile.close();
}

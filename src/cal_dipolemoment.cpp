#include "../include/cal_dipolemoment.h"

void
cal_dipole_moment(vector<double > & dip_mnt,
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
        read_dp_mono(monomers);
        read_dp_dimer(dimers,monomers,rcut_qm,rcut_sr,rcut_lr,map_monomer,map_dimer,
                      directlyDelta,dd_index,reduce_order,ro_index,partial_nobgc,mono_bgc,pbc);
        dip_mnt = mbe_dipole_moment(monomers,dimers,rcut_qm,rcut_sr,rcut_lr,stt_index);
    }
    else if (order == 3) {
        read_dp_mono(monomers);
        read_dp_dimer(dimers,monomers,rcut_qm,rcut_sr,rcut_lr,map_monomer,map_dimer,
                           directlyDelta,dd_index,reduce_order,ro_index,partial_nobgc,mono_bgc,pbc);
        read_dp_trimer(trimers,dimers,monomers,rcut_qm,rcut_sr,rcut_lr,map_monomer,map_dimer,
                            map_trimer,directlyDelta,dd_index,pbc);
        dip_mnt = mbe_dipole_moment(monomers,dimers,trimers,rcut_qm,rcut_sr,rcut_lr,directlyDelta,dd_index,stt_index);
    }
    write_output_dp(dip_mnt);
}

void
read_dp_mono(vector<Monomer > & monomers){
    int n_mono = monomers.size();
    cout << "\tReading dipole moments file for monomers...\n";
    for (int i=0; i<n_mono;i++) {
        Monomer &monomer = monomers[i];
        stringstream dpname;

        dpname << "monomer_" << monomer.index << ".dp";
        {
            ifstream dpfile;
            dpfile.open(dpname.str());
            if (!dpfile.is_open()) {
                cout << "Could not open the file " << dpname.str() << endl;
                cout << "Program terminating.\n";
                exit(EXIT_FAILURE);
            }
            string line; 
            vector<double > dipole_moment;
            dipole_moment.clear();
            dipole_moment.reserve(3);
            dipole_moment.resize(3);
            getline(dpfile,line);
            vector<string > data = split(line," ");
            for (int ii=0;ii<3;ii++){
                dipole_moment[ii] = atof(data[ii].c_str());
            }
            monomer.get_dipole_moment(dipole_moment);
            dpfile.close();
        }
    }
}

void
read_dp_dimer(vector<Dimer > & dimers,
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
    cout << "\tReading dipole moments files for dimers...\n";
    for (int i=0;i<n_dimer;i++) {
        Dimer &dimer = dimers[i];
        if (dimer.rij < rcut_qm) {
            stringstream dpname;
            dpname << "dimer_" << i << ".dp";
            {
                ifstream dpfile;
                dpfile.open(dpname.str());
                if (!dpfile.is_open()) {
                    cout << "Could not open the file " << dpname.str() << endl;
                    cout << "Program terminating.\n";
                    exit(EXIT_FAILURE);
                }
                int natom = dimer.ele.size();
                string line; 
                vector<double> dipole_moment;
                dipole_moment.clear();
                dipole_moment.reserve(3);
                dipole_moment.resize(3);
                getline(dpfile,line);
                vector<string > data = split(line," ");
                for (int ii=0;ii<3;ii++){
                    dipole_moment[ii] = atof(data[ii].c_str());
                }
                dimer.get_dipole_moment(dipole_moment);
                dpfile.close();
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
                            dimer.set_dd_dp_mnt();
                            already_done_mono = true;
                        }
                    }
                    if (not already_done_mono) {
                        if (partial_nobgc) {
                            // monomer_1
                            dimer.get_dp_m1(monomers[mono_i1].get_dipole_moment());
                            if (not is_element_in_vector(ro_index, monomers[mono_i1].species)) {
                                stringstream dpname_nobgc;
                                dpname_nobgc << "dimer_" << i << "_nobgc.dp";
                                vector<double > other_dp_mnt;
                                if (is_element_in_vector(ro_index, monomers[mono_i2].species)) {
                                    stringstream dpname_m2_nobgc;
                                    dpname_m2_nobgc << "monomer_" << mono_i2 << "_nobgc.dp";
                                    other_dp_mnt = read_mono_dpmnt_tmp(dpname_m2_nobgc);
                                    dimer.get_dipole_moment(other_dp_mnt + dimer.get_dp_m1());
                                }
                            }
                            // monomer_2
                            dimer.get_dp_m2(monomers[mono_i2].get_dipole_moment());
                            if (not is_element_in_vector(ro_index, monomers[mono_i2].species)) {
                                stringstream dpname_nobgc;
                                dpname_nobgc << "dimer_" << i << "_nobgc.dp";
                                vector<double > other_dp_mnt;
                                if (is_element_in_vector(ro_index, monomers[mono_i1].species)){
                                    stringstream dpname_m1_nobgc;
                                    dpname_m1_nobgc << "monomer_" << mono_i1 << "_nobgc.dp";
                                    other_dp_mnt = read_mono_dpmnt_tmp(dpname_m1_nobgc);
                                    dimer.get_dipole_moment(other_dp_mnt + dimer.get_dp_m2());
                                }
                            }
                        } else if (mono_bgc) {
                            // monomer_1
                            if (not is_element_in_vector(ro_index, monomers[mono_i1].species)) {
                                stringstream dpname_bgc;
                                dpname_bgc << "monomer_" << mono_i1 << "_bgc.dp";
                                dimer.get_dp_m1(read_mono_dpmnt_tmp(dpname_bgc));
                            } else {
                                dimer.get_dp_m1(monomers[mono_i1].get_dipole_moment());
                            }
                            // monomer_2
                            if (not is_element_in_vector(ro_index, monomers[mono_i2].species)) {
                                stringstream dpname_bgc;
                                dpname_bgc << "monomer_" << mono_i2 << "_bgc.dp";
                                dimer.get_dp_m2(read_mono_dpmnt_tmp(dpname_bgc));
                            } else {
                                dimer.get_dp_m2(monomers[mono_i2].get_dipole_moment());
                            }
                        } else {
                            
                            dimer.get_dp_m1(monomers[mono_i1].get_dipole_moment());
                            dimer.get_dp_m2(monomers[mono_i2].get_dipole_moment());
                        }
                    }
                }
            }
            else{ // Periodic boundary condition PBC is true.
                cout << "Dipole moment calculation for PBC condition is not implemented yet." << endl;
                exit(1);
            }
        }
    }
}

void
read_dp_trimer(vector<Trimer > & trimers,
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
    cout << "\tReading dipole moment files for trimers...\n";
    for (int i=0;i<n_trimer;i++){
        Trimer & trimer = trimers[i];
        if ( (trimer.rij < rcut_qm) && (trimer.rik < rcut_qm) && (trimer.rjk < rcut_qm) ){
            stringstream dpname;
            dpname << "trimer_" << i << ".dp";
            {// read dipole moment file, 3 nubmers per line, 1 line.
                ifstream dpfile;
                dpfile.open(dpname.str());
                if (!dpfile.is_open()) {
                    cout << "Could not open the file " << dpname.str() << endl;
                    cout << "Program terminating.\n";
                    exit(EXIT_FAILURE);
                }
                int natom = trimer.ele.size();
                string line; // tmp string for storing the data read from file.
                vector<double> dipole_moment;
                set_vector_zero(dipole_moment,3);
                getline(dpfile, line);
                vector<string> data = split(line, " ");
                for (int ii = 0; ii < 3; ii++) {
                    dipole_moment[ii] = atof(data[ii].c_str());
                }
                trimer.get_dipole_moment(dipole_moment);
                dpfile.close();
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
                            trimer.set_dd_dp_mnt();
                        } else {
                            // dipole moment
                            trimer.get_dp_m1(monomers[mono_i1].get_dipole_moment());
                            trimer.get_dp_m2(monomers[mono_i2].get_dipole_moment());
                            trimer.get_dp_m3(monomers[mono_i3].get_dipole_moment());
                            trimer.get_dp_d12(dimers[di_12].get_dipole_moment());
                            trimer.get_dp_d13(dimers[di_13].get_dipole_moment());
                            trimer.get_dp_d23(dimers[di_23].get_dipole_moment());
                        }
                    } else {
                        // dipole moment
                        trimer.get_dp_m1(monomers[mono_i1].get_dipole_moment());
                        trimer.get_dp_m2(monomers[mono_i2].get_dipole_moment());
                        trimer.get_dp_m3(monomers[mono_i3].get_dipole_moment());
                        trimer.get_dp_d12(dimers[di_12].get_dipole_moment());
                        trimer.get_dp_d13(dimers[di_13].get_dipole_moment());
                        trimer.get_dp_d23(dimers[di_23].get_dipole_moment());
                    }
                }
            }
            else { // pbc condition
                cout << "dipole moment calculation for PBC condition is not implemented yet." << endl;
                exit(1);
            }
        }
    }
}

vector<double > mbe_dipole_moment(
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

    vector<double > dp_mnt(ntot_atom*3,0);
    
    cal_mono_dp(dp_mnt,monomers,natom_set,natom_cum);
    
    cal_dimer_dp(dp_mnt,dimers,rcut_qm,rcut_sr,rcut_lr,natom_set,natom_cum,stt_index);
    
    return dp_mnt;
}

vector<double > mbe_dipole_moment(
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
    vector<double > dp_mnt(ntot_atom*3,0);
    
    cal_mono_dp(dp_mnt,monomers,natom_set,natom_cum);
    
    cal_dimer_dp(dp_mnt,dimers,rcut_qm,rcut_sr,rcut_lr,natom_set,natom_cum,stt_index);
    
    cal_trimer_dp(dp_mnt,monomers,trimers,rcut_qm,rcut_sr,rcut_lr,natom_set,natom_cum,directlyDelta,dd_index,stt_index);

    return dp_mnt;
}

void
cal_mono_dp(vector<double > & dip_mnt,
            const vector<Monomer > & monomers,
            const map<int, int> & natom_set,
            const map<int, int> & natom_cum){
    int nmol = monomers.size();
    for (int i=0;i<nmol;i++){
        const Monomer & monomer = monomers[i];
        int index = monomer.index;
        int nstart = natom_cum.at(index)*3; 
        int natom = natom_set.at(index);
        vector<double > dp_mono = monomer.get_dipole_moment();
        
        for (int k=0;k<3;k++) {
            dip_mnt[k] += dp_mono[k];
        }
    }
}

void
cal_dimer_dp(vector<double > & dip_mnt,
             vector<Dimer > & dimers,
             const double & rcut_qm,
             const double & rcut_sr,
             const double & rcut_lr,
             const map<int, int> & natom_set,
             const map<int, int> & natom_cum,
             const int & stt_index){
    bool pbc = false;
    if (stt_index >= 0){
        pbc = true;
    }
    if (not pbc) {
        int nmol = natom_set.size();
        int ndimer = dimers.size();
        for (int i = 0; i < ndimer; i++) {
            Dimer &dimer = dimers[i];
            int i1 = dimer.index1;
            int i2 = dimer.index2;
            int na1 = natom_set.at(i1);
            int na2 = natom_set.at(i2);
            if (dimer.rij < rcut_qm) {
                vector<double> dipole_moment = dimer.get_dipole_moment();
                vector<double> dp_m1 = dimer.get_dp_m1();
                vector<double> dp_m2 = dimer.get_dp_m2();
                dip_mnt = dip_mnt + (dipole_moment - dp_m1 - dp_m2);
            }
                
            else if (dimer.rij < rcut_sr) { 
                cout << "The calculations for short range dipole moment is not implemented. sorry." << endl;
                exit(1);
            }
            
        }
    }
    else {
        
        cout << "dipole moment calculation for PBC condition is not implemented yet." << endl;
        exit(1);
    }
}

void
cal_trimer_dp(vector<double > & dip_mnt,
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
    // for trimers.
    bool pbc = false;
    if (stt_index >= 0){
        pbc = true;
    }
    vector<double > dp_trimer_qm;
    set_vector_zero(dp_trimer_qm,3);
    if (not pbc) {
        for (int i = 0; i < trimers.size(); i++) {
            Trimer trimer = trimers[i];
            if ((trimer.rij < rcut_qm) && (trimer.rik < rcut_qm) && (trimer.rjk < rcut_qm)) {
                if (not directlyDelta) {
                    vector<double > dp_trim =
                            trimer.get_dipole_moment() - trimer.get_dp_d12() - trimer.get_dp_d13() - trimer.get_dp_d23() +
                            trimer.get_dp_m1() + trimer.get_dp_m2() + trimer.get_dp_m3();
                    dp_trimer_qm = dp_trimer_qm + dp_trim;
                } else {
                    
                    int i1 = trimer.index1;
                    int i2 = trimer.index2;
                    int i3 = trimer.index3;
                    int mono_i1 = i1;
                    int mono_i2 = i2;
                    int mono_i3 = i3;
                    bool bool_dd;
                    bool_dd = (is_element_in_vector(dd_index, monomers[mono_i1].species)
                               and is_element_in_vector(dd_index, monomers[mono_i2].species)
                               and is_element_in_vector(dd_index, monomers[mono_i3].species));
                    vector<double > dp_trim = trimer.get_dipole_moment();
                    if (not bool_dd) { 
                        dp_trim = trimer.get_dipole_moment();
                        bool bool_dd_dij = (is_element_in_vector(dd_index, monomers[mono_i1].species)
                                            and is_element_in_vector(dd_index, monomers[mono_i2].species));
                        bool bool_dd_dik = (is_element_in_vector(dd_index, monomers[mono_i1].species)
                                            and is_element_in_vector(dd_index, monomers[mono_i3].species));
                        bool bool_dd_djk = (is_element_in_vector(dd_index, monomers[mono_i2].species)
                                            and is_element_in_vector(dd_index, monomers[mono_i3].species));
                        if (bool_dd_dij) { dp_trim = dp_trim - trimer.get_dp_d12(); }
                        else { dp_trim = dp_trim - (trimer.get_dp_d12() - trimer.get_dp_m1() - trimer.get_dp_m2()); }
                        if (bool_dd_dik) { dp_trim = dp_trim - trimer.get_dp_d13(); }
                        else { dp_trim = dp_trim - (trimer.get_dp_d13() - trimer.get_dp_m1() - trimer.get_dp_m3()); }
                        if (bool_dd_djk) { dp_trim = dp_trim - trimer.get_dp_d23(); }
                        else { dp_trim = dp_trim - (trimer.get_dp_d23() - trimer.get_dp_m2() - trimer.get_dp_m3()); }
                        dp_trim = dp_trim - (trimer.get_dp_m1() + trimer.get_dp_m2() + trimer.get_dp_m3());
                    }
                    dp_trimer_qm = dp_trimer_qm + dp_trim;
                }
            }
        }
    }
    else {
        
        cout << "dipole moment calculation for PBC condition is not implemented yet." << endl;
        exit(1);
    }
    dip_mnt = dip_mnt + dp_trimer_qm;
}

vector<double >
read_mono_dpmnt_tmp(const stringstream & dpname){
    ifstream dpfile;
    dpfile.open(dpname.str());
    if (!dpfile.is_open()) { 
        cout << "Could not open the file " << dpname.str() << endl;
        cout << "Program terminating.\n";
        exit(EXIT_FAILURE);
    }
    string line; 
    getline(dpfile,line);
    vector<string> data = split(line, " ");
    vector<double > dipole_moment;
    set_vector_zero(dipole_moment, 3);
    for (int i=0; i<3; i++){
        dipole_moment[i] = atof(data[i].c_str());
    }
    dpfile.close();
    return dipole_moment;
}

void
write_output_dp(const vector<double > & dip_mnt){
    stringstream dpname;
    dpname << "tot_dip_mnt.dp";
    ofstream dpfile;
    dpfile.open(dpname.str());
    dpfile << fixed << right;
    for (int ii = 0; ii < 3; ii++) {
        dpfile << setw(18) << setprecision(9) << dip_mnt[ii];
    }
    dpfile << "\n";
    dpfile.close();
}

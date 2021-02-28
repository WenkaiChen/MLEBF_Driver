//
// Created by Wen-Kai on 2020/5/24.
//
#include "../include/cal_transition_dip_mnt.h"

void
cal_trans_dip_mnt(double & trans_dip_mnt,
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
        read_trans_dp_mono(monomers);
        read_trans_dp_dimer(dimers,monomers,rcut_qm,rcut_sr,rcut_lr,map_monomer,map_dimer,
                           directlyDelta,dd_index,reduce_order,ro_index,partial_nobgc,mono_bgc,pbc);
        trans_dip_mnt = mbe_trans_dip_mnt(monomers,dimers,rcut_qm,rcut_sr,rcut_lr,stt_index);
    }
    else if (order == 3) { 
        read_trans_dp_mono(monomers);
        read_trans_dp_dimer(dimers,monomers,rcut_qm,rcut_sr,rcut_lr,map_monomer,map_dimer,
                           directlyDelta,dd_index,reduce_order,ro_index,partial_nobgc,mono_bgc,pbc);
        read_trans_dp_trimer(trimers,dimers,monomers,rcut_qm,rcut_sr,rcut_lr,map_monomer,map_dimer,
                            map_trimer,directlyDelta,dd_index,pbc);
        trans_dip_mnt = mbe_trans_dip_mnt(monomers,dimers,trimers,rcut_qm,rcut_sr,rcut_lr,directlyDelta,dd_index,stt_index);
    }
    write_output_trans_dp(trans_dip_mnt);
}

void
read_trans_dp_mono(vector<Monomer > & monomers){
    int n_mono = monomers.size();
    cout << "\tReading transition dipole moment file for monomers...\n";
    for (int i=0; i<n_mono;i++){
        Monomer & monomer = monomers[i];
        stringstream tdpname;

        tdpname << "monomer_" << i << ".tdp";
        {
            ifstream tdpfile;
            tdpfile.open(tdpname.str());
            if (!tdpfile.is_open()) { // failed to open file
                cout << "Could not open the file " << tdpname.str() << endl;
                cout << "Program terminating.\n";
                exit(EXIT_FAILURE);
            }
            string line; 
            getline(tdpfile,line);
            vector<string> data = split(line, " ");
            double tdp;
            istringstream str_tdp(data[0]);
            str_tdp >> tdp;
            monomer.get_trans_dip_mnt(tdp);
            tdpfile.close();
        }
    }
}

void
read_trans_dp_dimer(vector<Dimer > & dimers,
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
                    const bool & pbc) {
    int n_dimer = dimers.size();
    cout << "\tReading transition dipole moment file for dimers...\n";
    for (int i=0;i<n_dimer;i++) {
        Dimer &dimer = dimers[i];
        if (dimer.rij < rcut_qm) {
            stringstream tdpname;

            tdpname << "dimer_" << i << ".tdp";
            {
                ifstream tdpfile;
                tdpfile.open(tdpname.str());
                if (!tdpfile.is_open()) { // failed to open file
                    cout << "Could not open the file " << tdpname.str() << endl;
                    cout << "Program terminating.\n";
                    exit(EXIT_FAILURE);
                }
                string line; 
                getline(tdpfile, line);
                vector<string> data = split(line, " ");
                double tdp;
                istringstream str_tdp(data[0]);
                str_tdp >> tdp;
                dimer.get_trans_dip_mnt(tdp);
                tdpfile.close();
            }
            {// deal with monomer and dimer energies and gradients.
                // indexes
                int i1 = dimer.index1;
                int i2 = dimer.index2;
                int mono_i1 = i1;
                int mono_i2 = i2;
                bool already_done_mono = false;
                if (directlyDelta) {
                    bool bool_dd;
                    bool_dd = is_element_in_vector(dd_index,monomers[mono_i1].species) and is_element_in_vector(dd_index, monomers[mono_i2].species);
                    if (bool_dd){
                        dimer.set_dd_tdp();
                        already_done_mono = true;
                    }
                }
                if (not already_done_mono) {
                    if (partial_nobgc) {
                        // monomer_1
                        dimer.get_trans_dp_m1(monomers[mono_i1].get_trans_dip_mnt());
                        if (not is_element_in_vector(ro_index,monomers[mono_i1].species)){
                            stringstream tdpname_nobgc;
                            tdpname_nobgc  << "dimer_" << i << "_nobgc.tdp";
                            double other_tdp;
                            if (is_element_in_vector(ro_index, monomers[mono_i2].species)){
                                stringstream tdpname_m2_nobgc;
                                tdpname_m2_nobgc << "monomer_" << mono_i2 << "_nobgc.tdp";
                                other_tdp = read_mono_trans_dpmnt_tmp(tdpname_m2_nobgc);
                                dimer.get_trans_dip_mnt(other_tdp+dimer.get_trans_dp_m1());
                            }
                        }
                        // monomer_2
                        dimer.get_trans_dp_m2(monomers[mono_i2].get_trans_dip_mnt());
                        if (not is_element_in_vector(ro_index,monomers[mono_i2].species)){
                            stringstream tdpname_nobgc;
                            tdpname_nobgc  << "dimer_" << i << "_nobgc.tdp";
                            int natom = dimer.ele.size();
                            double other_tdp;
                            if (is_element_in_vector(ro_index, monomers[mono_i1].species)){
                                stringstream tdpname_m1_nobgc;
                                tdpname_m1_nobgc << "monomer_" << mono_i1 << "_nobgc.tdp";
                                other_tdp = read_mono_trans_dpmnt_tmp(tdpname_m1_nobgc);
                                dimer.get_trans_dip_mnt(other_tdp+dimer.get_trans_dp_m2());
                            }
                        }
                    }
                    else if (mono_bgc) {
                        // monomer_1
                        if (not is_element_in_vector(ro_index,monomers[mono_i1].species)){
                            stringstream tdpname_bgc;
                            tdpname_bgc  << "monomer_" << mono_i1 << "_bgc.tdp";
                            dimer.get_trans_dp_m1(read_mono_trans_dpmnt_tmp(tdpname_bgc));
                        } else {
                            dimer.get_trans_dp_m1(monomers[mono_i1].get_trans_dip_mnt());
                        }
                        // monomer_2
                        if (not is_element_in_vector(ro_index,monomers[mono_i2].species)){
                            stringstream tdpname_bgc;
                            tdpname_bgc  << "monomer_" << mono_i2 << "_bgc.tdp";
                            dimer.get_trans_dp_m2(read_mono_trans_dpmnt_tmp(tdpname_bgc));
                        } else {
                            dimer.get_trans_dp_m2(monomers[mono_i2].get_trans_dip_mnt());
                        }
                    }
                    else {
                        // transition dipole moment.
                        dimer.get_trans_dp_m1(monomers[mono_i1].get_trans_dip_mnt());
                        dimer.get_trans_dp_m2(monomers[mono_i2].get_trans_dip_mnt());
                    }
                }
            }
        }
    }
}

void
read_trans_dp_trimer(vector<Trimer > & trimers,
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
                     const bool & pbc) {
    int n_trimer = trimers.size();
    cout << "\tReading transition dipole moment file for trimers...\n";
    for (int i=0;i<n_trimer;i++){
        Trimer & trimer = trimers[i];
        if ( (trimer.rij < rcut_qm) && (trimer.rik < rcut_qm) && (trimer.rjk < rcut_qm) ){
            stringstream tdpname;

            tdpname << "trimer_" << i << ".tdp";
            {
                ifstream tdpfile;
                tdpfile.open(tdpname.str());
                if (!tdpfile.is_open()) { // failed to open file
                    cout << "Could not open the file " << tdpname.str() << endl;
                    cout << "Program terminating.\n";
                    exit(EXIT_FAILURE);
                }
                string line; 
                getline(tdpfile, line);
                vector<string> data = split(line, " ");
                double tdp;
                istringstream str_tdp(data[0]);
                str_tdp >> tdp;
                trimer.get_trans_dip_mnt(tdp);
                tdpfile.close();
            }
            {// deal with monomer and dimer energies.
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
                    bool_dd = ( is_element_in_vector(dd_index, monomers[mono_i1].species)
                                and is_element_in_vector(dd_index, monomers[mono_i2].species)
                                and is_element_in_vector(dd_index, monomers[mono_i3].species));
                    if (bool_dd){
                        trimer.set_dd_tdp();
                    }
                    else {
                        // transition dipole moment
                        trimer.get_trans_dp_m1(monomers[mono_i1].get_trans_dip_mnt());
                        trimer.get_trans_dp_m2(monomers[mono_i2].get_trans_dip_mnt());
                        trimer.get_trans_dp_m3(monomers[mono_i3].get_trans_dip_mnt());
                        trimer.get_trans_dp_d12(dimers[di_12].get_trans_dip_mnt());
                        trimer.get_trans_dp_d13(dimers[di_13].get_trans_dip_mnt());
                        trimer.get_trans_dp_d23(dimers[di_23].get_trans_dip_mnt());
                    }
                }
                else {
                    // transition dipole moment
                    trimer.get_trans_dp_m1(monomers[mono_i1].get_trans_dip_mnt());
                    trimer.get_trans_dp_m2(monomers[mono_i2].get_trans_dip_mnt());
                    trimer.get_trans_dp_m3(monomers[mono_i3].get_trans_dip_mnt());
                    trimer.get_trans_dp_d12(dimers[di_12].get_trans_dip_mnt());
                    trimer.get_trans_dp_d13(dimers[di_13].get_trans_dip_mnt());
                    trimer.get_trans_dp_d23(dimers[di_23].get_trans_dip_mnt());
                }
            }
        }
    }
}

double mbe_trans_dip_mnt(
        const vector<Monomer > & monomers,
        vector<Dimer > & dimers,
        const double & rcut_qm,
        const double & rcut_sr,
        const double & rcut_lr,
        const int & stt_index) {
    bool pbc = false;
    if (stt_index >= 0){
        pbc = true;
    }
    double trans_dip_mnt = 0.0;
    // monomer
    double tdp_monomer = 0.0;
    for (int i=0;i<monomers.size();i++){
        double tdp_mono = monomers[i].get_trans_dip_mnt();
        tdp_monomer += tdp_mono;
    }
    // dimer
    double tdp_dimer_qm = 0.0;
    double tdp_dimer_sr = 0.0;
    double tdp_dimer_lr = 0.0;
    if (not pbc) {
        for (int i = 0; i < dimers.size(); i++) {
            const Dimer &dimer = dimers[i];
            if (dimer.rij < rcut_qm) {
                double tdp_dim = dimer.get_trans_dip_mnt() - dimer.get_trans_dp_m1() - dimer.get_trans_dp_m2();
                tdp_dimer_qm += tdp_dim;
            }
        }
    }
    else { // pbc condition
        //not finished.
        cout << "Transition dipole moment for PBC condition is not finished." << endl;
        exit(1);
    }
    cout << setw(18) << setprecision(9) << "tdp_monomer = " << tdp_monomer << endl;
    cout << setw(18) << setprecision(9) << "tdp_dimer_qm= " << tdp_dimer_qm << endl;
    trans_dip_mnt = tdp_monomer + tdp_dimer_qm;
    return trans_dip_mnt;
}

double mbe_trans_dip_mnt(
        const vector<Monomer > & monomers,
        vector<Dimer > & dimers,
        const vector<Trimer > & trimers,
        const double & rcut_qm,
        const double & rcut_sr,
        const double & rcut_lr,
        const bool & directlyDelta,
        const vector<int> & dd_index,
        const int & stt_index) {
    // for trimers.
    bool pbc = false;
    if (stt_index >= 0){
        pbc = true;
    }
    double trans_dip_mnt = 0.0;
    // monomer
    double tdp_monomer = 0.0;
    for (int i=0;i<monomers.size();i++){
        double tdp_mono = monomers[i].get_trans_dip_mnt();
        tdp_monomer += tdp_mono;
    }
    // dimer
    double tdp_dimer_qm = 0.0;
    if (not pbc) {
        for (int i = 0; i < dimers.size(); i++) {
            const Dimer &dimer = dimers[i];
            if (dimer.rij < rcut_qm) {
                double tdp_dim = dimer.get_trans_dip_mnt() - dimer.get_trans_dp_m1() - dimer.get_trans_dp_m2();
                tdp_dimer_qm += tdp_dim;
            }
        }
    }
    else { // pbc condition
        //not finished.
        cout << "Transition dipole moment for PBC condition is not finished." << endl;
        exit(1);
    }

    // trimer
    double tdp_trimer_qm = 0.0;
    if (not pbc) {
        for (int i = 0; i < trimers.size(); i++) {
            const Trimer &trimer = trimers[i];
            if ((trimer.rij < rcut_qm) && (trimer.rik < rcut_qm) && (trimer.rjk < rcut_qm)) {
                if (not directlyDelta) {
                    double tdp_trim =
                            trimer.get_trans_dip_mnt() - trimer.get_trans_dp_d12() - trimer.get_trans_dp_d13() -
                            trimer.get_trans_dp_d23() + trimer.get_trans_dp_m1() +
                            trimer.get_trans_dp_m2() + trimer.get_trans_dp_m3();
                    tdp_trimer_qm += tdp_trim;
                } else {
                    // indexes
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
                    double tdp_trim = trimer.get_trans_dip_mnt();
                    if (not bool_dd) { 
                        tdp_trim = trimer.get_trans_dip_mnt();
                        bool bool_dd_dij = (is_element_in_vector(dd_index, monomers[mono_i1].species)
                                            and is_element_in_vector(dd_index, monomers[mono_i2].species));
                        bool bool_dd_dik = (is_element_in_vector(dd_index, monomers[mono_i1].species)
                                            and is_element_in_vector(dd_index, monomers[mono_i3].species));
                        bool bool_dd_djk = (is_element_in_vector(dd_index, monomers[mono_i2].species)
                                            and is_element_in_vector(dd_index, monomers[mono_i3].species));
                        if (bool_dd_dij) { tdp_trim -= trimer.get_trans_dp_d12(); }
                        else { tdp_trim -= (trimer.get_trans_dp_d12() - trimer.get_trans_dp_m1() - trimer.get_trans_dp_m2()); }
                        if (bool_dd_dik) { tdp_trim -= trimer.get_trans_dp_d13(); }
                        else { tdp_trim -= (trimer.get_trans_dp_d13() - trimer.get_trans_dp_m1() - trimer.get_trans_dp_m3()); }
                        if (bool_dd_djk) { tdp_trim -= trimer.get_trans_dp_d23(); }
                        else { tdp_trim -= (trimer.get_trans_dp_d23() - trimer.get_trans_dp_m2() - trimer.get_trans_dp_m3()); }
                        tdp_trim -= (trimer.get_trans_dp_m1() + trimer.get_trans_dp_m2() + trimer.get_trans_dp_m3());
                    }
                    tdp_trimer_qm += tdp_trim;
                }
            }
        }
    }
    else { // pbc condition
        //not finished.
        cout << "Transition dipole moment for PBC condition is not finished." << endl;
        exit(1);
    }
    cout << setw(18) << setprecision(9) << "tdp_monomer = " << tdp_monomer << endl;
    cout << setw(18) << setprecision(9) << "tdp_dimer_qm= " << tdp_dimer_qm << endl;
    cout << setw(18) << setprecision(9) << "tdp_trimerqm= " << tdp_trimer_qm << endl;
    trans_dip_mnt = tdp_monomer + tdp_dimer_qm + tdp_trimer_qm;
    return trans_dip_mnt;
}

double
read_mono_trans_dpmnt_tmp(const stringstream & tdpname) {
    ifstream tdpfile;
    tdpfile.open(tdpname.str());
    if (!tdpfile.is_open()) { // failed to open file
        cout << "Could not open the file " << tdpname.str() << endl;
        cout << "Program terminating.\n";
        exit(EXIT_FAILURE);
    }
    string line;
    getline(tdpfile,line);
    vector<string> data = split(line, " ");
    double tdp;
    istringstream str_tdp(data[0]);
    str_tdp >> tdp;
    tdpfile.close();
    return tdp;
}

void
write_output_trans_dp(const double & trans_dip_mnt) {
    stringstream tdpname;
    tdpname << "tot_trans_dip_mnt.tdp";
    cout << "Writing output file ...\n";
    // energy file
    ofstream tdpfile;
    tdpfile.open(tdpname.str());
    tdpfile << fixed << right;
    tdpfile << setw(18) << setprecision(9) << trans_dip_mnt << "\n";
    tdpfile.close();
}

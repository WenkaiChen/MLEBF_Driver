#include "../include/cal_eandg.h"

using namespace std;

void
cal_eandg(double & ene,
          vector<vector<double > > & grad,
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
    // external sp calculation
    cout << "Calculating SP energy and gradients (maybe hessian / dipole moment) from external program package...\n";
    int shell_status = -1;
    if (not use_state) {
        shell_status = system("./MBE_External");
    }
    else {
        char cmd[80];
        sprintf(cmd,"./MBE_External%10d%10d",tot_state,current_state);
        //cout << cmd << endl;
        shell_status = system(cmd);
    }
    
    cout << "Reading energy for monomers and dimers...\n";
    if (cal_grad) {
        read_eg_mono(monomers);
        
        read_eg_dimer(dimers, monomers, rcut_qm, rcut_sr, rcut_lr, map_monomer, map_dimer,
                      directlyDelta, dd_index, reduce_order, ro_index, partial_nobgc, mono_bgc, pbc);
        if (order > 2) {
            read_eg_trimer(trimers, dimers, monomers, rcut_qm, rcut_sr, rcut_lr, map_monomer, map_dimer, map_trimer,
                           directlyDelta, dd_index, pbc);
        }
    }
    else {
        read_only_e_mono(monomers);
        read_only_e_dimer(dimers, monomers, rcut_qm, rcut_sr, rcut_lr, map_monomer, map_dimer,
                          directlyDelta, dd_index, reduce_order, ro_index, partial_nobgc, mono_bgc);
        if (order > 2) {
            read_only_e_trimer(trimers, dimers, monomers, rcut_qm, rcut_sr, rcut_lr, map_monomer, map_dimer, map_trimer,
                               directlyDelta, dd_index);
        }
    }
    cout << "Done.\n";
    
    cout << "Calculating total energy and total gradient...\n";
    if (order == 2) {
        ene = mbe_energy(monomers, dimers, rcut_qm, rcut_sr, rcut_lr, stt_index);
        if (cal_grad) {
            grad = mbe_gradient(monomers, dimers, rcut_qm, rcut_sr, rcut_lr, stt_index);
        }
    }
    else if (order == 3){ // no short range and long range.
        ene = mbe_energy(monomers, dimers, trimers, rcut_qm, rcut_sr, rcut_lr, directlyDelta, dd_index, stt_index);
        if (cal_grad) {
            grad = mbe_gradient(monomers, dimers, trimers, rcut_qm, rcut_sr, rcut_lr, directlyDelta, dd_index, stt_index);
        }
    }
    if (has_freeze_mol and cal_grad) {
        set_freeze_grad(monomers,grad);
    }
    cout << "Done.\n";
}

void
read_eg_mono(vector<Monomer > & monomers){
    int n_mono = monomers.size();
    cout << "\tReading energy file and gradients file for monomers...\n";
    for (int i=0; i<n_mono;i++){
        Monomer & monomer = monomers[i];
        stringstream enename;
        stringstream gradname;

        enename << "monomer_" << monomer.index << ".ene";
        gradname<< "monomer_" << monomer.index << ".grad";
        {
            ifstream enefile;
            enefile.open(enename.str());
            if (!enefile.is_open()) { // failed to open file
                cout << "Could not open the file " << enename.str() << endl;
                cout << "Program terminating.\n";
                exit(EXIT_FAILURE);
            }
            string line; 
            getline(enefile,line);
            vector<string> data = split(line, " ");
            double ene;
            istringstream str_ene(data[0]);
            str_ene >> ene;
            monomer.get_energy(ene);
            enefile.close();
        }
        {
            ifstream gradfile;
            gradfile.open(gradname.str());
            if (!gradfile.is_open()){
                cout << "Could not open the file " << gradname.str() << endl;
                cout << "Program terminating.\n";
                exit(EXIT_FAILURE);
            }
            int natom = monomer.ele.size();
            string line; 
            vector<vector<double > > grad;
            grad.clear();
            for (int ii=0;ii<natom;ii++){
                getline(gradfile,line);
                vector<string > data = split(line," ");
                vector<double > sg;
                sg.push_back(atof(data[0].c_str()));
                sg.push_back(atof(data[1].c_str()));
                sg.push_back(atof(data[2].c_str()));
                grad.push_back(sg);
            }
            monomer.get_gradient(grad);
            gradfile.close();
        }
    }
    
}

void
read_eg_dimer(vector<Dimer > & dimers,
        const double & rcut_qm,
        const double & rcut_sr,
        const double & rcut_lr){
    int n_dimer = dimers.size();
    cout << "\tReading energy file and gradients file for dimers...\n";
    for (int i=0;i<n_dimer;i++){
        Dimer & dimer = dimers[i];
        if (dimer.rij < rcut_qm) {
            stringstream enename;
            stringstream gradname;

            enename << "dimer_" << i << ".ene";
            gradname << "dimer_" << i << ".grad";
            {
                ifstream enefile;
                enefile.open(enename.str());
                if (!enefile.is_open()) { // failed to open file
                    cout << "Could not open the file " << enename.str() << endl;
                    cout << "Program terminating.\n";
                    exit(EXIT_FAILURE);
                }
                string line; 
                getline(enefile, line);
                vector<string> data = split(line, " ");
                double ene;
                istringstream str_ene(data[0]);
                str_ene >> ene;
                dimer.get_energy(ene);
                enefile.close();
            }
            {
                ifstream gradfile;
                gradfile.open(gradname.str());
                if (!gradfile.is_open()) {
                    cout << "Could not open the file " << enename.str() << endl;
                    cout << "Program terminating.\n";
                    exit(EXIT_FAILURE);
                }
                int natom = dimer.ele.size();
                string line; 
                vector<vector<double> > grad;
                grad.clear();
                for (int ii = 0; ii < natom; ii++) {
                    getline(gradfile, line);
                    vector<string> data = split(line, " ");
                    vector<double> sg;
                    sg.push_back(atof(data[0].c_str()));
                    sg.push_back(atof(data[1].c_str()));
                    sg.push_back(atof(data[2].c_str()));
                    grad.push_back(sg);
                }
                dimer.get_gradient(grad);
                gradfile.close();
            }
            {// read dimer_m1
                stringstream enename;
                stringstream gradname;

                enename  << "dimer_" << i << "_m1.ene";
                gradname << "dimer_" << i << "_m1.grad";
                {
                    ifstream enefile;
                    enefile.open(enename.str());
                    if (!enefile.is_open()) { // failed to open file
                        cout << "Could not open the file " << enename.str() << endl;
                        cout << "Program terminating.\n";
                        exit(EXIT_FAILURE);
                    }
                    string line; 
                    getline(enefile, line);
                    vector<string> data = split(line, " ");
                    double ene;
                    istringstream str_ene(data[0]);
                    str_ene >> ene;
                    dimer.get_ene_m1(ene);
                    enefile.close();
                }
                {
                    ifstream gradfile;
                    gradfile.open(gradname.str());
                    if (!gradfile.is_open()) {
                        cout << "Could not open the file " << enename.str() << endl;
                        cout << "Program terminating.\n";
                        exit(EXIT_FAILURE);
                    }
                    int natom = dimer.get_natom_m1();
                    string line; 
                    vector<vector<double> > grad;
                    grad.clear();
                    for (int ii = 0; ii < natom; ii++) {
                        getline(gradfile, line);
                        vector<string> data = split(line, " ");
                        vector<double> sg;
                        sg.push_back(atof(data[0].c_str()));
                        sg.push_back(atof(data[1].c_str()));
                        sg.push_back(atof(data[2].c_str()));
                        grad.push_back(sg);
                    }
                    dimer.get_grad_m1(grad);
                    gradfile.close();
                }
            }
            {// read dimer_m2
                stringstream enename;
                stringstream gradname;

                enename  << "dimer_" << i << "_m2.ene";
                gradname << "dimer_" << i << "_m2.grad";
                {
                    ifstream enefile;
                    enefile.open(enename.str());
                    if (!enefile.is_open()) { // failed to open file
                        cout << "Could not open the file " << enename.str() << endl;
                        cout << "Program terminating.\n";
                        exit(EXIT_FAILURE);
                    }
                    string line; 
                    getline(enefile, line);
                    vector<string> data = split(line, " ");
                    double ene;
                    istringstream str_ene(data[0]);
                    str_ene >> ene;
                    dimer.get_ene_m2(ene);
                    enefile.close();
                }
                {
                    ifstream gradfile;
                    gradfile.open(gradname.str());
                    if (!gradfile.is_open()) {
                        cout << "Could not open the file " << enename.str() << endl;
                        cout << "Program terminating.\n";
                        exit(EXIT_FAILURE);
                    }
                    int natom = dimer.get_natom_m2();
                    string line; 
                    vector<vector<double> > grad;
                    grad.clear();
                    for (int ii = 0; ii < natom; ii++) {
                        getline(gradfile, line);
                        vector<string> data = split(line, " ");
                        vector<double> sg;
                        sg.push_back(atof(data[0].c_str()));
                        sg.push_back(atof(data[1].c_str()));
                        sg.push_back(atof(data[2].c_str()));
                        grad.push_back(sg);
                    }
                    dimer.get_grad_m2(grad);
                    gradfile.close();
                }
            }
            {// calculate the interaction energy
                double interaction_ene;
                interaction_ene = dimer.get_energy() - dimer.get_ene_m1() - dimer.get_ene_m2();
                stringstream dim_int_ene_fn;
                dim_int_ene_fn << "dimer_" << i << ".int";
                ofstream dim_int_f;
                dim_int_f.open(dim_int_ene_fn.str());
                if (dim_int_f.is_open()) {
                    dim_int_f << setw(18) << setprecision(9) << interaction_ene;
                    dim_int_f.close();
                }
            }

        }
        else if (dimer.rij < rcut_sr){
            cal_dimer_gmm(dimer, -1.0);
        }
        else if (dimer.rij < rcut_lr){
            cal_dimer_gmm(dimer,  1.0);
        }
    }
    // cout << "\tDone.\n";
}

void
read_eg_dimer(vector<Dimer > & dimers,
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
    cout << "\tReading energy file and gradients file for dimers...\n";
    for (int i=0;i<n_dimer;i++) {
        Dimer &dimer = dimers[i];
        if (dimer.rij < rcut_qm) {
            stringstream enename;
            stringstream gradname;

            enename << "dimer_" << i << ".ene";
            gradname << "dimer_" << i << ".grad";
            {
                ifstream enefile;
                enefile.open(enename.str());
                if (!enefile.is_open()) { // failed to open file
                    cout << "Could not open the file " << enename.str() << endl;
                    cout << "Program terminating.\n";
                    exit(EXIT_FAILURE);
                }
                string line; 
                getline(enefile, line);
                vector<string> data = split(line, " ");
                double ene;
                istringstream str_ene(data[0]);
                str_ene >> ene;
                dimer.get_energy(ene);
                enefile.close();
            }
            {
                ifstream gradfile;
                gradfile.open(gradname.str());
                if (!gradfile.is_open()) {
                    cout << "Could not open the file " << enename.str() << endl;
                    cout << "Program terminating.\n";
                    exit(EXIT_FAILURE);
                }
                int natom = dimer.ele.size();
                string line; 
                vector<vector<double> > grad;
                grad.clear();
                for (int ii = 0; ii < natom; ii++) {
                    getline(gradfile, line);
                    vector<string> data = split(line, " ");
                    vector<double> sg;
                    sg.push_back(atof(data[0].c_str()));
                    sg.push_back(atof(data[1].c_str()));
                    sg.push_back(atof(data[2].c_str()));
                    grad.push_back(sg);
                }
                dimer.get_gradient(grad);
                gradfile.close();
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
                            dimer.set_dd_ene_grad();
                            already_done_mono = true;
                        }
                    }
                    if (not already_done_mono) {
                        if (partial_nobgc) {
                            // monomer_1
                            dimer.get_ene_m1(monomers[mono_i1].get_energy());
                            dimer.get_grad_m1(monomers[mono_i1].get_gradient());
                            if (not is_element_in_vector(ro_index, monomers[mono_i1].species)) {
                                stringstream gradname_nobgc;
                                stringstream enename_nobgc;
                                enename_nobgc << "dimer_" << i << "_nobgc.ene";
                                gradname_nobgc << "dimer_" << i << "_nobgc.grad";
                                int natom = dimer.ele.size();
                                double other_ene;
                                if (is_element_in_vector(ro_index, monomers[mono_i2].species)) {
                                    stringstream enename_m2_nobgc;
                                    enename_m2_nobgc << "monomer_" << mono_i2 << "_nobgc.ene";
                                    other_ene = read_mono_ene_tmp(enename_m2_nobgc);
                                    dimer.get_energy(other_ene + dimer.get_ene_m1());
                                }
                                vector<vector<double> > temp_grad = read_mono_grad_from_dimer(gradname_nobgc, natom,
                                                                                              0,
                                                                                              monomers[mono_i1].ele.size());
                                dimer.set_partial_grad(temp_grad, 0, monomers[mono_i1].ele.size());
                            }
                            // monomer_2
                            dimer.get_ene_m2(monomers[mono_i2].get_energy());
                            dimer.get_grad_m2(monomers[mono_i2].get_gradient());
                            if (not is_element_in_vector(ro_index, monomers[mono_i2].species)) {
                                stringstream gradname_nobgc;
                                stringstream enename_nobgc;
                                enename_nobgc << "dimer_" << i << "_nobgc.ene";
                                gradname_nobgc << "dimer_" << i << "_nobgc.grad";
                                int natom = dimer.ele.size();
                                double other_ene;
                                if (is_element_in_vector(ro_index, monomers[mono_i1].species)) {
                                    stringstream enename_m1_nobgc;
                                    enename_m1_nobgc << "monomer_" << mono_i1 << "_nobgc.ene";
                                    other_ene = read_mono_ene_tmp(enename_m1_nobgc);
                                    dimer.get_energy(other_ene + dimer.get_ene_m2());
                                }
                                vector<vector<double> > temp_grad = read_mono_grad_from_dimer(gradname_nobgc, natom,
                                                                                              monomers[mono_i1].ele.size(),
                                                                                              natom);
                                dimer.set_partial_grad(temp_grad, monomers[mono_i1].ele.size(), natom);
                            }
                        } else if (mono_bgc) {
                            // monomer_1
                            if (not is_element_in_vector(ro_index, monomers[mono_i1].species)) {
                                stringstream gradname_bgc;
                                stringstream enename_bgc;
                                enename_bgc << "monomer_" << mono_i1 << "_bgc.ene";
                                gradname_bgc << "monomer_" << mono_i1 << "_bgc.grad";
                                int natom = monomers[mono_i1].ele.size();
                                dimer.get_ene_m1(read_mono_ene_tmp(enename_bgc));
                                dimer.get_grad_m1(read_mono_grad_tmp(gradname_bgc, natom));
                            } else {
                                dimer.get_ene_m1(monomers[mono_i1].get_energy());
                                dimer.get_grad_m1(monomers[mono_i1].get_gradient());
                            }
                            // monomer_2
                            if (not is_element_in_vector(ro_index, monomers[mono_i2].species)) {
                                stringstream gradname_bgc;
                                stringstream enename_bgc;
                                enename_bgc << "monomer_" << mono_i2 << "_bgc.ene";
                                gradname_bgc << "monomer_" << mono_i2 << "_bgc.grad";
                                int natom = monomers[mono_i2].ele.size();
                                dimer.get_ene_m2(read_mono_ene_tmp(enename_bgc));
                                dimer.get_grad_m2(read_mono_grad_tmp(gradname_bgc, natom));
                            } else {
                                dimer.get_ene_m2(monomers[mono_i2].get_energy());
                                dimer.get_grad_m2(monomers[mono_i2].get_gradient());
                            }
                        } else {
                            // energies
                            dimer.get_ene_m1(monomers[mono_i1].get_energy());
                            dimer.get_ene_m2(monomers[mono_i2].get_energy());
                            // gradients
                            dimer.get_grad_m1(monomers[mono_i1].get_gradient());
                            dimer.get_grad_m2(monomers[mono_i2].get_gradient());
                        }
                    }
                }
            }
            else{ 
                {// read dimer_m1
                    stringstream enename;
                    stringstream gradname;

                    enename  << "monomer_" << dimer.index1 << ".ene";
                    gradname << "monomer_" << dimer.index1 << ".grad";
                    {
                        ifstream enefile;
                        enefile.open(enename.str());
                        if (!enefile.is_open()) { // failed to open file
                            cout << "Could not open the file " << enename.str() << endl;
                            cout << "Program terminating.\n";
                            exit(EXIT_FAILURE);
                        }
                        string line; 
                        getline(enefile, line);
                        vector<string> data = split(line, " ");
                        double ene;
                        istringstream str_ene(data[0]);
                        str_ene >> ene;
                        dimer.get_ene_m1(ene);
                        enefile.close();
                    }
                    {
                        ifstream gradfile;
                        gradfile.open(gradname.str());
                        if (!gradfile.is_open()) {
                            cout << "Could not open the file " << enename.str() << endl;
                            cout << "Program terminating.\n";
                            exit(EXIT_FAILURE);
                        }
                        int natom = dimer.get_natom_m1();
                        string line; 
                        vector<vector<double> > grad;
                        grad.clear();
                        for (int ii = 0; ii < natom; ii++) {
                            getline(gradfile, line);
                            vector<string> data = split(line, " ");
                            vector<double> sg;
                            sg.push_back(atof(data[0].c_str()));
                            sg.push_back(atof(data[1].c_str()));
                            sg.push_back(atof(data[2].c_str()));
                            grad.push_back(sg);
                        }
                        dimer.get_grad_m1(grad);
                        gradfile.close();
                    }
                }
                {// read dimer_m2
                    stringstream enename;
                    stringstream gradname;

                    enename  << "monomer_" << dimer.index2 << ".ene";
                    gradname << "monomer_" << dimer.index2 << ".grad";
                    {
                        ifstream enefile;
                        enefile.open(enename.str());
                        if (!enefile.is_open()) { // failed to open file
                            cout << "Could not open the file " << enename.str() << endl;
                            cout << "Program terminating.\n";
                            exit(EXIT_FAILURE);
                        }
                        string line; 
                        getline(enefile, line);
                        vector<string> data = split(line, " ");
                        double ene;
                        istringstream str_ene(data[0]);
                        str_ene >> ene;
                        dimer.get_ene_m2(ene);
                        enefile.close();
                    }
                    {
                        ifstream gradfile;
                        gradfile.open(gradname.str());
                        if (!gradfile.is_open()) {
                            cout << "Could not open the file " << enename.str() << endl;
                            cout << "Program terminating.\n";
                            exit(EXIT_FAILURE);
                        }
                        int natom = dimer.get_natom_m2();
                        string line; 
                        vector<vector<double> > grad;
                        grad.clear();
                        for (int ii = 0; ii < natom; ii++) {
                            getline(gradfile, line);
                            vector<string> data = split(line, " ");
                            vector<double> sg;
                            sg.push_back(atof(data[0].c_str()));
                            sg.push_back(atof(data[1].c_str()));
                            sg.push_back(atof(data[2].c_str()));
                            grad.push_back(sg);
                        }
                        dimer.get_grad_m2(grad);
                        gradfile.close();
                    }
                }
            }
        }
    }
}

void
read_eg_trimer(vector<Trimer > & trimers,
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
    cout << "\tReading energy file and gradients file for trimers...\n";
    for (int i=0;i<n_trimer;i++){
        Trimer & trimer = trimers[i];
        if ( (trimer.rij < rcut_qm) && (trimer.rik < rcut_qm) && (trimer.rjk < rcut_qm) ){
            stringstream enename;
            stringstream gradname;

            enename << "trimer_" << i << ".ene";
            gradname << "trimer_" << i << ".grad";
            {
                ifstream enefile;
                enefile.open(enename.str());
                if (!enefile.is_open()) { // failed to open file
                    cout << "Could not open the file " << enename.str() << endl;
                    cout << "Program terminating.\n";
                    exit(EXIT_FAILURE);
                }
                string line; 
                getline(enefile, line);
                vector<string> data = split(line, " ");
                double ene;
                istringstream str_ene(data[0]);
                str_ene >> ene;
                trimer.get_energy(ene);
                enefile.close();
            }
            {
                ifstream gradfile;
                gradfile.open(gradname.str());
                if (!gradfile.is_open()) {
                    cout << "Could not open the file " << enename.str() << endl;
                    cout << "Program terminating.\n";
                    exit(EXIT_FAILURE);
                }
                int natom = trimer.ele.size();
                string line; 
                vector<vector<double> > grad;
                grad.clear();
                for (int ii = 0; ii < natom; ii++) {
                    getline(gradfile, line);
                    vector<string> data = split(line, " ");
                    vector<double> sg;
                    sg.push_back(atof(data[0].c_str()));
                    sg.push_back(atof(data[1].c_str()));
                    sg.push_back(atof(data[2].c_str()));
                    grad.push_back(sg);
                }
                trimer.get_gradient(grad);
                gradfile.close();
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
                            trimer.set_dd_ene_grad();
                        } else {
                            // energies
                            trimer.get_ene_m1(monomers[mono_i1].get_energy());
                            trimer.get_ene_m2(monomers[mono_i2].get_energy());
                            trimer.get_ene_m3(monomers[mono_i3].get_energy());
                            trimer.get_ene_d12(dimers[di_12].get_energy());
                            trimer.get_ene_d13(dimers[di_13].get_energy());
                            trimer.get_ene_d23(dimers[di_23].get_energy());
                            // gradients
                            trimer.get_grad_m1(monomers[mono_i1].get_gradient());
                            trimer.get_grad_m2(monomers[mono_i2].get_gradient());
                            trimer.get_grad_m3(monomers[mono_i3].get_gradient());
                            trimer.get_grad_d12(dimers[di_12].get_gradient());
                            trimer.get_grad_d13(dimers[di_13].get_gradient());
                            trimer.get_grad_d23(dimers[di_23].get_gradient());
                        }
                    } else {
                        // energies
                        trimer.get_ene_m1(monomers[mono_i1].get_energy());
                        trimer.get_ene_m2(monomers[mono_i2].get_energy());
                        trimer.get_ene_m3(monomers[mono_i3].get_energy());
                        trimer.get_ene_d12(dimers[di_12].get_energy());
                        trimer.get_ene_d13(dimers[di_13].get_energy());
                        trimer.get_ene_d23(dimers[di_23].get_energy());
                        // gradients
                        trimer.get_grad_m1(monomers[mono_i1].get_gradient());
                        trimer.get_grad_m2(monomers[mono_i2].get_gradient());
                        trimer.get_grad_m3(monomers[mono_i3].get_gradient());
                        trimer.get_grad_d12(dimers[di_12].get_gradient());
                        trimer.get_grad_d13(dimers[di_13].get_gradient());
                        trimer.get_grad_d23(dimers[di_23].get_gradient());
                    }
                }
            }
            else { // pbc condition
                {// indexes
                    int i1 = trimer.index1;
                    int i2 = trimer.index2;
                    int i3 = trimer.index3;
                    int mono_i1 = i1;
                    int mono_i2 = i2;
                    int mono_i3 = i3;
                    int num_dimer = dimers.size();
                    stringstream di_ij_name;
                    stringstream di_ik_name;
                    stringstream di_jk_name;
                    di_ij_name << i1 << "_" << i2;
                    di_ik_name << i1 << "_" << i3;
                    di_jk_name << i2 << "_" << i3;
                    int di_12 = map_dimer.at(di_ij_name.str());
                    int di_13 = map_dimer.at(di_ik_name.str());
                    int di_23 = map_dimer.at(di_jk_name.str());
                    bool in_dimers_12 = (di_12 < num_dimer);
                    bool in_dimers_13 = (di_13 < num_dimer);
                    bool in_dimers_23 = (di_23 < num_dimer);
                    if (directlyDelta) {
                        bool bool_dd;
                        bool_dd = (is_element_in_vector(dd_index, monomers[mono_i1].species)
                                   and is_element_in_vector(dd_index, monomers[mono_i2].species)
                                   and is_element_in_vector(dd_index, monomers[mono_i3].species));
                        if (bool_dd) {
                            trimer.set_dd_ene_grad();
                        }
                        else {
                            // monomer --- energies
                            trimer.get_ene_m1(read_monomer_ene_not_in_monomers(mono_i1));
                            trimer.get_ene_m2(read_monomer_ene_not_in_monomers(mono_i2));
                            trimer.get_ene_m3(read_monomer_ene_not_in_monomers(mono_i3));
                            // monomer --- gradients
                            trimer.get_grad_m1(read_monomer_grad_not_in_monomers(mono_i1, trimer.get_natom_m1()));
                            trimer.get_grad_m2(read_monomer_grad_not_in_monomers(mono_i2, trimer.get_natom_m2()));
                            trimer.get_grad_m3(read_monomer_grad_not_in_monomers(mono_i3, trimer.get_natom_m3()));
                            // dimer 1
                            if (in_dimers_12) {
                                trimer.get_ene_d12(dimers[di_12].get_energy());
                                trimer.get_grad_d12(dimers[di_12].get_gradient());
                            } else {
                                int natom_dimer = trimer.ele1.size() + trimer.ele2.size();
                                trimer.get_ene_d12 (read_dimer_ene_not_in_dimers(di_12));
                                trimer.get_grad_d12(read_dimer_grad_not_in_dimers(di_12,natom_dimer));
                            }
                            // dimer 2
                            if (in_dimers_13) {
                                trimer.get_ene_d13(dimers[di_13].get_energy());
                                trimer.get_grad_d13(dimers[di_13].get_gradient());
                            } else {
                                int natom_dimer = trimer.ele1.size() + trimer.ele3.size();
                                trimer.get_ene_d13 (read_dimer_ene_not_in_dimers(di_13));
                                trimer.get_grad_d13(read_dimer_grad_not_in_dimers(di_13,natom_dimer));
                            }
                            // dimer 3
                            if (in_dimers_23) {
                                trimer.get_ene_d23(dimers[di_23].get_energy());
                                trimer.get_grad_d23(dimers[di_23].get_gradient());
                            } else {
                                int natom_dimer = trimer.ele2.size() + trimer.ele3.size();
                                trimer.get_ene_d23 (read_dimer_ene_not_in_dimers(di_23));
                                trimer.get_grad_d23(read_dimer_grad_not_in_dimers(di_23,natom_dimer));
                            }
                        }
                    }
                    else {
                        // monomer --- energies
                        trimer.get_ene_m1(read_monomer_ene_not_in_monomers(mono_i1));
                        trimer.get_ene_m2(read_monomer_ene_not_in_monomers(mono_i2));
                        trimer.get_ene_m3(read_monomer_ene_not_in_monomers(mono_i3));
                        // monomer --- gradients
                        trimer.get_grad_m1(read_monomer_grad_not_in_monomers(mono_i1, trimer.get_natom_m1()));
                        trimer.get_grad_m2(read_monomer_grad_not_in_monomers(mono_i2, trimer.get_natom_m2()));
                        trimer.get_grad_m3(read_monomer_grad_not_in_monomers(mono_i3, trimer.get_natom_m3()));
                        // dimer 1
                        if (in_dimers_12) {
                            trimer.get_ene_d12(dimers[di_12].get_energy());
                            trimer.get_grad_d12(dimers[di_12].get_gradient());
                        } else {
                            int natom_dimer = trimer.ele1.size() + trimer.ele2.size();
                            trimer.get_ene_d12 (read_dimer_ene_not_in_dimers(di_12));
                            trimer.get_grad_d12(read_dimer_grad_not_in_dimers(di_12,natom_dimer));
                        }
                        // dimer 2
                        if (in_dimers_13) {
                            trimer.get_ene_d13(dimers[di_13].get_energy());
                            trimer.get_grad_d13(dimers[di_13].get_gradient());
                        } else {
                            int natom_dimer = trimer.ele1.size() + trimer.ele3.size();
                            trimer.get_ene_d13 (read_dimer_ene_not_in_dimers(di_13));
                            trimer.get_grad_d13(read_dimer_grad_not_in_dimers(di_13,natom_dimer));
                        }
                        // dimer 3
                        if (in_dimers_23) {
                            trimer.get_ene_d23(dimers[di_23].get_energy());
                            trimer.get_grad_d23(dimers[di_23].get_gradient());
                        } else {
                            int natom_dimer = trimer.ele2.size() + trimer.ele3.size();
                            trimer.get_ene_d23 (read_dimer_ene_not_in_dimers(di_23));
                            trimer.get_grad_d23(read_dimer_grad_not_in_dimers(di_23,natom_dimer));
                        }
                    }
                }
            }
        }
    }
}

void
cal_dimer_gmm(Dimer & dimer,
        const double & prefactor){
    int na1 = dimer.ele1.size();
    int na2 = dimer.ele2.size();
    int natom=na1+na2;
    vector<vector<double > > grad_tot(natom, vector<double >(3,0));
    vector<vector<double > > grad_m1(na1,vector<double >(3,0));
    vector<vector<double > > grad_m2(na2,vector<double >(3,0));
    
    for (int j=0;j<na1;j++){
        vector<double > crd_j = dimer.coord1[j] * 1.889726; // convert unit to Bohr
        double charge_j = dimer.charge1[j];
        for (int k=0;k<na2;k++){
            vector<double > crd_k = dimer.coord2[k] * 1.889726; // convert unit to Bohr
            double charge_k = dimer.charge2[k];

            vector<double > r_jk_vec = crd_k - crd_j;
            double r_jk = sqrt(norm2(r_jk_vec));
            double r3 = r_jk * r_jk * r_jk;
            vector<double > grad = r_jk_vec * (1./r3);
            grad = grad * charge_j;
            grad = grad * charge_k;
            grad_m1[j][0] += prefactor * -1. * grad[0];
            grad_m1[j][1] += prefactor * -1. * grad[1];
            grad_m1[j][2] += prefactor * -1. * grad[2];

            grad_m2[k][0] += prefactor * grad[0];
            grad_m2[k][1] += prefactor * grad[1];
            grad_m2[k][2] += prefactor * grad[2];

            grad_tot[j][0] += prefactor * -1. * grad[0];
            grad_tot[j][1] += prefactor * -1. * grad[1];
            grad_tot[j][2] += prefactor * -1. * grad[2];
            grad_tot[k+na1][0] += prefactor * grad[0];
            grad_tot[k+na1][1] += prefactor * grad[1];
            grad_tot[k+na1][2] += prefactor * grad[2];
        }
    }
    dimer.get_gradient(grad_tot);
    dimer.get_grad_m1(grad_m1);
    dimer.get_grad_m2(grad_m2);
}

double mbe_energy(const vector<Monomer > & monomers,
                  const vector<Dimer > & dimers,
                  const double & rcut_qm,
                  const double & rcut_sr,
                  const double & rcut_lr,
                  const int & stt_index){
    bool pbc = false;
    if (stt_index >= 0){
        pbc = true;
    }
    double ene = 0.0;
    // monomer
    double e_monomer = 0.0;
    for (int i=0;i<monomers.size();i++){
        double e_mono = monomers[i].get_energy();
        e_monomer += e_mono;
    }
    // dimer
    double e_dimer_qm = 0.0;
    double e_dimer_sr = 0.0;
    double e_dimer_lr = 0.0;
    if (not pbc) {
        for (int i = 0; i < dimers.size(); i++) {
            const Dimer &dimer = dimers[i];
            if (dimer.rij < rcut_qm) {
                double e_dim = dimer.get_energy() - dimer.get_ene_m1() - dimer.get_ene_m2();
                e_dimer_qm += e_dim;
            } else if (dimer.rij < rcut_sr) {
                double e_dim = cal_dimer_emm(dimer, -1.0);
                e_dimer_sr += e_dim;
            } else if (dimer.rij < rcut_lr) {
                double e_dim = cal_dimer_emm(dimer, 1.0);
                e_dimer_lr += e_dim;
            }
        }
    }
    else { // pbc condition
        int nmol_uc = monomers.size();
        for (int i = 0; i < dimers.size(); i++) {
            const Dimer &dimer = dimers[i];
            bool both_in_uc = false;
            if ( (dimer.index1 >= stt_index and dimer.index1 < (stt_index+nmol_uc)) and
                 (dimer.index2 >= stt_index and dimer.index2 < (stt_index+nmol_uc)) ){
                both_in_uc = true;
            }
            if (dimer.rij < rcut_qm) {
                double e_dim = dimer.get_energy() - dimer.get_ene_m1() - dimer.get_ene_m2();
                if (not both_in_uc) {
                    e_dim = e_dim * 0.5;
                }
                if (abs(e_dim) > 0.002) {
                    cout << "e_dim = " << e_dim << " for dimer_" << i << ", mono is " << dimer.index1 << " "
                         << dimer.index2 << endl;
                }
                e_dimer_qm += e_dim;
            }
        }
    }
    cout << setw(18) << setprecision(9) << "e_monomer = " << e_monomer << endl;
    cout << setw(18) << setprecision(9) << "e_dimer_qm= " << e_dimer_qm << endl;
    ene = e_monomer + e_dimer_qm + e_dimer_sr + e_dimer_lr;
    return ene;
}

double cal_dimer_emm(const Dimer & dimer,
        const double & prefactor){
    double ene_clmb_intractn = 0.0;
    for (int i=0;i<dimer.ele1.size();i++){
        vector<double > crd_i = dimer.coord1[i] * 1.889726; // convert unit to Bohr
        for (int j=0;j<dimer.ele2.size();j++){
            vector<double > crd_j = dimer.coord2[j] * 1.889726; // convert unit to Bohr
            vector<double > r_ij_vec = crd_j - crd_i;
            double rij = sqrt(norm2(r_ij_vec));
            double qa = dimer.charge1[i];
            double qb = dimer.charge2[j];
            ene_clmb_intractn += qa*qb/rij;
        }
    }
    ene_clmb_intractn = prefactor * ene_clmb_intractn;
    return ene_clmb_intractn;
}

double cal_dimer_esr(const Dimer & dimer){ // useless now, may be have bugs
    double rij = dimer.rij;
    double ene_clmb_intractn = 0.0;
    for (int i=0;i<dimer.ele1.size();i++){
        for (int j=0;j<dimer.ele2.size();j++){
            double qa = dimer.charge1[i];
            double qb = dimer.charge2[j];
            ene_clmb_intractn += qa*qb/rij;
        }
    }
    ene_clmb_intractn = -1.0 * ene_clmb_intractn;
    return ene_clmb_intractn;
}

double cal_dimer_elr(const Dimer & dimer){ // useless now, may be have bugs
    double rij = dimer.rij;
    double ene_long_range = 0.0;
    for (int i=0;i<dimer.ele1.size();i++){
        for (int j=0;j<dimer.ele2.size();j++){
            double qa = dimer.charge1[i];
            double qb = dimer.charge2[j];
            ene_long_range += qa*qb/rij;
        }
    }
    return ene_long_range;
}

vector<vector<double > > mbe_gradient(
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

    vector<vector<double > > grad(ntot_atom, vector<double > (3,0));
    
    cal_mono_grad(grad,monomers,natom_set,natom_cum);
    
    cal_dimer_grad(grad,dimers,rcut_qm,rcut_sr,rcut_lr,natom_set,natom_cum,stt_index);

    return grad;
}

double mbe_energy(const vector<Monomer > & monomers,
                  const vector<Dimer > & dimers,
                  const vector<Trimer > & trimers,
                  const double & rcut_qm,
                  const double & rcut_sr,
                  const double & rcut_lr,
                  const bool & directlyDelta,
                  const vector<int> & dd_index,
                  const int & stt_index){
    // for trimers.
    bool pbc = false;
    if (stt_index >= 0){
        pbc = true;
    }
    double ene = 0.0;
    // monomer
    double e_monomer = 0.0;
    for (int i=0;i<monomers.size();i++){
        double e_mono = monomers[i].get_energy();
        e_monomer += e_mono;
    }
    // dimer
    double e_dimer_qm = 0.0;
    double e_dimer_sr = 0.0;
    double e_dimer_lr = 0.0;
    if (not pbc) {
        for (int i = 0; i < dimers.size(); i++) {
            const Dimer &dimer = dimers[i];
            if (dimer.rij < rcut_qm) {
                double e_dim = dimer.get_energy() - dimer.get_ene_m1() - dimer.get_ene_m2();
                e_dimer_qm += e_dim;
            } else if (dimer.rij < rcut_sr) {
                double e_dim = cal_dimer_emm(dimer, -1.0);
                e_dimer_sr += e_dim;
            } else if (dimer.rij < rcut_lr) {
                double e_dim = cal_dimer_emm(dimer, 1.0);
                e_dimer_lr += e_dim;
            }
        }
    }
    else { 
        int nmol_uc = monomers.size();
        for (int i = 0; i < dimers.size(); i++) {
            const Dimer &dimer = dimers[i];
            bool both_in_uc = false;
            if ( (dimer.index1 >= stt_index and dimer.index1 < (stt_index+nmol_uc)) and
                 (dimer.index2 >= stt_index and dimer.index2 < (stt_index+nmol_uc)) ){
                both_in_uc = true;
            }
            if (dimer.rij < rcut_qm) {
                double e_dim = dimer.get_energy() - dimer.get_ene_m1() - dimer.get_ene_m2();
                if (not both_in_uc) {
                    e_dim = e_dim * 0.5;
                }
                e_dimer_qm += e_dim;
            }
        }
    }

    // trimer
    double e_trimer_qm = 0.0;
    if (not pbc) {
        for (int i = 0; i < trimers.size(); i++) {
            const Trimer &trimer = trimers[i];
            if ((trimer.rij < rcut_qm) && (trimer.rik < rcut_qm) && (trimer.rjk < rcut_qm)) {
                if (not directlyDelta) {
                    double e_trim =
                            trimer.get_energy() - trimer.get_ene_d12() - trimer.get_ene_d13() - trimer.get_ene_d23() +
                            trimer.get_ene_m1() + trimer.get_ene_m2() + trimer.get_ene_m3();
                    e_trimer_qm += e_trim;
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
                    double e_trim = trimer.get_energy();
                    if (not bool_dd) { 
                        e_trim = trimer.get_energy();
                        bool bool_dd_dij = (is_element_in_vector(dd_index, monomers[mono_i1].species)
                                            and is_element_in_vector(dd_index, monomers[mono_i2].species));
                        bool bool_dd_dik = (is_element_in_vector(dd_index, monomers[mono_i1].species)
                                            and is_element_in_vector(dd_index, monomers[mono_i3].species));
                        bool bool_dd_djk = (is_element_in_vector(dd_index, monomers[mono_i2].species)
                                            and is_element_in_vector(dd_index, monomers[mono_i3].species));
                        if (bool_dd_dij) { e_trim -= trimer.get_ene_d12(); }
                        else { e_trim -= (trimer.get_ene_d12() - trimer.get_ene_m1() - trimer.get_ene_m2()); }
                        if (bool_dd_dik) { e_trim -= trimer.get_ene_d13(); }
                        else { e_trim -= (trimer.get_ene_d13() - trimer.get_ene_m1() - trimer.get_ene_m3()); }
                        if (bool_dd_djk) { e_trim -= trimer.get_ene_d23(); }
                        else { e_trim -= (trimer.get_ene_d23() - trimer.get_ene_m2() - trimer.get_ene_m3()); }
                        e_trim -= (trimer.get_ene_m1() + trimer.get_ene_m2() + trimer.get_ene_m3());
                    }
                    e_trimer_qm += e_trim;
                }
            }
        }
    }
    else {
        int nmol_uc = monomers.size();
        for (int i = 0; i < trimers.size(); i++) {
            const Trimer &trimer = trimers[i];
            if ((trimer.rij < rcut_qm) && (trimer.rik < rcut_qm) && (trimer.rjk < rcut_qm)) {
                double num_in_uc = 0.0;
                if (trimer.index1 >= stt_index and trimer.index1 < (stt_index+nmol_uc)) {
                    num_in_uc += 1.0;
                }
                if (trimer.index2 >= stt_index and trimer.index2 < (stt_index+nmol_uc)) {
                    num_in_uc += 1.0;
                }
                if (trimer.index3 >= stt_index and trimer.index3 < (stt_index+nmol_uc)) {
                    num_in_uc += 1.0;
                }
                if (not directlyDelta) {
                    double e_trim =
                            trimer.get_energy() - trimer.get_ene_d12() - trimer.get_ene_d13() - trimer.get_ene_d23() +
                            trimer.get_ene_m1() + trimer.get_ene_m2() + trimer.get_ene_m3();
                    e_trim = e_trim * num_in_uc / 3.0;
                    e_trimer_qm += e_trim;
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
                    double e_trim = trimer.get_energy();
                    if (not bool_dd) { 
                        e_trim = trimer.get_energy();
                        bool bool_dd_dij = (is_element_in_vector(dd_index, monomers[mono_i1].species)
                                            and is_element_in_vector(dd_index, monomers[mono_i2].species));
                        bool bool_dd_dik = (is_element_in_vector(dd_index, monomers[mono_i1].species)
                                            and is_element_in_vector(dd_index, monomers[mono_i3].species));
                        bool bool_dd_djk = (is_element_in_vector(dd_index, monomers[mono_i2].species)
                                            and is_element_in_vector(dd_index, monomers[mono_i3].species));
                        if (bool_dd_dij) { e_trim -= trimer.get_ene_d12(); }
                        else { e_trim -= (trimer.get_ene_d12() - trimer.get_ene_m1() - trimer.get_ene_m2()); }
                        if (bool_dd_dik) { e_trim -= trimer.get_ene_d13(); }
                        else { e_trim -= (trimer.get_ene_d13() - trimer.get_ene_m1() - trimer.get_ene_m3()); }
                        if (bool_dd_djk) { e_trim -= trimer.get_ene_d23(); }
                        else { e_trim -= (trimer.get_ene_d23() - trimer.get_ene_m2() - trimer.get_ene_m3()); }
                        e_trim -= (trimer.get_ene_m1() + trimer.get_ene_m2() + trimer.get_ene_m3());
                    }
                    e_trim = e_trim * num_in_uc / 3.0;
                    e_trimer_qm += e_trim;
                }
            }
        }
    }
    cout << setw(18) << setprecision(9) << "e_monomer = " << e_monomer << endl;
    cout << setw(18) << setprecision(9) << "e_dimer_qm= " << e_dimer_qm << endl;
    cout << setw(18) << setprecision(9) << "e_trimerqm= " << e_trimer_qm << endl;
    ene = e_monomer + e_dimer_qm + e_dimer_sr + e_dimer_lr + e_trimer_qm;
    return ene;
}

vector<vector<double > > mbe_gradient(
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
    vector<vector<double > > grad(ntot_atom, vector<double > (3,0));
    cal_mono_grad(grad,monomers,natom_set,natom_cum);
    cal_dimer_grad(grad,dimers,rcut_qm,rcut_sr,rcut_lr,natom_set,natom_cum,stt_index);
    cal_trimer_grad(grad,monomers,trimers,rcut_qm,rcut_sr,rcut_lr,natom_set,natom_cum,directlyDelta,dd_index,stt_index);

    return grad;
}

void
cal_mono_grad(vector<vector<double > > & grad,
        const vector<Monomer > & monomers,
        const map<int, int> & natom_set,
        const map<int, int> & natom_cum){
    int nmol = monomers.size();
    for (int i=0;i<nmol;i++){
        const Monomer & monomer = monomers[i];
        int index = monomer.index;
        int nstart = natom_cum.at(index);
        int natom = natom_set.at(index);
        vector<vector<double > > grad_mono = monomer.get_gradient();
        // record into gradient.
        for (int j=0;j<natom;j++){
            grad[j+nstart][0] += grad_mono[j][0];
            grad[j+nstart][1] += grad_mono[j][1];
            grad[j+nstart][2] += grad_mono[j][2];
        }
    }
}

void
cal_dimer_grad(vector<vector<double > > & grad,
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
            vector<vector<double> > gradient = dimer.get_gradient();
            vector<vector<double> > grad_m1 = dimer.get_grad_m1();
            vector<vector<double> > grad_m2 = dimer.get_grad_m2();
            
            if (monomer1_in_uc or (stt_index < 0) ) { 
                for (int j = 0; j < na1; j++) {
                    grad[j + ns1][0] += (gradient[j][0] - grad_m1[j][0]);
                    grad[j + ns1][1] += (gradient[j][1] - grad_m1[j][1]);
                    grad[j + ns1][2] += (gradient[j][2] - grad_m1[j][2]);
                }
            } else { 
                // do nothing.
            }
            if (monomer2_in_uc or (stt_index < 0) ) { 
                for (int j = 0; j < na2; j++) {
                    grad[j + ns2][0] += (gradient[j + na1][0] - grad_m2[j][0]);
                    grad[j + ns2][1] += (gradient[j + na1][1] - grad_m2[j][1]);
                    grad[j + ns2][2] += (gradient[j + na1][2] - grad_m2[j][2]);
                }
            } else{ 
                // do nothing.
            }
        }
    }
}

void
cal_trimer_grad(vector<vector<double > > & grad,
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
            vector<vector<double > > gradient = trimer.get_gradient();
            vector<vector<double > > grad_d12 = trimer.get_grad_d12();
            vector<vector<double > > grad_d13 = trimer.get_grad_d13();
            vector<vector<double > > grad_d23 = trimer.get_grad_d23();
            vector<vector<double > > grad_m1  = trimer.get_grad_m1();
            vector<vector<double > > grad_m2  = trimer.get_grad_m2();
            vector<vector<double > > grad_m3  = trimer.get_grad_m3();
            if (not directlyDelta) {
                if (monomer1_in_uc or (stt_index < 0) ) { 
                    for (int j = 0; j < na1; j++) {
                        grad[j + ns1][0] += (gradient[j][0] - grad_d12[j][0] - grad_d13[j][0] + grad_m1[j][0]);
                        grad[j + ns1][1] += (gradient[j][1] - grad_d12[j][1] - grad_d13[j][1] + grad_m1[j][1]);
                        grad[j + ns1][2] += (gradient[j][2] - grad_d12[j][2] - grad_d13[j][2] + grad_m1[j][2]);
                    }
                }
                if (monomer2_in_uc or (stt_index < 0) ) { 
                    for (int j = 0; j < na2; j++) {
                        grad[j + ns2][0] += (gradient[j + na1][0] - grad_d12[j + na1][0] - grad_d23[j][0] +
                                             grad_m2[j][0]);
                        grad[j + ns2][1] += (gradient[j + na1][1] - grad_d12[j + na1][1] - grad_d23[j][1] +
                                             grad_m2[j][1]);
                        grad[j + ns2][2] += (gradient[j + na1][2] - grad_d12[j + na1][2] - grad_d23[j][2] +
                                             grad_m2[j][2]);
                    }
                }
                if (monomer3_in_uc or (stt_index < 0) ) { 
                    for (int j = 0; j < na3; j++) {
                        grad[j + ns3][0] += (gradient[j + na1 + na2][0] - grad_d13[j + na1][0] - grad_d23[j + na2][0] +
                                             grad_m3[j][0]);
                        grad[j + ns3][1] += (gradient[j + na1 + na2][1] - grad_d13[j + na1][1] - grad_d23[j + na2][1] +
                                             grad_m3[j][1]);
                        grad[j + ns3][2] += (gradient[j + na1 + na2][2] - grad_d13[j + na1][2] - grad_d23[j + na2][2] +
                                             grad_m3[j][2]);
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
                if (bool_dd){
                    if (monomer1_in_uc or (stt_index < 0) ) { 
                        for (int j = 0; j < na1; j++) {
                            grad[j + ns1][0] += (gradient[j][0]);
                            grad[j + ns1][1] += (gradient[j][1]);
                            grad[j + ns1][2] += (gradient[j][2]);
                        }
                    }
                    if (monomer2_in_uc or (stt_index < 0) ) { 
                        for (int j = 0; j < na2; j++) {
                            grad[j + ns2][0] += (gradient[j + na1][0]);
                            grad[j + ns2][1] += (gradient[j + na1][1]);
                            grad[j + ns2][2] += (gradient[j + na1][2]);
                        }
                    }
                    if (monomer3_in_uc or (stt_index < 0) ) { 
                        for (int j = 0; j < na3; j++) {
                            grad[j + ns3][0] += (gradient[j + na1 + na2][0]);
                            grad[j + ns3][1] += (gradient[j + na1 + na2][1]);
                            grad[j + ns3][2] += (gradient[j + na1 + na2][2]);
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
                        vector<vector<double > > dd_grad_d12(na1+na2, vector<double> (3,0.0));
                        vector<vector<double > > dd_grad_d13(na1+na3, vector<double> (3,0.0));
                        vector<vector<double > > dd_grad_d23(na2+na3, vector<double> (3,0.0));
                        bool bool_dd_dij = ( is_element_in_vector(dd_index, monomers[mono_i1].species)
                                             and is_element_in_vector(dd_index, monomers[mono_i2].species) );
                        bool bool_dd_dik = ( is_element_in_vector(dd_index, monomers[mono_i1].species)
                                             and is_element_in_vector(dd_index, monomers[mono_i3].species) );
                        bool bool_dd_djk = ( is_element_in_vector(dd_index, monomers[mono_i2].species)
                                             and is_element_in_vector(dd_index, monomers[mono_i3].species) );
                        if (bool_dd_dij) {
                            for (int j = 0; j < na1+na2; j++) {
                                dd_grad_d12[j][0] = grad_d12[j][0];
                                dd_grad_d12[j][1] = grad_d12[j][1];
                                dd_grad_d12[j][2] = grad_d12[j][2];
                            }
                        } else {
                            for (int j = 0; j < na1; j++) {
                                dd_grad_d12[j][0] = grad_d12[j][0]-grad_m1[j][0];
                                dd_grad_d12[j][1] = grad_d12[j][1]-grad_m1[j][1];
                                dd_grad_d12[j][2] = grad_d12[j][2]-grad_m1[j][2];
                            }
                            for (int j = 0; j < na2; j++) {
                                dd_grad_d12[j + na1][0] = grad_d12[j + na1][0]-grad_m2[j][0];
                                dd_grad_d12[j + na1][1] = grad_d12[j + na1][1]-grad_m2[j][1];
                                dd_grad_d12[j + na1][2] = grad_d12[j + na1][2]-grad_m2[j][2];
                            }
                        }
                        if (bool_dd_dik) {
                            for (int j = 0; j < na1+na3; j++){
                                dd_grad_d13[j][0] = grad_d13[j][0];
                                dd_grad_d13[j][1] = grad_d13[j][1];
                                dd_grad_d13[j][2] = grad_d13[j][2];
                            }
                        } else {
                            for (int j = 0; j < na1; j++) {
                                dd_grad_d13[j][0] = grad_d13[j][0]-grad_m1[j][0];
                                dd_grad_d13[j][1] = grad_d13[j][1]-grad_m1[j][1];
                                dd_grad_d13[j][2] = grad_d13[j][2]-grad_m1[j][2];
                            }
                            for (int j = 0; j < na3; j++) {
                                dd_grad_d13[j + na1][0] = grad_d13[j + na1][0]-grad_m3[j][0];
                                dd_grad_d13[j + na1][1] = grad_d13[j + na1][1]-grad_m3[j][1];
                                dd_grad_d13[j + na1][2] = grad_d13[j + na1][2]-grad_m3[j][2];
                            }
                        }
                        if (bool_dd_djk) {
                            for (int j = 0; j < na2+na3; j++){
                                dd_grad_d23[j][0] = grad_d23[j][0];
                                dd_grad_d23[j][1] = grad_d23[j][1];
                                dd_grad_d23[j][2] = grad_d23[j][2];
                            }
                        } else {
                            for (int j = 0; j < na2; j++) {
                                dd_grad_d23[j][0] = grad_d23[j][0]-grad_m2[j][0];
                                dd_grad_d23[j][1] = grad_d23[j][0]-grad_m2[j][1];
                                dd_grad_d23[j][2] = grad_d23[j][0]-grad_m2[j][2];
                            }
                            for (int j = 0; j < na3; j++) {
                                dd_grad_d23[j + na2][0] = grad_d23[j + na2][0]-grad_m3[j][0];
                                dd_grad_d23[j + na2][1] = grad_d23[j + na2][0]-grad_m3[j][0];
                                dd_grad_d23[j + na2][2] = grad_d23[j + na2][0]-grad_m3[j][0];
                            }
                        }
                        // calculate total gradient
                        if (monomer1_in_uc or (stt_index < 0) ) { 
                            for (int j = 0; j < na1; j++) {
                                grad[j + ns1][0] += (gradient[j][0] - dd_grad_d12[j][0] - dd_grad_d13[j][0] -
                                                     grad_m1[j][0]);
                                grad[j + ns1][1] += (gradient[j][1] - dd_grad_d12[j][1] - dd_grad_d13[j][1] -
                                                     grad_m1[j][1]);
                                grad[j + ns1][2] += (gradient[j][2] - dd_grad_d12[j][2] - dd_grad_d13[j][2] -
                                                     grad_m1[j][2]);
                            }
                        }
                        if (monomer2_in_uc or (stt_index < 0) ) { 
                            for (int j = 0; j < na2; j++) {
                                grad[j + ns2][0] += (gradient[j + na1][0] - dd_grad_d12[j + na1][0] -
                                                     dd_grad_d23[j][0] - grad_m2[j][0]);
                                grad[j + ns2][1] += (gradient[j + na1][1] - dd_grad_d12[j + na1][1] -
                                                     dd_grad_d23[j][1] - grad_m2[j][1]);
                                grad[j + ns2][2] += (gradient[j + na1][2] - dd_grad_d12[j + na1][2] -
                                                     dd_grad_d23[j][2] - grad_m2[j][2]);
                            }
                        }
                        if (monomer3_in_uc or (stt_index < 0) ) { 
                            for (int j = 0; j < na3; j++) {
                                grad[j + ns3][0] += (gradient[j + na1 + na2][0] - dd_grad_d13[j + na1][0] -
                                                     dd_grad_d23[j + na2][0] -
                                                     grad_m3[j][0]);
                                grad[j + ns3][1] += (gradient[j + na1 + na2][1] - dd_grad_d13[j + na1][1] -
                                                     dd_grad_d23[j + na2][1] -
                                                     grad_m3[j][1]);
                                grad[j + ns3][2] += (gradient[j + na1 + na2][2] - dd_grad_d13[j + na1][2] -
                                                     dd_grad_d23[j + na2][2] -
                                                     grad_m3[j][2]);
                            }
                        }
                    }
                }
            }
        }
    }
}

vector<vector<double > >
read_mono_grad_tmp(const stringstream & gradname, const int & natom){
    ifstream gradfile;
    gradfile.open(gradname.str());
    if (!gradfile.is_open()){
        cout << "Could not open the file " << gradname.str() << endl;
        cout << "Program terminating.\n";
        exit(EXIT_FAILURE);
    }
    string line; 
    vector<vector<double > > grad;
    grad.clear();
    for (int ii=0;ii<natom;ii++){
        getline(gradfile,line);
        vector<string > data = split(line," ");
        vector<double > sg;
        sg.push_back(atof(data[0].c_str()));
        sg.push_back(atof(data[1].c_str()));
        sg.push_back(atof(data[2].c_str()));
        grad.push_back(sg);
    }
    gradfile.close();
    return grad;
}

double
read_mono_ene_tmp(const stringstream & enename) {
    ifstream enefile;
    enefile.open(enename.str());
    if (!enefile.is_open()) { // failed to open file
        cout << "Could not open the file " << enename.str() << endl;
        cout << "Program terminating.\n";
        exit(EXIT_FAILURE);
    }
    string line; 
    getline(enefile,line);
    vector<string> data = split(line, " ");
    double ene;
    istringstream str_ene(data[0]);
    str_ene >> ene;
    enefile.close();
    return ene;
}

vector<vector<double > >
read_mono_grad_from_dimer(const stringstream & gradname,
        const int & natom,
        const int & start_index,
        const int & final_index) {
    ifstream gradfile;
    gradfile.open(gradname.str());
    if (!gradfile.is_open()){
        cout << "Could not open the file " << gradname.str() << endl;
        cout << "Program terminating.\n";
        exit(EXIT_FAILURE);
    }
    string line; 
    vector<vector<double > > grad;
    grad.clear();
    for (int ii=0;ii<natom;ii++){
        getline(gradfile,line);
        if ((ii >= start_index) and (ii<final_index)) {
            vector<string> data = split(line, " ");
            vector<double> sg;
            sg.push_back(atof(data[0].c_str()));
            sg.push_back(atof(data[1].c_str()));
            sg.push_back(atof(data[2].c_str()));
            grad.push_back(sg);
        }
    }
    gradfile.close();
    return grad;
}

void
read_only_e_mono(vector<Monomer > & monomers){
    int n_mono = monomers.size();
    cout << "\tReading energy file for monomers...\n";
    for (int i=0; i<n_mono;i++){
        Monomer & monomer = monomers[i];
        stringstream enename;

        enename << "monomer_" << i << ".ene";
        {
            ifstream enefile;
            enefile.open(enename.str());
            if (!enefile.is_open()) { // failed to open file
                cout << "Could not open the file " << enename.str() << endl;
                cout << "Program terminating.\n";
                exit(EXIT_FAILURE);
            }
            string line; 
            getline(enefile,line);
            vector<string> data = split(line, " ");
            double ene;
            istringstream str_ene(data[0]);
            str_ene >> ene;
            monomer.get_energy(ene);
            enefile.close();
        }
    }
    // cout << "\tDone.\n";
}

void
read_only_e_dimer(vector<Dimer > & dimers,
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
                  const bool & mono_bgc){
    int n_dimer = dimers.size();
    cout << "\tReading energy file for dimers...\n";
    for (int i=0;i<n_dimer;i++) {
        Dimer &dimer = dimers[i];
        if (dimer.rij < rcut_qm) {
            stringstream enename;

            enename << "dimer_" << i << ".ene";
            {
                ifstream enefile;
                enefile.open(enename.str());
                if (!enefile.is_open()) { // failed to open file
                    cout << "Could not open the file " << enename.str() << endl;
                    cout << "Program terminating.\n";
                    exit(EXIT_FAILURE);
                }
                string line; 
                getline(enefile, line);
                vector<string> data = split(line, " ");
                double ene;
                istringstream str_ene(data[0]);
                str_ene >> ene;
                dimer.get_energy(ene);
                enefile.close();
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
                        dimer.set_dd_ene_grad();
                        already_done_mono = true;
                    }
                }
                if (not already_done_mono) {
                    if (partial_nobgc) {
                        // monomer_1
                        dimer.get_ene_m1(monomers[mono_i1].get_energy());
                        if (not is_element_in_vector(ro_index,monomers[mono_i1].species)){
                            stringstream enename_nobgc;
                            enename_nobgc  << "dimer_" << i << "_nobgc.ene";
                            double other_ene;
                            if (is_element_in_vector(ro_index, monomers[mono_i2].species)){
                                stringstream enename_m2_nobgc;
                                enename_m2_nobgc << "monomer_" << mono_i2 << "_nobgc.ene";
                                other_ene = read_mono_ene_tmp(enename_m2_nobgc);
                                dimer.get_energy(other_ene+dimer.get_ene_m1());
                            }
                        }
                        // monomer_2
                        dimer.get_ene_m2(monomers[mono_i2].get_energy());
                        if (not is_element_in_vector(ro_index,monomers[mono_i2].species)){
                            stringstream enename_nobgc;
                            enename_nobgc  << "dimer_" << i << "_nobgc.ene";
                            int natom = dimer.ele.size();
                            double other_ene;
                            if (is_element_in_vector(ro_index, monomers[mono_i1].species)){
                                stringstream enename_m1_nobgc;
                                enename_m1_nobgc << "monomer_" << mono_i1 << "_nobgc.ene";
                                other_ene = read_mono_ene_tmp(enename_m1_nobgc);
                                dimer.get_energy(other_ene+dimer.get_ene_m2());
                            }
                        }
                    }
                    else if (mono_bgc) {
                        // monomer_1
                        if (not is_element_in_vector(ro_index,monomers[mono_i1].species)){
                            stringstream enename_bgc;
                            enename_bgc  << "monomer_" << mono_i1 << "_bgc.ene";
                            dimer.get_ene_m1(read_mono_ene_tmp(enename_bgc));
                        } else {
                            dimer.get_ene_m1(monomers[mono_i1].get_energy());
                        }
                        // monomer_2
                        if (not is_element_in_vector(ro_index,monomers[mono_i2].species)){
                            stringstream enename_bgc;
                            enename_bgc  << "monomer_" << mono_i2 << "_bgc.ene";
                            dimer.get_ene_m2(read_mono_ene_tmp(enename_bgc));
                        } else {
                            dimer.get_ene_m2(monomers[mono_i2].get_energy());
                        }
                    }
                    else {
                        // energies
                        dimer.get_ene_m1(monomers[mono_i1].get_energy());
                        dimer.get_ene_m2(monomers[mono_i2].get_energy());
                    }
                }
            }
        }
    }
}

void
read_only_e_trimer(vector<Trimer > & trimers,
                   const vector<Dimer > & dimers,
                   const vector<Monomer > & monomers,
                   const double & rcut_qm,
                   const double & rcut_sr,
                   const double & rcut_lr,
                   const map<string, int> & map_monomer,
                   const map<string, int> & map_dimer,
                   const map<string, int> & map_trimer,
                   const bool & directlyDelta,
                   const vector<int> & dd_index){
    int n_trimer = trimers.size();
    cout << "\tReading energy file for trimers...\n";
    for (int i=0;i<n_trimer;i++){
        Trimer & trimer = trimers[i];
        if ( (trimer.rij < rcut_qm) && (trimer.rik < rcut_qm) && (trimer.rjk < rcut_qm) ){
            stringstream enename;

            enename << "trimer_" << i << ".ene";
            {
                ifstream enefile;
                enefile.open(enename.str());
                if (!enefile.is_open()) { // failed to open file
                    cout << "Could not open the file " << enename.str() << endl;
                    cout << "Program terminating.\n";
                    exit(EXIT_FAILURE);
                }
                string line; 
                getline(enefile, line);
                vector<string> data = split(line, " ");
                double ene;
                istringstream str_ene(data[0]);
                str_ene >> ene;
                trimer.get_energy(ene);
                enefile.close();
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
                        trimer.set_dd_ene_grad();
                    }
                    else {
                        // energies
                        trimer.get_ene_m1(monomers[mono_i1].get_energy());
                        trimer.get_ene_m2(monomers[mono_i2].get_energy());
                        trimer.get_ene_m3(monomers[mono_i3].get_energy());
                        trimer.get_ene_d12(dimers[di_12].get_energy());
                        trimer.get_ene_d13(dimers[di_13].get_energy());
                        trimer.get_ene_d23(dimers[di_23].get_energy());
                    }
                }
                else {
                    // energies
                    trimer.get_ene_m1(monomers[mono_i1].get_energy());
                    trimer.get_ene_m2(monomers[mono_i2].get_energy());
                    trimer.get_ene_m3(monomers[mono_i3].get_energy());
                    trimer.get_ene_d12(dimers[di_12].get_energy());
                    trimer.get_ene_d13(dimers[di_13].get_energy());
                    trimer.get_ene_d23(dimers[di_23].get_energy());
                }
                {// calculate the interaction energy
                    double interaction_ene;
                    interaction_ene = trimer.get_energy() + trimer.get_ene_m1() + trimer.get_ene_m2() + trimer.get_ene_m3()
                                      - trimer.get_ene_d12() - trimer.get_ene_d13() - trimer.get_ene_d23();
                    stringstream tri_int_ene_fn;
                    tri_int_ene_fn << "trimer_" << i << ".int";
                    ofstream tri_int_f;
                    tri_int_f.open(tri_int_ene_fn.str());
                    if (tri_int_f.is_open()) {
                        tri_int_f << setw(18) << setprecision(9) << interaction_ene;
                        tri_int_f.close();
                    }
                }
            }
        }
    }
}


vector<vector<double > >
read_dimer_grad_not_in_dimers(const int & index, const int & natom){
    stringstream dimer_gradname;
    dimer_gradname << "dimer_" << index << ".grad";
    ifstream gradfile;
    gradfile.open(dimer_gradname.str());
    if (!gradfile.is_open()) {
        cout << "Could not open the file " << dimer_gradname.str() << endl;
        cout << "Program terminating.\n";
        exit(EXIT_FAILURE);
    }
    string line; 
    vector<vector<double> > grad;
    grad.clear();
    grad.reserve(natom);
    for (int ii = 0; ii < natom; ii++) {
        getline(gradfile, line);
        vector<string> data = split(line, " ");
        vector<double> sg;
        sg.push_back(atof(data[0].c_str()));
        sg.push_back(atof(data[1].c_str()));
        sg.push_back(atof(data[2].c_str()));
        grad.push_back(sg);
    }
    gradfile.close();
    return grad;
}

double
read_dimer_ene_not_in_dimers(const int & index){
    
    stringstream dimer_enename;
    dimer_enename << "dimer_" << index << ".ene";
    ifstream enefile;
    enefile.open(dimer_enename.str());
    if (!enefile.is_open()) { // failed to open file
        cout << "Could not open the file " << dimer_enename.str() << endl;
        cout << "Program terminating.\n";
        exit(EXIT_FAILURE);
    }
    string line; 
    getline(enefile, line);
    vector<string> data = split(line, " ");
    double ene;
    istringstream str_ene(data[0]);
    str_ene >> ene;
    enefile.close();
    return ene;
}

vector<vector<double > >
read_monomer_grad_not_in_monomers(const int & index, const int & natom){
    stringstream monomer_gradname;
    monomer_gradname << "monomer_" << index << ".grad";
    ifstream gradfile;
    gradfile.open(monomer_gradname.str());
    if (!gradfile.is_open()) {
        cout << "Could not open the file " << monomer_gradname.str() << endl;
        cout << "Program terminating.\n";
        exit(EXIT_FAILURE);
    }
    string line; 
    vector<vector<double> > grad;
    grad.clear();
    grad.reserve(natom);
    for (int ii = 0; ii < natom; ii++) {
        getline(gradfile, line);
        vector<string> data = split(line, " ");
        vector<double> sg;
        sg.push_back(atof(data[0].c_str()));
        sg.push_back(atof(data[1].c_str()));
        sg.push_back(atof(data[2].c_str()));
        grad.push_back(sg);
    }
    gradfile.close();
    return grad;
}

double
read_monomer_ene_not_in_monomers(const int & index){
    
    stringstream monomer_enename;
    monomer_enename << "monomer_" << index << ".ene";
    ifstream enefile;
    enefile.open(monomer_enename.str());
    if (!enefile.is_open()) { // failed to open file
        cout << "Could not open the file " << monomer_enename.str() << endl;
        cout << "Program terminating.\n";
        exit(EXIT_FAILURE);
    }
    string line; 
    getline(enefile, line);
    vector<string> data = split(line, " ");
    double ene;
    istringstream str_ene(data[0]);
    str_ene >> ene;
    enefile.close();
    return ene;
}

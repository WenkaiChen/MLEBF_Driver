#include "../include/qm_sp_driver.h"

using namespace std;
void
qm_sp_driver(const int & method,
        const vector<Monomer > & monomers,
        const vector<Dimer > & dimers,
        const vector<Trimer > & trimers,
        map<string, int> & map_monomer,
        map<string, int> & map_dimer,
        map<string, int> & map_trimer,
        const double & rcut_qm,
        const bool & reduce_order,
        const vector<int> &  ro_index,
        const bool & partial_nobgc,
        const bool & mono_bgc,
        int order,
        bool pbc
        ){
    if (method != 9999) { // just write the coord file, then processed by external script for qm/dnn calculation.
        if (order >= 2) {
            cout << "\nOnly xyz file and bgc file is written.\n"
                 << "You may need to do qm calculation with external scripts\n";
            int n_mono = monomers.size();
            int n_dimer = dimers.size();
            cout << "\tWriting Monomer files (xyz and bgc file)...\n";
            write_monomer(monomers, n_mono, reduce_order, ro_index, partial_nobgc, mono_bgc);
            cout << "\tDone.\n";
            cout << "\tWriting Dimer files (xyz and bgc file)...\n";
            write_dimer(dimers, n_dimer, rcut_qm, reduce_order, ro_index, partial_nobgc, pbc);
            cout << "\tDone.\n";
        }
        if (order > 2) {
            int n_trimer = trimers.size();
            cout << "\tWriting Trimer files (xyz and bgc file)...\n";
            write_trimer(trimers, n_trimer, rcut_qm, pbc, map_monomer, map_dimer, map_trimer);
        }
        cout << "QM_SP_FILE has been writen.\n";
    }
}

void
write_monomer(const vector<Monomer > & monomers,
        const int & n_mono,
        const bool & reduce_order,
        const vector<int> & ro_index,
        const bool & partial_nobgc,
        const bool & mono_bgc){
    for (int i=0; i<n_mono; i++){
        const Monomer & monomer = monomers[i];
        stringstream xyzname;
        stringstream bgcname;
        xyzname << "monomer_" << monomer.index << ".xyz";
        bgcname << "monomer_" << monomer.index << ".bgc";
        
        ofstream xyzfile;
        xyzfile.open(xyzname.str());
        xyzfile << fixed << right;
        xyzfile << setw(3) << monomer.ele.size() << "\n";
        
        xyzfile << left << setw(5) << monomer.tot_charge
                << setw(3) << monomer.tot_spin << right << "\n";
        for (int ii=0; ii<monomer.ele.size(); ii++){
            xyzfile << setw(3) << atmnb2ele.at(monomer.ele[ii])
                    << setw(18) << setprecision(9) << monomer.coord[ii][0]
                    << setw(18) << setprecision(9) << monomer.coord[ii][1]
                    << setw(18) << setprecision(9) << monomer.coord[ii][2] << "\n";
        }
        xyzfile.close();
        
        ofstream bgcfile;
        bgcfile.open(bgcname.str());
        bgcfile << fixed << right;
        for (int ii=0; ii<monomer.bgc_chg.size();ii++){
            bgcfile << setw(18) << setprecision(9) << monomer.bgc_crd[ii][0]
                    << setw(18) << setprecision(9) << monomer.bgc_crd[ii][1]
                    << setw(18) << setprecision(9) << monomer.bgc_crd[ii][2]
                    << setw(20) << setprecision(5) << monomer.bgc_chg[ii] << "\n";
        }
        bgcfile.close();
        if (partial_nobgc and (is_element_in_vector(ro_index, monomer.species))){
            stringstream xyzname2;
            stringstream bgcname2;
            xyzname2 << "monomer_" << monomer.index << "_nobgc.xyz";
            bgcname2 << "monomer_" << monomer.index << "_nobgc.bgc";
            
            ofstream xyzfile2;
            xyzfile2.open(xyzname2.str());
            xyzfile2 << fixed << right;
            xyzfile2 << setw(3) << monomer.ele.size() << "\n";
            
            xyzfile2 << left << setw(5) << monomer.tot_charge
                     << setw(3) << monomer.tot_spin << right << "\n";
            for (int ii=0; ii<monomer.ele.size(); ii++){
                xyzfile2 << setw(3) << atmnb2ele.at(monomer.ele[ii])
                         << setw(18) << setprecision(9) << monomer.coord[ii][0]
                         << setw(18) << setprecision(9) << monomer.coord[ii][1]
                         << setw(18) << setprecision(9) << monomer.coord[ii][2] << "\n";
            }
            xyzfile2.close();
            
            ofstream bgcfile2;
            bgcfile2.open(bgcname2.str());
            bgcfile2 << fixed << right;
            for (int ii=0; ii<monomer.bgc_chg.size();ii++){
                bgcfile2 << setw(18) << setprecision(9) << monomer.bgc_crd[ii][0]
                         << setw(18) << setprecision(9) << monomer.bgc_crd[ii][1]
                         << setw(18) << setprecision(9) << monomer.bgc_crd[ii][2]
                         << setw(20) << setprecision(5) << monomer.bgc_chg[ii] << "\n";
            }
            bgcfile2.close();
        }
        if (mono_bgc and (not is_element_in_vector(ro_index, monomer.species))){
            stringstream xyzname2;
            stringstream bgcname2;
            xyzname2 << "monomer_" << monomer.index << "_bgc.xyz";
            bgcname2 << "monomer_" << monomer.index << "_bgc.bgc";
            
            ofstream xyzfile2;
            xyzfile2.open(xyzname2.str());
            xyzfile2 << fixed << right;
            xyzfile2 << setw(3) << monomer.ele.size() << "\n";
            
            xyzfile2 << left << setw(5) << monomer.tot_charge
                     << setw(3) << monomer.tot_spin << right << "\n";
            for (int ii=0; ii<monomer.ele.size(); ii++){
                xyzfile2 << setw(3) << atmnb2ele.at(monomer.ele[ii])
                         << setw(18) << setprecision(9) << monomer.coord[ii][0]
                         << setw(18) << setprecision(9) << monomer.coord[ii][1]
                         << setw(18) << setprecision(9) << monomer.coord[ii][2] << "\n";
            }
            xyzfile2.close();
            
            ofstream bgcfile2;
            bgcfile2.open(bgcname2.str());
            bgcfile2 << fixed << right;
            for (int ii=0; ii<monomer.bgc_chg.size();ii++){
                bgcfile2 << setw(18) << setprecision(9) << monomer.bgc_crd[ii][0]
                         << setw(18) << setprecision(9) << monomer.bgc_crd[ii][1]
                         << setw(18) << setprecision(9) << monomer.bgc_crd[ii][2]
                         << setw(20) << setprecision(5) << monomer.bgc_chg[ii] << "\n";
            }
            bgcfile2.close();
        }
    }
}

void
write_dimer(const vector<Dimer > & dimers,
        const int & n_dimer,
        const double & rcut_qm,
        const bool & reduce_order,
        const vector<int> &  ro_index,
        const bool & partial_nobgc,
        const bool & pbc){
    for (int i=0; i<n_dimer; i++){
        Dimer dimer = dimers[i];
        if (dimer.rij < rcut_qm) {
            stringstream xyzname;
            stringstream bgcname;
            xyzname << "dimer_" << i << ".xyz";
            bgcname << "dimer_" << i << ".bgc";
            
            ofstream xyzfile;
            xyzfile.open(xyzname.str());
            xyzfile << fixed << right;
            xyzfile << setw(3) << dimer.ele.size() << "\n";
            
            xyzfile << left << setw(5) << dimer.tot_charge
                    << setw(3) << dimer.tot_spin << right << "\n";
            for (int ii = 0; ii < dimer.ele.size(); ii++) {
                xyzfile << setw(3) << atmnb2ele.at(dimer.ele[ii])
                        << setw(18) << setprecision(9) << dimer.coord[ii][0]
                        << setw(18) << setprecision(9) << dimer.coord[ii][1]
                        << setw(18) << setprecision(9) << dimer.coord[ii][2] << "\n";
            }
            xyzfile.close();
            
            ofstream bgcfile;
            bgcfile.open(bgcname.str());
            bgcfile << fixed << right;
            for (int ii = 0; ii < dimer.bgc_chg.size(); ii++) {
                bgcfile << setw(18) << setprecision(9) << dimer.bgc_crd[ii][0]
                        << setw(18) << setprecision(9) << dimer.bgc_crd[ii][1]
                        << setw(18) << setprecision(9) << dimer.bgc_crd[ii][2]
                        << setw(20) << setprecision(5) << dimer.bgc_chg[ii] << "\n";
            }
            bgcfile.close();
            if (pbc) {
                {
                    stringstream xyzname;
                    stringstream bgcname;
                    xyzname << "monomer_" << dimer.index1 << ".xyz";
                    bgcname << "monomer_" << dimer.index1 << ".bgc";
                    if(access(xyzname.str().c_str(), F_OK)!=0) { 
                        Monomer m1 = dimer.get_m1();
                        
                        ofstream xyzfile;
                        xyzfile.open(xyzname.str());
                        xyzfile << fixed << right;
                        xyzfile << setw(3) << m1.ele.size() << "\n";
                        
                        xyzfile << left << setw(5) << m1.tot_charge
                                << setw(3) << m1.tot_spin << right << "\n";
                        for (int ii = 0; ii < m1.ele.size(); ii++) {
                            xyzfile << setw(3) << atmnb2ele.at(m1.ele[ii])
                                    << setw(18) << setprecision(9) << m1.coord[ii][0]
                                    << setw(18) << setprecision(9) << m1.coord[ii][1]
                                    << setw(18) << setprecision(9) << m1.coord[ii][2] << "\n";
                        }
                        xyzfile.close();
                        
                        ofstream bgcfile;
                        bgcfile.open(bgcname.str());
                        bgcfile << fixed << right;
                        for (int ii = 0; ii < m1.bgc_chg.size(); ii++) {
                            bgcfile << setw(18) << setprecision(9) << m1.bgc_crd[ii][0]
                                    << setw(18) << setprecision(9) << m1.bgc_crd[ii][1]
                                    << setw(18) << setprecision(9) << m1.bgc_crd[ii][2]
                                    << setw(20) << setprecision(5) << m1.bgc_chg[ii] << "\n";
                        }
                        bgcfile.close();
                    }
                }
                {
                    stringstream xyzname;
                    stringstream bgcname;
                    xyzname << "monomer_" << dimer.index2 << ".xyz";
                    bgcname << "monomer_" << dimer.index2 << ".bgc";
                    if(access(xyzname.str().c_str(), F_OK)!=0) { 
                        Monomer m2 = dimer.get_m2();
                        
                        ofstream xyzfile;
                        xyzfile.open(xyzname.str());
                        xyzfile << fixed << right;
                        xyzfile << setw(3) << m2.ele.size() << "\n";
                        
                        xyzfile << left << setw(5) << m2.tot_charge
                                << setw(3) << m2.tot_spin << right << "\n";
                        for (int ii = 0; ii < m2.ele.size(); ii++) {
                            xyzfile << setw(3) << atmnb2ele.at(m2.ele[ii])
                                    << setw(18) << setprecision(9) << m2.coord[ii][0]
                                    << setw(18) << setprecision(9) << m2.coord[ii][1]
                                    << setw(18) << setprecision(9) << m2.coord[ii][2] << "\n";
                        }
                        xyzfile.close();
                        
                        ofstream bgcfile;
                        bgcfile.open(bgcname.str());
                        bgcfile << fixed << right;
                        for (int ii = 0; ii < m2.bgc_chg.size(); ii++) {
                            bgcfile << setw(18) << setprecision(9) << m2.bgc_crd[ii][0]
                                    << setw(18) << setprecision(9) << m2.bgc_crd[ii][1]
                                    << setw(18) << setprecision(9) << m2.bgc_crd[ii][2]
                                    << setw(20) << setprecision(5) << m2.bgc_chg[ii] << "\n";
                        }
                        bgcfile.close();
                    }
                }
            }
            if (partial_nobgc and (   is_element_in_vector(ro_index, dimer.species1)
                                   or is_element_in_vector(ro_index, dimer.species2))){
                stringstream xyzname2;
                stringstream bgcname2;
                xyzname2 << "dimer_" << i << "_nobgc.xyz";
                bgcname2 << "dimer_" << i << "_nobgc.bgc";
                
                ofstream xyzfile2;
                xyzfile2.open(xyzname2.str());
                xyzfile2 << fixed << right;
                xyzfile2 << setw(3) << dimer.ele.size() << "\n";
                
                xyzfile2 << left << setw(5) << dimer.tot_charge
                         << setw(3) << dimer.tot_spin << right << "\n";
                for (int ii=0; ii<dimer.ele.size(); ii++){
                    xyzfile2 << setw(3) << atmnb2ele.at(dimer.ele[ii])
                             << setw(18) << setprecision(9) << dimer.coord[ii][0]
                             << setw(18) << setprecision(9) << dimer.coord[ii][1]
                             << setw(18) << setprecision(9) << dimer.coord[ii][2] << "\n";
                }
                xyzfile2.close();
                
                ofstream bgcfile2;
                bgcfile2.open(bgcname2.str());
                bgcfile2 << fixed << right;
                for (int ii=0; ii<dimer.bgc_chg.size();ii++){
                    bgcfile2 << setw(18) << setprecision(9) << dimer.bgc_crd[ii][0]
                             << setw(18) << setprecision(9) << dimer.bgc_crd[ii][1]
                             << setw(18) << setprecision(9) << dimer.bgc_crd[ii][2]
                             << setw(20) << setprecision(5) << dimer.bgc_chg[ii] << "\n";
                }
                bgcfile2.close();
            }
        }
    }
}

void
write_trimer(const vector<Trimer > & trimers,
             const int & n_trimer,
             const double & rcut_qm,
             const bool & pbc,
             map<string, int> & map_monomer,
             map<string, int> & map_dimer,
             map<string, int> & map_trimer){
    // Caution:: Now rcut_sr and rcut_lr is not implenmented.
    int count_dimer = map_dimer.size();
    for (int i=0; i<n_trimer; i++){
        Trimer trimer = trimers[i];
        if ( (trimer.rij < rcut_qm) && (trimer.rik < rcut_qm) && (trimer.rij < rcut_qm) ){
            stringstream trim_xyzname;
            stringstream trim_bgcname;
            trim_xyzname << "trimer_" << i << ".xyz";
            trim_bgcname << "trimer_" << i << ".bgc";
            
            ofstream trim_xyzfile;
            trim_xyzfile.open(trim_xyzname.str());
            trim_xyzfile << fixed << right;
            trim_xyzfile << setw(3) << trimer.ele.size() << "\n";
            
            trim_xyzfile << left << setw(5) << trimer.tot_charge
                         << setw(3) << trimer.tot_spin << right << "\n";
            for (int ii = 0; ii < trimer.ele.size(); ii++) {
                trim_xyzfile << setw(3) << atmnb2ele.at(trimer.ele[ii])
                             << setw(18) << setprecision(9) << trimer.coord[ii][0]
                             << setw(18) << setprecision(9) << trimer.coord[ii][1]
                             << setw(18) << setprecision(9) << trimer.coord[ii][2] << "\n";
            }
            trim_xyzfile.close();
            
            ofstream trim_bgcfile;
            trim_bgcfile.open(trim_bgcname.str());
            trim_bgcfile << fixed << right;
            for (int ii = 0; ii < trimer.bgc_chg.size(); ii++) {
                trim_bgcfile << setw(18) << setprecision(9) << trimer.bgc_crd[ii][0]
                             << setw(18) << setprecision(9) << trimer.bgc_crd[ii][1]
                             << setw(18) << setprecision(9) << trimer.bgc_crd[ii][2]
                             << setw(20) << setprecision(5) << trimer.bgc_chg[ii] << "\n";
            }
            trim_bgcfile.close();
            if (pbc) { // for pbc condition
                {
                    stringstream dimer_xyzname;
                    stringstream dimer_bgcname;
                    stringstream dimer_indexes;
                    dimer_indexes << trimer.index1 << "_" << trimer.index2;
                    auto findkey = map_dimer.find(dimer_indexes.str());
                    if (findkey != map_dimer.end()) { 
                        dimer_xyzname << "dimer_" << findkey->second << ".xyz";
                        dimer_bgcname << "dimer_" << findkey->second << ".bgc";
                    }
                    else{
                        map_dimer[dimer_indexes.str()] = count_dimer;
                        dimer_xyzname << "dimer_" << count_dimer << ".xyz";
                        dimer_bgcname << "dimer_" << count_dimer << ".bgc";
                        count_dimer = count_dimer+1;
                    }
                    if(access(dimer_xyzname.str().c_str(), F_OK)!=0) { 
                        int natom_dimer = trimer.get_natom_m1() + trimer.get_natom_m2();
                        int tc_dimer = trimer.tc1 + trimer.tc2;
                        int ts_dimer = 2 * (tc_dimer) + 1;
                        if (trimer.ts2 != 1 or trimer.ts1 != 1){
                            ts_dimer = max(trimer.ts1,trimer.ts2);
                        }
                        ofstream dimer_xyzfile;
                        dimer_xyzfile.open(dimer_xyzname.str());
                        dimer_xyzfile << fixed << right;
                        dimer_xyzfile << setw(3) << natom_dimer << "\n";
                        
                        dimer_xyzfile << left << setw(5) << tc_dimer
                                      << setw(3) << ts_dimer << right << "\n";
                        for (int ii = 0; ii < trimer.ele1.size(); ii++) {
                            dimer_xyzfile << setw(3) << atmnb2ele.at(trimer.ele1[ii])
                                          << setw(18) << setprecision(9) << trimer.coord1[ii][0]
                                          << setw(18) << setprecision(9) << trimer.coord1[ii][1]
                                          << setw(18) << setprecision(9) << trimer.coord1[ii][2] << "\n";
                        }
                        for (int ii = 0; ii < trimer.ele2.size(); ii++) {
                            dimer_xyzfile << setw(3) << atmnb2ele.at(trimer.ele2[ii])
                                          << setw(18) << setprecision(9) << trimer.coord2[ii][0]
                                          << setw(18) << setprecision(9) << trimer.coord2[ii][1]
                                          << setw(18) << setprecision(9) << trimer.coord2[ii][2] << "\n";
                        }
                        dimer_xyzfile.close();
                        
                        ofstream dimer_bgcfile;
                        dimer_bgcfile.open(dimer_bgcname.str());
                        dimer_bgcfile << fixed << right;
                        for (int ii = 0; ii < trimer.bgc_chg.size(); ii++) {
                            dimer_bgcfile << setw(18) << setprecision(9) << trimer.bgc_crd[ii][0]
                                          << setw(18) << setprecision(9) << trimer.bgc_crd[ii][1]
                                          << setw(18) << setprecision(9) << trimer.bgc_crd[ii][2]
                                          << setw(20) << setprecision(5) << trimer.bgc_chg[ii] << "\n";
                        }
                        for (int ii = 0; ii < trimer.ele3.size(); ii++) {
                            dimer_bgcfile << setw(18) << setprecision(9) << trimer.coord3[ii][0]
                                          << setw(18) << setprecision(9) << trimer.coord3[ii][1]
                                          << setw(18) << setprecision(9) << trimer.coord3[ii][2]
                                          << setw(20) << setprecision(5) << trimer.charge3[ii] << "\n";
                        }
                        dimer_bgcfile.close();
                        {
                            stringstream m1_xyzname;
                            stringstream m1_bgcname;
                            m1_xyzname << "monomer_" << trimer.index1 << ".xyz";
                            m1_bgcname << "monomer_" << trimer.index1 << ".bgc";
                            if(access(m1_xyzname.str().c_str(), F_OK)!=0) { 
                                Monomer m1 = trimer.get_m1();
                                
                                ofstream m1_xyzfile;
                                m1_xyzfile.open(m1_xyzname.str());
                                m1_xyzfile << fixed << right;
                                m1_xyzfile << setw(3) << m1.ele.size() << "\n";
                                
                                m1_xyzfile << left << setw(5) << m1.tot_charge
                                           << setw(3) << m1.tot_spin << right << "\n";
                                for (int ii = 0; ii < m1.ele.size(); ii++) {
                                    m1_xyzfile << setw(3) << atmnb2ele.at(m1.ele[ii])
                                               << setw(18) << setprecision(9) << m1.coord[ii][0]
                                               << setw(18) << setprecision(9) << m1.coord[ii][1]
                                               << setw(18) << setprecision(9) << m1.coord[ii][2] << "\n";
                                }
                                m1_xyzfile.close();
                                
                                ofstream m1_bgcfile;
                                m1_bgcfile.open(m1_bgcname.str());
                                m1_bgcfile << fixed << right;
                                for (int ii = 0; ii < m1.bgc_chg.size(); ii++) {
                                    m1_bgcfile << setw(18) << setprecision(9) << m1.bgc_crd[ii][0]
                                               << setw(18) << setprecision(9) << m1.bgc_crd[ii][1]
                                               << setw(18) << setprecision(9) << m1.bgc_crd[ii][2]
                                               << setw(20) << setprecision(5) << m1.bgc_chg[ii] << "\n";
                                }
                                m1_bgcfile.close();
                            }
                        }
                        {
                            stringstream m2_xyzname;
                            stringstream m2_bgcname;
                            m2_xyzname << "monomer_" << trimer.index2 << ".xyz";
                            m2_bgcname << "monomer_" << trimer.index2 << ".bgc";
                            if(access(m2_xyzname.str().c_str(), F_OK)!=0) { 
                                Monomer m2 = trimer.get_m2();
                                
                                ofstream m2_xyzfile;
                                m2_xyzfile.open(m2_xyzname.str());
                                m2_xyzfile << fixed << right;
                                m2_xyzfile << setw(3) << m2.ele.size() << "\n";
                                
                                m2_xyzfile << left << setw(5) << m2.tot_charge
                                           << setw(3) << m2.tot_spin << right << "\n";
                                for (int ii = 0; ii < m2.ele.size(); ii++) {
                                    m2_xyzfile << setw(3) << atmnb2ele.at(m2.ele[ii])
                                               << setw(18) << setprecision(9) << m2.coord[ii][0]
                                               << setw(18) << setprecision(9) << m2.coord[ii][1]
                                               << setw(18) << setprecision(9) << m2.coord[ii][2] << "\n";
                                }
                                m2_xyzfile.close();
                                
                                ofstream m2_bgcfile;
                                m2_bgcfile.open(m2_bgcname.str());
                                m2_bgcfile << fixed << right;
                                for (int ii = 0; ii < m2.bgc_chg.size(); ii++) {
                                    m2_bgcfile << setw(18) << setprecision(9) << m2.bgc_crd[ii][0]
                                               << setw(18) << setprecision(9) << m2.bgc_crd[ii][1]
                                               << setw(18) << setprecision(9) << m2.bgc_crd[ii][2]
                                               << setw(20) << setprecision(5) << m2.bgc_chg[ii] << "\n";
                                }
                                m2_bgcfile.close();
                            }
                        }
                    }
                }

                {
                    stringstream dimer_xyzname;
                    stringstream dimer_bgcname;
                    stringstream dimer_indexes;
                    dimer_indexes << trimer.index1 << "_" << trimer.index3;
                    auto findkey = map_dimer.find(dimer_indexes.str());
                    if (findkey != map_dimer.end()) { 
                        dimer_xyzname << "dimer_" << findkey->second << ".xyz";
                        dimer_bgcname << "dimer_" << findkey->second << ".bgc";
                    }
                    else{
                        map_dimer[dimer_indexes.str()] = count_dimer;
                        dimer_xyzname << "dimer_" << count_dimer << ".xyz";
                        dimer_bgcname << "dimer_" << count_dimer << ".bgc";
                        count_dimer = count_dimer+1;
                    }
                    if(access(dimer_xyzname.str().c_str(), F_OK)!=0) { 
                        int natom_dimer = trimer.get_natom_m1() + trimer.get_natom_m3();
                        int tc_dimer = trimer.tc1 + trimer.tc3;
                        int ts_dimer = 2 * (tc_dimer) + 1;
                        if (trimer.ts3 != 1 or trimer.ts1 != 1){
                            ts_dimer = max(trimer.ts1,trimer.ts3);
                        }
                        ofstream dimer_xyzfile;
                        dimer_xyzfile.open(dimer_xyzname.str());
                        dimer_xyzfile << fixed << right;
                        dimer_xyzfile << setw(3) << natom_dimer << "\n";
                        
                        dimer_xyzfile << left << setw(5) << tc_dimer
                                      << setw(3) << ts_dimer << right << "\n";
                        for (int ii = 0; ii < trimer.ele1.size(); ii++) {
                            dimer_xyzfile << setw(3) << atmnb2ele.at(trimer.ele1[ii])
                                          << setw(18) << setprecision(9) << trimer.coord1[ii][0]
                                          << setw(18) << setprecision(9) << trimer.coord1[ii][1]
                                          << setw(18) << setprecision(9) << trimer.coord1[ii][2] << "\n";
                        }
                        for (int ii = 0; ii < trimer.ele3.size(); ii++) {
                            dimer_xyzfile << setw(3) << atmnb2ele.at(trimer.ele3[ii])
                                          << setw(18) << setprecision(9) << trimer.coord3[ii][0]
                                          << setw(18) << setprecision(9) << trimer.coord3[ii][1]
                                          << setw(18) << setprecision(9) << trimer.coord3[ii][2] << "\n";
                        }
                        dimer_xyzfile.close();
                        
                        ofstream dimer_bgcfile;
                        dimer_bgcfile.open(dimer_bgcname.str());
                        dimer_bgcfile << fixed << right;
                        for (int ii = 0; ii < trimer.bgc_chg.size(); ii++) {
                            dimer_bgcfile << setw(18) << setprecision(9) << trimer.bgc_crd[ii][0]
                                          << setw(18) << setprecision(9) << trimer.bgc_crd[ii][1]
                                          << setw(18) << setprecision(9) << trimer.bgc_crd[ii][2]
                                          << setw(20) << setprecision(5) << trimer.bgc_chg[ii] << "\n";
                        }
                        for (int ii = 0; ii < trimer.ele2.size(); ii++) {
                            dimer_bgcfile << setw(18) << setprecision(9) << trimer.coord2[ii][0]
                                          << setw(18) << setprecision(9) << trimer.coord2[ii][1]
                                          << setw(18) << setprecision(9) << trimer.coord2[ii][2]
                                          << setw(20) << setprecision(5) << trimer.charge2[ii] << "\n";
                        }
                        dimer_bgcfile.close();
                        {
                            stringstream m1_xyzname;
                            stringstream m1_bgcname;
                            m1_xyzname << "monomer_" << trimer.index1 << ".xyz";
                            m1_bgcname << "monomer_" << trimer.index1 << ".bgc";
                            if(access(m1_xyzname.str().c_str(), F_OK)!=0) { 
                                Monomer m1 = trimer.get_m1();
                                
                                ofstream m1_xyzfile;
                                m1_xyzfile.open(m1_xyzname.str());
                                m1_xyzfile << fixed << right;
                                m1_xyzfile << setw(3) << m1.ele.size() << "\n";
                                
                                m1_xyzfile << left << setw(5) << m1.tot_charge
                                           << setw(3) << m1.tot_spin << right << "\n";
                                for (int ii = 0; ii < m1.ele.size(); ii++) {
                                    m1_xyzfile << setw(3) << atmnb2ele.at(m1.ele[ii])
                                               << setw(18) << setprecision(9) << m1.coord[ii][0]
                                               << setw(18) << setprecision(9) << m1.coord[ii][1]
                                               << setw(18) << setprecision(9) << m1.coord[ii][2] << "\n";
                                }
                                m1_xyzfile.close();
                                
                                ofstream m1_bgcfile;
                                m1_bgcfile.open(m1_bgcname.str());
                                m1_bgcfile << fixed << right;
                                for (int ii = 0; ii < m1.bgc_chg.size(); ii++) {
                                    m1_bgcfile << setw(18) << setprecision(9) << m1.bgc_crd[ii][0]
                                               << setw(18) << setprecision(9) << m1.bgc_crd[ii][1]
                                               << setw(18) << setprecision(9) << m1.bgc_crd[ii][2]
                                               << setw(20) << setprecision(5) << m1.bgc_chg[ii] << "\n";
                                }
                                m1_bgcfile.close();
                            }
                        }
                        {
                            stringstream m2_xyzname;
                            stringstream m2_bgcname;
                            m2_xyzname << "monomer_" << trimer.index3 << ".xyz";
                            m2_bgcname << "monomer_" << trimer.index3 << ".bgc";
                            if(access(m2_xyzname.str().c_str(), F_OK)!=0) { 
                                Monomer m2 = trimer.get_m3();
                               
                                ofstream m2_xyzfile;
                                m2_xyzfile.open(m2_xyzname.str());
                                m2_xyzfile << fixed << right;
                                m2_xyzfile << setw(3) << m2.ele.size() << "\n";
                                
                                m2_xyzfile << left << setw(5) << m2.tot_charge
                                           << setw(3) << m2.tot_spin << right << "\n";
                                for (int ii = 0; ii < m2.ele.size(); ii++) {
                                    m2_xyzfile << setw(3) << atmnb2ele.at(m2.ele[ii])
                                               << setw(18) << setprecision(9) << m2.coord[ii][0]
                                               << setw(18) << setprecision(9) << m2.coord[ii][1]
                                               << setw(18) << setprecision(9) << m2.coord[ii][2] << "\n";
                                }
                                m2_xyzfile.close();
                                
                                ofstream m2_bgcfile;
                                m2_bgcfile.open(m2_bgcname.str());
                                m2_bgcfile << fixed << right;
                                for (int ii = 0; ii < m2.bgc_chg.size(); ii++) {
                                    m2_bgcfile << setw(18) << setprecision(9) << m2.bgc_crd[ii][0]
                                               << setw(18) << setprecision(9) << m2.bgc_crd[ii][1]
                                               << setw(18) << setprecision(9) << m2.bgc_crd[ii][2]
                                               << setw(20) << setprecision(5) << m2.bgc_chg[ii] << "\n";
                                }
                                m2_bgcfile.close();
                            }
                        }
                    }
                }

                {// for dimer 3
                    stringstream dimer_xyzname;
                    stringstream dimer_bgcname;
                    stringstream dimer_indexes;
                    dimer_indexes << trimer.index2 << "_" << trimer.index3;
                    auto findkey = map_dimer.find(dimer_indexes.str());
                    if (findkey != map_dimer.end()) { 
                        dimer_xyzname << "dimer_" << findkey->second << ".xyz";
                        dimer_bgcname << "dimer_" << findkey->second << ".bgc";
                    }
                    else{
                        map_dimer[dimer_indexes.str()] = count_dimer;
                        dimer_xyzname << "dimer_" << count_dimer << ".xyz";
                        dimer_bgcname << "dimer_" << count_dimer << ".bgc";
                        count_dimer = count_dimer+1;
                    }
                    if(access(dimer_xyzname.str().c_str(), F_OK)!=0) { 
                        int natom_dimer = trimer.get_natom_m2() + trimer.get_natom_m3();
                        int tc_dimer = trimer.tc2 + trimer.tc3;
                        int ts_dimer = 2 * (tc_dimer) + 1;
                        if (trimer.ts3 != 1 or trimer.ts2 != 1){
                            ts_dimer = max(trimer.ts2,trimer.ts3);
                        }
                        ofstream dimer_xyzfile;
                        dimer_xyzfile.open(dimer_xyzname.str());
                        dimer_xyzfile << fixed << right;
                        dimer_xyzfile << setw(3) << natom_dimer << "\n";
                        
                        dimer_xyzfile << left << setw(5) << tc_dimer
                                      << setw(3) << ts_dimer << right << "\n";
                        for (int ii = 0; ii < trimer.ele2.size(); ii++) {
                            dimer_xyzfile << setw(3) << atmnb2ele.at(trimer.ele2[ii])
                                          << setw(18) << setprecision(9) << trimer.coord2[ii][0]
                                          << setw(18) << setprecision(9) << trimer.coord2[ii][1]
                                          << setw(18) << setprecision(9) << trimer.coord2[ii][2] << "\n";
                        }
                        for (int ii = 0; ii < trimer.ele3.size(); ii++) {
                            dimer_xyzfile << setw(3) << atmnb2ele.at(trimer.ele3[ii])
                                          << setw(18) << setprecision(9) << trimer.coord3[ii][0]
                                          << setw(18) << setprecision(9) << trimer.coord3[ii][1]
                                          << setw(18) << setprecision(9) << trimer.coord3[ii][2] << "\n";
                        }
                        dimer_xyzfile.close();
                        
                        ofstream dimer_bgcfile;
                        dimer_bgcfile.open(dimer_bgcname.str());
                        dimer_bgcfile << fixed << right;
                        for (int ii = 0; ii < trimer.bgc_chg.size(); ii++) {
                            dimer_bgcfile << setw(18) << setprecision(9) << trimer.bgc_crd[ii][0]
                                          << setw(18) << setprecision(9) << trimer.bgc_crd[ii][1]
                                          << setw(18) << setprecision(9) << trimer.bgc_crd[ii][2]
                                          << setw(20) << setprecision(5) << trimer.bgc_chg[ii] << "\n";
                        }
                        for (int ii = 0; ii < trimer.ele1.size(); ii++) {
                            dimer_bgcfile << setw(18) << setprecision(9) << trimer.coord1[ii][0]
                                          << setw(18) << setprecision(9) << trimer.coord1[ii][1]
                                          << setw(18) << setprecision(9) << trimer.coord1[ii][2]
                                          << setw(20) << setprecision(5) << trimer.charge1[ii] << "\n";
                        }
                        dimer_bgcfile.close();
                        {
                            stringstream m1_xyzname;
                            stringstream m1_bgcname;
                            m1_xyzname << "monomer_" << trimer.index1 << ".xyz";
                            m1_bgcname << "monomer_" << trimer.index1 << ".bgc";
                            if(access(m1_xyzname.str().c_str(), F_OK)!=0) {
                                Monomer m1 = trimer.get_m1();
                                
                                ofstream m1_xyzfile;
                                m1_xyzfile.open(m1_xyzname.str());
                                m1_xyzfile << fixed << right;
                                m1_xyzfile << setw(3) << m1.ele.size() << "\n";
                                
                                m1_xyzfile << left << setw(5) << m1.tot_charge
                                           << setw(3) << m1.tot_spin << right << "\n";
                                for (int ii = 0; ii < m1.ele.size(); ii++) {
                                    m1_xyzfile << setw(3) << atmnb2ele.at(m1.ele[ii])
                                               << setw(18) << setprecision(9) << m1.coord[ii][0]
                                               << setw(18) << setprecision(9) << m1.coord[ii][1]
                                               << setw(18) << setprecision(9) << m1.coord[ii][2] << "\n";
                                }
                                m1_xyzfile.close();
                                
                                ofstream m1_bgcfile;
                                m1_bgcfile.open(m1_bgcname.str());
                                m1_bgcfile << fixed << right;
                                for (int ii = 0; ii < m1.bgc_chg.size(); ii++) {
                                    m1_bgcfile << setw(18) << setprecision(9) << m1.bgc_crd[ii][0]
                                               << setw(18) << setprecision(9) << m1.bgc_crd[ii][1]
                                               << setw(18) << setprecision(9) << m1.bgc_crd[ii][2]
                                               << setw(20) << setprecision(5) << m1.bgc_chg[ii] << "\n";
                                }
                                m1_bgcfile.close();
                            }
                        }
                        {
                            stringstream m2_xyzname;
                            stringstream m2_bgcname;
                            m2_xyzname << "monomer_" << trimer.index3 << ".xyz";
                            m2_bgcname << "monomer_" << trimer.index3 << ".bgc";
                            if(access(m2_xyzname.str().c_str(), F_OK)!=0) { 
                                Monomer m2 = trimer.get_m3();
                                
                                ofstream m2_xyzfile;
                                m2_xyzfile.open(m2_xyzname.str());
                                m2_xyzfile << fixed << right;
                                m2_xyzfile << setw(3) << m2.ele.size() << "\n";
                                
                                m2_xyzfile << left << setw(5) << m2.tot_charge
                                           << setw(3) << m2.tot_spin << right << "\n";
                                for (int ii = 0; ii < m2.ele.size(); ii++) {
                                    m2_xyzfile << setw(3) << atmnb2ele.at(m2.ele[ii])
                                               << setw(18) << setprecision(9) << m2.coord[ii][0]
                                               << setw(18) << setprecision(9) << m2.coord[ii][1]
                                               << setw(18) << setprecision(9) << m2.coord[ii][2] << "\n";
                                }
                                m2_xyzfile.close();
                                
                                ofstream m2_bgcfile;
                                m2_bgcfile.open(m2_bgcname.str());
                                m2_bgcfile << fixed << right;
                                for (int ii = 0; ii < m2.bgc_chg.size(); ii++) {
                                    m2_bgcfile << setw(18) << setprecision(9) << m2.bgc_crd[ii][0]
                                               << setw(18) << setprecision(9) << m2.bgc_crd[ii][1]
                                               << setw(18) << setprecision(9) << m2.bgc_crd[ii][2]
                                               << setw(20) << setprecision(5) << m2.bgc_chg[ii] << "\n";
                                }
                                m2_bgcfile.close();
                            }
                        }
                    }
                } 


            }
        }
    }
}

void
write_output(const double & tot_ene,
        vector<vector<double > > & tot_grad,
        const bool & cal_grad){
    stringstream enename;
    enename << "tot_ene.ene";
    cout << "Writing output file ...\n";
    // energy file
    ofstream enefile;
    enefile.open(enename.str());
    enefile << fixed << right;
    enefile << setw(18) << setprecision(9) << tot_ene << "\n";
    enefile.close();
    if (cal_grad) {
        // gradient file
        stringstream gradname;
        gradname << "tot_grad.grad";
        ofstream gradfile;
        gradfile.open(gradname.str());
        gradfile << fixed << right;
        for (int ii = 0; ii < tot_grad.size(); ii++) {
            gradfile << setw(18) << setprecision(9) << tot_grad[ii][0]
                     << setw(18) << setprecision(9) << tot_grad[ii][1]
                     << setw(18) << setprecision(9) << tot_grad[ii][2] << "\n";
        }
        gradfile.close();
    }
    cout << "Done.\n";
}

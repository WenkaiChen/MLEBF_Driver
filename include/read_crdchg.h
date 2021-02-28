#ifndef MBE_DRIVER_READ_CRDCHG_H
#define MBE_DRIVER_READ_CRDCHG_H
#include <fstream>
#include <iostream>
#include <string>
#include <strings.h>
#include <vector>
#include <map>
#include "classes.h"
#include "element.h"

using namespace std;

inline
vector<string> split(const string &s, const string &seperator){
    vector<string> result;
    typedef string::size_type string_size;
    string_size i = 0;

    while(i != s.size()){
        
        int flag = 0;
        while(i != s.size() && flag == 0){
            flag = 1;
            for(string_size x = 0; x < seperator.size(); ++x) {
                if (s[i] == seperator[x]) {
                    ++i;
                    flag = 0;
                    break;
                }
            }
        }

        
        flag = 0;
        string_size j = i;
        while(j != s.size() && flag == 0){
            for(string_size x = 0; x < seperator.size(); ++x)
                if(s[j] == seperator[x]){
                    flag = 1;
                    break;
                }
            if(flag == 0)
                ++j;
        }
        if(i != j){
            result.push_back(s.substr(i, j-i));
            i = j;
        }
    }
    return result;
}

inline
void
read_crdchg(int & n_monomer,
        int & n_atoms,
        vector<SingleMol> & smols,
        const string & fn_cc,
        const double & rcut_boxcenter,
        bool & has_freeze_mol){ 
    ifstream infile;
    infile.open(fn_cc.c_str());
    if (!infile.is_open()){ 
        cout << "Could not open the file " << fn_cc.c_str() << endl;
        cout << "Program terminating.\n";
        exit (EXIT_FAILURE);
    }

    cout << "Reading crdchg file..." << endl << "\tPlease Wait for a moment" << endl;
    string line;
    
    getline(infile, line);
    vector<string> data = split(line,  " ");
    n_monomer = atoi(data[0].c_str());
    n_atoms = atoi(data[1].c_str());
    smols.reserve(n_monomer);
    
    getline(infile, line); 
    int i_mol = 0;
    int tc ;
    int ts ;
    int s  ;
    bool frozen;
    int com; 
    vector<int> ele;
    vector<vector<double > > coord;
    vector<double > charge;
    ele.clear(); coord.clear(); charge.clear();
    for (int i=1;i <= n_atoms;){
        getline(infile, line);
        vector<string> data = split(line, " ");
        if (data.size() == 3){ 
            if (i_mol != 0){ 
                SingleMol smol = SingleMol(s,tc,ts,ele,coord,charge,rcut_boxcenter);
                smols.push_back(smol);
                ele.clear(); coord.clear(); charge.clear();
                if (smol.is_frozen()) {
                    has_freeze_mol = true;
                }
            }
            tc = atoi(data[0].c_str());
            ts = atoi(data[1].c_str());
            s  = atoi(data[2].c_str());
            i_mol += 1;
        }
        else if (data.size() == 4){ 
            if (i_mol != 0){ 
                SingleMol smol = SingleMol(s,tc,ts,ele,coord,charge,frozen);
                smols.push_back(smol);
                ele.clear(); coord.clear(); charge.clear();
            }
            tc = atoi(data[0].c_str());
            ts = atoi(data[1].c_str());
            s  = atoi(data[2].c_str());
            if ((strcasecmp(data[3].c_str(),"A")==0) or (strcasecmp(data[3].c_str(),"active")==0)) {
                frozen = false;
            }
            else if ((strcasecmp(data[3].c_str(),"F")==0) or (strcasecmp(data[3].c_str(),"freeze")==0)) {
                frozen = true;
                has_freeze_mol = true;
            }
            else {
                cout << "Your crdchg file has format error.\nPlease check you crdchg file.\n"
                     << "Program will exit with code 235.\n" << endl;
                exit(235);
            }
            i_mol += 1;
        }
        else if (data.size() == 5){
            int e = ele2atmnb.at(data[0]);
            double x = atof(data[1].c_str());
            double y = atof(data[2].c_str());
            double z = atof(data[3].c_str());
            double c = atof(data[4].c_str());
            vector<double> crd_one_atom;
            crd_one_atom.push_back(x);
            crd_one_atom.push_back(y);
            crd_one_atom.push_back(z);
            coord.push_back(crd_one_atom);
            ele.push_back(e);
            charge.push_back(c);
            i += 1;
        }
    }
    
    if ( rcut_boxcenter == 0.0 ) {
        SingleMol smol = SingleMol(s, tc, ts, ele, coord, charge, frozen);
        smols.push_back(smol);
        ele.clear(); coord.clear(); charge.clear();
    }
    else {
        SingleMol smol = SingleMol(s, tc, ts, ele, coord, charge, rcut_boxcenter);
        smols.push_back(smol);
        ele.clear(); coord.clear(); charge.clear();
        if (smol.is_frozen()) {
            has_freeze_mol = true;
        }
    }

    cout << "Done." << endl;
    infile.close();
};

inline
void
read_crdchg(int & n_monomer,
            int & n_atoms,
            vector<SingleMol> & smols,
            vector<vector<double> >  & lattice_vector,
            const string & fn_cc,
            const double & rcut_boxcenter){ 
    ifstream infile;
    infile.open(fn_cc.c_str());
    if (!infile.is_open()){ 
        cout << "Could not open the file " << fn_cc.c_str() << endl;
        cout << "Program terminating.\n";
        exit (EXIT_FAILURE);
    }

    cout << "Reading crdchg file..." << endl << "\tPlease Wait for a moment" << endl;
    string line;
    
    getline(infile, line);
    vector<string> data = split(line,  " ");
    int size = data.size();
    if (size == 11) {
        cout << "Reading PBC crdchg file...\n";
    } else {
        cout << "You need to write the write lattice vector data at the very first line of the crdchg file.\n"
             << "Check your crdchg file please.\n"
             << "Program will exit with code 266.\n";
        exit(266);
    }
    n_monomer = atoi(data[0].c_str());
    n_atoms = atoi(data[1].c_str());
    lattice_vector.resize(3);
    for (int i=0; i<3; i++){
        lattice_vector[i].resize(3);
        for (int j=0; j<3; j++) {
            lattice_vector[i][j] = atof(data[i*3 +j +2].c_str());
            
        }
    }
    smols.reserve(n_monomer);
    
    getline(infile, line); 
    int i_mol = 0;
    int tc ;
    int ts ;
    int s  ;
    int com; 
    vector<int> ele;
    vector<vector<double > > coord;
    vector<double > charge;
    ele.clear(); coord.clear(); charge.clear();
    for (int i=1;i <= n_atoms;){
        getline(infile, line);
        vector<string> data = split(line, " ");
        if (data.size() == 3){
            if (i_mol != 0){ 
                SingleMol smol = SingleMol(s,tc,ts,ele,coord,charge,rcut_boxcenter);
                smols.push_back(smol);
                ele.clear(); coord.clear(); charge.clear();
            }
            tc = atoi(data[0].c_str());
            ts = atoi(data[1].c_str());
            s  = atoi(data[2].c_str());
            i_mol += 1;
        }
        else if (data.size() == 5){
            int e = ele2atmnb.at(data[0]);
            double x = atof(data[1].c_str());
            double y = atof(data[2].c_str());
            double z = atof(data[3].c_str());
            double c = atof(data[4].c_str());
            vector<double> crd_one_atom;
            crd_one_atom.push_back(x);
            crd_one_atom.push_back(y);
            crd_one_atom.push_back(z);
            coord.push_back(crd_one_atom);
            ele.push_back(e);
            charge.push_back(c);
            i += 1;
        }
    }
    
    SingleMol smol = SingleMol(s,tc,ts,ele,coord,charge,rcut_boxcenter);
    smols.push_back(smol);
    ele.clear(); coord.clear(); charge.clear();

    cout << "Done." << endl;
    infile.close();
};
#endif //MBE_DRIVER_READ_CRDCHG_H

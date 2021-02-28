//
// Created by wenkai on 18-11-29.
//

#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <ctime>
#include "./include/cmdline.h"
#include "./include/read_crdchg.h"
#include "./include/banner.h"
#include "./include/build_minimer.h"
#include "./include/qm_sp_driver.h"
#include "./include/cal_eandg.h"
#include "./include/cal_hessian.h"
#include "./include/dis_mat.h"
#include "./include/pbc_tools.h"
#include "./include/ewald_wk.h"
#include "./include/cal_dipolemoment.h"
#include "./include/cal_dip_mnt_der.h"
#include "./include/cal_transition_dip_mnt.h"

int main(int argc, char *argv[]){
    using namespace std;
    time_t start_time,end_time;
    time(& start_time);
    BANNER();
    // create a parser
    cmdline::parser argparser;
    // add parser
    //argparser.add<string>("input",'i',"input file name",false," ");
    argparser.add<string>("crdchg",'c',"file name of coordinates and charge",true,"coord.crdchg");
    //argparser.add<string>("out",'o',"output file name, default named results.out",false,"results.out");
    argparser.add<double>("rcut_qm",'r',"cut off radius of QM, default 6.0",false,6.0);
    argparser.add<double>("rcut_sr",'s',"cut off radius of short range electrostatic interactions, default 8.0",false,8.0);
    argparser.add<double>("rcut_lr",'l',"cut off radius of long range electrostatic interactions, default 12.0",false,12.0);
    argparser.add<int   >("qm_method",'m',"method for qm/dnn calculation, default 0, Normal Calculation.\n"
                     "                         options:  0 --- Normal calculation.\n"
                     "                                   1 --- Don't cal but just output.\n"
                     "                                9999 --- Debug release. skip sp calculation.",false,0);
    argparser.add<int   >("order",'o',"expansion order, default 2.",false,2);
    argparser.add<string>("box",'b',"filename of box information.",false," ");
    argparser.add<double>("rcut_box_center",'B',"cut off radius of mol from box center. Only useful when --box or -b is used.\nUseless now.\n"
                     "                          default:  15.0", false,15.0);
    argparser.add<int   >("tot_state",'\0',"Only useful for excitation calculations. \n"
                     "                          How many electronic states in total you want to involved.\n"
                     "                          No default value.",false,0);
    argparser.add<int   >("current_state",'\0',"Only useful for excitation calculations. \n"
                     "                          How many electronic states in total you want to involved.\n"
                     "                          No default value.",false,0);
    argparser.add<int   >("directly_delta",'\0',"How many kinds of species will be read delta energy and gradient directly.\n"
                     "                           default: 0 --- All kinds of species will be read with full data.\n"
                     "                                    1 --- One kind of species will be read with delta energy and gradient.\n"
                     "                                    2 --- Two kinds of species will be read with delta energy and gradient.\n"
                     "                                    At the current version, this option must be less than 2,\n"
                     "                                    if you want to add more, please contact the author.",false,0);
    argparser.add<int   >("dd1",'\0',"index 1 for Derectly reading Delta energy and gradient",false,1);
    argparser.add<int   >("dd2",'\0',"index 2 for Derectly reading Delta energy and gradient",false,2);
    argparser.add<int   >("reduce_order",'\0',"If you use this option, some of species could reduce to truncated at the second order.\n"
                     "                           Only Useful When order is equal to 3\n"
                     "                           default: 0 --- All kinds of species will be truncated at the third order.\n"
                     "                                    1 --- One kind of species will be truncated at the third order.\n"
                     "                                    2 --- Two kind of species will be truncated at the third order.\n"
                     "                                    At the current version, this option must be less than 2,\n"
                     "                                    if you want to add more, please contact the author.", false,0);
    argparser.add<int   >("ro1",'\0',"species index 1 for reduce order.",false,0);
    argparser.add<int   >("ro2",'\0',"species index 2 for reduce order.",false,0);
    argparser.add("mono_bgc",'\0',"This must be used for Hybrid EE-2 & NoEE-3 order (--reduce_order must be used).\n"
                     "                           If you use this, then the non-reduced species will be done with an \n"
                     "                                            external background point charge calculation.\n"
                     "                           ------- mono_bgc is the method for simplified MLEBF\n");
    
    argparser.add("partial_nobgc",'\0',"This must be used for Hybrid EE-2 & NoEE-3 order (--reduce_order must be used).\n"
                     "                           If you use this, then the inter-layer interaction will be done with an\n"
                     "                                            external calculation without background charge.\n"
                     "                           This result is not good. We will give up this method.\n");

    argparser.add("judge_frozen",'\0',"if you use this, then you must have rcut_box_center. useless now.");
    argparser.add("no_grad",'\0',"if you use this, no gradients will be calculated.");
    argparser.add("do_hessian",'\0',"if you use this, hessian matreix will be calculated.");
    argparser.add("pbc",'\0', "if you use this, the system will be regarded as periodic system.");
    argparser.add("do_dipole_moment",'\0',"if you use this, dipole moment will be calculated.");
    argparser.add("do_dp_der",'\0', "if you use this, dipole moment derivative will be calculated.");
    argparser.add("do_tans_dp",'\0', "if you use this, transition dipole moment will be calculated.\n"
                     "                                         Please set the value of --excited_fragment_index");
    argparser.add("excited_fragment_index",'\0',"fragment index of the excited fragment. \n"
                     "                           Must be exist in calculating transition dipole moment.",false,0);
    // argparser.add<int   >("electrostatically embedded",'e',)
    argparser.parse_check(argc, argv);

    //string input = argparser.get<string>("input");
    string fn_cc = argparser.get<string>("crdchg");
    //string fn_out= argparser.get<string>("out");
    string fn_box= argparser.get<string>("box");
    double rcut_qm = argparser.get<double>("rcut_qm");
    double rcut_sr = argparser.get<double>("rcut_sr");
    double rcut_lr = argparser.get<double>("rcut_lr");
    int qm_method  = argparser.get<int   >("qm_method");
    int order      = argparser.get<int   >("order");
    double rcut_boxcenter = 0.0;
    if (argparser.exist("judge_frozen")) rcut_boxcenter = argparser.get<double>("rcut_box_center");
    // deal with dd
    bool directlyDelta = argparser.exist("directly_delta") or argparser.exist("dd1") or argparser.exist("dd2");
    vector<int> dd_index;
    if (directlyDelta){
        dd_index.clear();
        for (int i=1;i<argparser.get<int>("directly_delta")+1;i++){
            stringstream ddi_ss;
            ddi_ss << "dd" << i;
            dd_index.push_back(argparser.get<int>(ddi_ss.str()));
        }
    }
    bool reduce_order = argparser.exist("reduce_order") or argparser.exist("ro1") or argparser.exist("ro2");
    bool mono_bgc     = argparser.exist("mono_bgc");
    bool partial_nobgc= argparser.exist("partial_nobgc");
    if (partial_nobgc){
        cout << "This result is not good. We will give up this method.\n"
                "       Program will exit with code 234." << endl;
        exit(234);
    }
    vector<int > ro_index;
    if (reduce_order){
        if (order <= 2){
            cout << "Warning: Order is equal to 2, so your option reduce_order will be useless.\n";
            if (mono_bgc){
                cout << "Error: If you don't used the hybrid 2 & 3 method, the --mono_bgc is not available.\n"
                        "       Program will exit with code 231." << endl;
                exit(231);
            }
            if (partial_nobgc){
                cout << "Error: If you don't used the hybrid 2 & 3 method, the --partial_nobgc is not available.\n"
                        "       Program will exit with code 232." << endl;
                exit(232);
            }
        }
        ro_index.clear();
        for (int i=1;i<argparser.get<int>("reduce_order")+1;i++){
            stringstream roi_ss;
            roi_ss << "ro" << i;
            ro_index.push_back(argparser.get<int>(roi_ss.str()));
        }
        if (mono_bgc){
            cout << "\n\nCaution: You are using --mono_bgc option.\n\n" << endl;
        }
        if (partial_nobgc){
            cout << "\n\nCaution: You are using --partial_bgc option.\n\n" << endl;
        }
    }
    if (not reduce_order){
        if (mono_bgc){
            cout << "Error: If you don't used the hybrid 2 & 3 method, the --mono_bgc is not available.\n"
                    "       Program will exit with code 231." << endl;
            exit(231);
        }
        if (partial_nobgc){
            cout << "Error: If you don't used the hybrid 2 & 3 method, the --partial_nobgc is not available.\n"
                    "       Program will exit with code 232." << endl;
            exit(232);
        }
    }
    if (mono_bgc and partial_nobgc){
        cout << "Error: You are using mono_bgc and partial_nobgc in the same time.\n"
             << "       Program will exit with code 233";
        exit(233);
    }
    bool cal_grad = not argparser.exist("no_grad");
    bool do_hessian = argparser.exist("do_hessian");
    bool do_dp = argparser.exist("do_dipole_moment");
    bool do_dp_der = argparser.exist("do_dp_der");
    if (do_dp_der) {
        do_dp = true;
    }
    bool do_trans_dp = argparser.exist("do_tans_dp");
    int excit_frag_index;
    if (do_trans_dp) {
        excit_frag_index = argparser.get<int   >("excited_fragment_index");
    }
    bool use_pbc  = argparser.exist("pbc");
    // read coordinates and charges and elements from crdchg file.
    vector<SingleMol> smols;
    int n_monomer;
    int n_atoms;
    /*
    system("rm monomer*xyz > /dev/null 2>&1");
    system("rm monomer*bgc > /dev/null 2>&1");
    system("rm dimer*xyz > /dev/null 2>&1");
    system("rm dimer*bgc > /dev/null 2>&1");
    system("rm trimer*bgc > /dev/null 2>&1");
    system("rm trimer*bgc > /dev/null 2>&1");
     */
    if (not use_pbc) {
        bool has_freeze_mol = false;
        read_crdchg(n_monomer, n_atoms, smols, fn_cc, rcut_boxcenter, has_freeze_mol); // read crdchg file and form smols
        if (n_monomer < order) {
            order = n_monomer;
            cout << "Order is larger than the number of fragment, now order has been set to " << order << endl;
        }
        vector<vector<double> > dis_mat;
        compute_dis_mat(dis_mat, smols, fn_box);
        // build monomer
        vector<Monomer> monomers;
        map<string, int> map_monomer;
        build_monomers(monomers, smols, dis_mat, rcut_sr, map_monomer);
        // build dimer
        vector<Dimer> dimers;
        map<string, int> map_dimer;
        build_dimers(dimers, smols, dis_mat, rcut_qm, rcut_sr, map_dimer);
        vector<Trimer> trimers;
        map<string, int> map_trimer;
        if (order == 2) {
            // do qm/dnn calculation with external packages.
            qm_sp_driver(qm_method, monomers, dimers, trimers, map_monomer, map_dimer, map_trimer, rcut_qm, reduce_order, ro_index, partial_nobgc,
                         mono_bgc);
            if (qm_method == 1) {
                cout << "MBE monomer & dimer sp files has been writen. Program will exit with code 0.\n";
            } else if (qm_method == 0) {
                // after sp calculation, calculate the energy and gradients of the system.
                double tot_ene;
                vector<vector<double> > tot_grad;
                if (argparser.exist("tot_state") and argparser.exist("current_state")) {
                    int tot_state = argparser.get<int>("tot_state");
                    int current_state = argparser.get<int>("current_state");
                    cal_eandg(tot_ene, tot_grad, monomers, dimers, trimers, has_freeze_mol, rcut_qm, rcut_sr, rcut_lr,
                              map_monomer, map_dimer, map_trimer, directlyDelta, dd_index,
                              reduce_order, ro_index, partial_nobgc, mono_bgc, cal_grad,
                              2, true, tot_state, current_state);
                    if (do_hessian){
                        vector<vector<double > > hess;
                        cal_hessian(hess,monomers,dimers,trimers,has_freeze_mol,rcut_qm,rcut_sr,rcut_lr,
                                    map_monomer, map_dimer, map_trimer, directlyDelta, dd_index,
                                    reduce_order, ro_index, partial_nobgc, mono_bgc, cal_grad,
                                    2, true, tot_state, current_state);
                    }
                    if (do_dp){
                        vector<double > dip_mnt;
                        cal_dipole_moment(dip_mnt,monomers,dimers,trimers,rcut_qm,rcut_sr,rcut_lr,
                                          map_monomer, map_dimer, map_trimer, directlyDelta, dd_index,
                                          reduce_order, ro_index, partial_nobgc, mono_bgc, cal_grad,
                                          2, true, tot_state, current_state);
                        if (do_dp_der){
                            vector<vector<double > > dipole_moment_der;
                            cal_dipole_moment_der(dipole_moment_der,monomers,dimers,trimers,rcut_qm,rcut_sr,rcut_lr,
                                                  map_monomer, map_dimer, map_trimer, directlyDelta, dd_index,
                                                  reduce_order, ro_index, partial_nobgc, mono_bgc, cal_grad,
                                                  2, true, tot_state, current_state);
                        }
                    }
                    if (do_trans_dp){
//                        cal_trans_dip_mnt();
                    }
                } else if (!(argparser.exist("tot_state")) and !(argparser.exist("current_state"))) {
                    cal_eandg(tot_ene, tot_grad, monomers, dimers, trimers, has_freeze_mol, rcut_qm, rcut_sr, rcut_lr,
                              map_monomer, map_dimer, map_trimer, directlyDelta, dd_index, \
                              reduce_order, ro_index, partial_nobgc, mono_bgc, cal_grad);
                    if (do_hessian){
                        vector<vector<double > > hess;
                        cal_hessian(hess,monomers,dimers,trimers,has_freeze_mol,rcut_qm, rcut_sr, rcut_lr,
                                    map_monomer, map_dimer, map_trimer, directlyDelta, dd_index, \
                                    reduce_order, ro_index, partial_nobgc, mono_bgc, cal_grad);
                    }
                    if (do_dp){
                        vector<double > dip_mnt;
                        cal_dipole_moment(dip_mnt,monomers, dimers, trimers, rcut_qm, rcut_sr, rcut_lr,
                                          map_monomer, map_dimer, map_trimer, directlyDelta, dd_index, \
                                          reduce_order, ro_index, partial_nobgc, mono_bgc, cal_grad);
                        if (do_dp_der){
                            vector<vector<double > > dipole_moment_der;
                            cal_dipole_moment_der(dipole_moment_der,monomers,dimers,trimers,rcut_qm,rcut_sr,rcut_lr,
                                    map_monomer,map_dimer,map_trimer,directlyDelta,dd_index,
                                    reduce_order,ro_index,partial_nobgc,mono_bgc,cal_grad);
                        }
                    }
                    if (do_trans_dp){
//                        cal_trans_dip_mnt();
                    }
                } else {
                    cout
                            << "Command Error: tot_state and current_state must be exist together nor don't exist together.\n"
                            << "Program will exit with code 100";
                    exit(100);
                }
                // output energy and gradient.
                write_output(tot_ene, tot_grad, cal_grad);
                cout << "MBE calculation done. Program will exit with code 0\n";
            } else if (qm_method == 9999) {
                // after sp calculation, calculate the energy and gradients of the system.
                double tot_ene;
                vector<vector<double> > tot_grad;
                cal_eandg(tot_ene, tot_grad, monomers, dimers, trimers, has_freeze_mol, rcut_qm, rcut_sr, rcut_lr,
                          map_monomer, map_dimer, map_trimer, directlyDelta, dd_index,
                          reduce_order, ro_index, partial_nobgc, mono_bgc, cal_grad);
                if (do_hessian){
                    vector<vector<double > > hess;
                    cal_hessian(hess,monomers,dimers,trimers, has_freeze_mol, rcut_qm, rcut_sr, rcut_lr,
                                map_monomer, map_dimer, map_trimer, directlyDelta, dd_index,
                                reduce_order, ro_index, partial_nobgc, mono_bgc, cal_grad);
                }
                if (do_dp){
                    vector<double > dip_mnt;
                    cal_dipole_moment(dip_mnt, monomers, dimers, trimers, rcut_qm, rcut_sr, rcut_lr,
                                      map_monomer, map_dimer, map_trimer, directlyDelta, dd_index,
                                      reduce_order, ro_index, partial_nobgc, mono_bgc, cal_grad);
                    if (do_dp_der){
                        vector<vector<double > > dipole_moment_der;
                        cal_dipole_moment_der(dipole_moment_der,monomers, dimers, trimers, rcut_qm, rcut_sr, rcut_lr,
                                              map_monomer, map_dimer, map_trimer, directlyDelta, dd_index,
                                              reduce_order, ro_index, partial_nobgc, mono_bgc, cal_grad);
                    }
                }
                if (do_trans_dp){
//                    cal_trans_dip_mnt();
                }
                // output energy and gradient.
                write_output(tot_ene, tot_grad, cal_grad);
                cout << "MBE calculation done. Program will exit with code 0\n";
            }
        } else if (order > 2 and order <= 3) {
            // build trimer
            build_trimers(trimers, smols, dis_mat, rcut_qm, rcut_sr, map_trimer, reduce_order, ro_index);
            // do qm/dnn calculation with external packages.
            qm_sp_driver(qm_method, monomers, dimers, trimers, map_monomer, map_dimer, map_trimer, rcut_qm,
                         reduce_order, ro_index, partial_nobgc, mono_bgc, order);
            // after sp calculation, calculate the energy and gradients of the system.
            double tot_ene;
            vector<vector<double> > tot_grad;
            if (argparser.exist("tot_state") and argparser.exist("current_state")) {
                int tot_state = argparser.get<int>("tot_state");
                int current_state = argparser.get<int>("current_state");
                cal_eandg(tot_ene, tot_grad, monomers, dimers, trimers, has_freeze_mol, rcut_qm, rcut_sr, rcut_lr,
                          map_monomer, map_dimer, map_trimer, directlyDelta, dd_index,
                          reduce_order, ro_index, partial_nobgc, mono_bgc, cal_grad,
                          3, true, tot_state, current_state);
                if (do_hessian){
                    vector<vector<double > > hess;
                    cal_hessian(hess, monomers, dimers, trimers, has_freeze_mol, rcut_qm, rcut_sr, rcut_lr,
                                map_monomer, map_dimer, map_trimer, directlyDelta, dd_index,
                                reduce_order, ro_index, partial_nobgc, mono_bgc, cal_grad,
                                3, true, tot_state, current_state);
                }
                if (do_dp){
                    vector<double > dip_mnt;
                    cal_dipole_moment(dip_mnt,monomers,dimers,trimers,rcut_qm,rcut_sr,rcut_lr,
                            map_monomer,map_dimer,map_trimer,directlyDelta,dd_index,
                            reduce_order,ro_index,partial_nobgc,mono_bgc,cal_grad,
                            3,true,tot_state,current_state);
                    if (do_dp_der){
                        vector<vector<double > > dipole_moment_der;
                        cal_dipole_moment_der(dipole_moment_der,monomers,dimers,trimers,rcut_qm,rcut_sr,rcut_lr,
                                              map_monomer,map_dimer,map_trimer,directlyDelta,dd_index,
                                              reduce_order,ro_index,partial_nobgc,mono_bgc,cal_grad,
                                              3,true,tot_state,current_state);
                    }
                }
                if (do_trans_dp){
//                    cal_trans_dip_mnt();
                }
            } else if (!(argparser.exist("tot_state")) and !(argparser.exist("current_state"))) {
                cal_eandg(tot_ene, tot_grad, monomers, dimers, trimers, has_freeze_mol, rcut_qm, rcut_sr, rcut_lr,
                          map_monomer, map_dimer, map_trimer, directlyDelta, dd_index,
                          reduce_order, ro_index, partial_nobgc, mono_bgc, cal_grad,
                          3);
                if (do_hessian){
                    vector<vector<double > > hess;
                    cal_hessian(hess, monomers, dimers, trimers, has_freeze_mol, rcut_qm, rcut_sr, rcut_lr,
                                map_monomer, map_dimer, map_trimer, directlyDelta, dd_index,
                                reduce_order, ro_index, partial_nobgc, mono_bgc, cal_grad,
                                3);
                }
                if (do_dp){
                    vector<double > dip_mnt;
                    cal_dipole_moment(dip_mnt,monomers, dimers, trimers, rcut_qm, rcut_sr, rcut_lr,
                                      map_monomer, map_dimer, map_trimer, directlyDelta, dd_index,
                                      reduce_order, ro_index, partial_nobgc, mono_bgc, cal_grad,
                                      3);
                    if (do_dp_der){
                        vector<vector<double > > dipole_moment_der;
                        cal_dipole_moment_der(dipole_moment_der,monomers, dimers, trimers, rcut_qm, rcut_sr, rcut_lr,
                                              map_monomer, map_dimer, map_trimer, directlyDelta, dd_index,
                                              reduce_order, ro_index, partial_nobgc, mono_bgc, cal_grad,
                                              3);
                    }
                }
                if (do_trans_dp){
//                    cal_trans_dip_mnt();
                }
            } else {
                cout << "Command Error: tot_state and current_state must be exist together nor don't exist together.\n"
                     << "Program will exit with code 100";
                exit(100);
            }
            // output energy and gradient.
            write_output(tot_ene, tot_grad, cal_grad);
            cout << "MBE calculation done. Program will exit with code 0\n";

        }
    }
    else { // use pbc condition.
        bool has_freeze_mol=false;
        vector<double > pbc_box_info;
        vector<vector<double > > lattice_vector;
        read_crdchg(n_monomer, n_atoms, smols, lattice_vector, fn_cc, rcut_boxcenter); // read crdchg file and form smols
        cell2cellpar(lattice_vector, pbc_box_info);
        vector<int> num_uc;
        num_uc.resize(3);
        num_uc[0] = 1 + int(rcut_sr/pbc_box_info[0]);
        num_uc[1] = 1 + int(rcut_sr/pbc_box_info[1]);
        num_uc[2] = 1 + int(rcut_sr/pbc_box_info[2]);
        Pbc_cell_smols unit_cell(pbc_box_info, lattice_vector, smols, false);
        unit_cell.set_external_charges();
        Pbc_cell_smols super_cell = unit_cell.make_supercell(num_uc);
        vector<vector<double> > dis_mat;
        compute_dis_mat(dis_mat, super_cell.smols, fn_box);
        int stt_index_in_scsmols = ((num_uc[0] * (2 * num_uc[1]+1) * (2 * num_uc[2]+1)) + (num_uc[1] * (2 * num_uc[2]+1)) + num_uc[2])*unit_cell.smols.size();
        //int stt_index_in_scsmols = 0; // for debug
        // build monomer
        vector<Monomer > monomers;
        map<string, int> map_monomer;
        build_monomers(monomers,unit_cell,super_cell,stt_index_in_scsmols,num_uc,dis_mat,rcut_sr,map_monomer);
        // build dimer
        vector<Dimer > dimers;
        map<string, int> map_dimer;
        build_dimers(dimers,unit_cell,super_cell,stt_index_in_scsmols,num_uc,dis_mat,rcut_qm,rcut_sr,map_dimer);
        vector<Trimer> trimers;
        map<string, int> map_trimer;
        if (order == 2) {
            // do qm/dnn calculation with external packages.
            qm_sp_driver(qm_method, monomers, dimers, trimers, map_monomer, map_dimer, map_trimer, rcut_qm, reduce_order, ro_index, partial_nobgc,
                         mono_bgc, order, use_pbc);
            if (qm_method == 1) {
                cout << "MBE monomer & dimer sp files has been writen. Program will exit with code 0.\n";
            } else if (qm_method == 0) {
                // after sp calculation, calculate the energy and gradients of the system.
                double tot_ene;
                vector<vector<double> > tot_grad;
                if (argparser.exist("tot_state") and argparser.exist("current_state")) {
                    int tot_state = argparser.get<int>("tot_state");
                    int current_state = argparser.get<int>("current_state");
                    cal_eandg(tot_ene, tot_grad, monomers, dimers, trimers, has_freeze_mol, rcut_qm, rcut_sr, rcut_lr, map_monomer,
                              map_dimer, map_trimer,
                              directlyDelta, dd_index, reduce_order, ro_index, partial_nobgc, mono_bgc, cal_grad, order,
                              true, tot_state, current_state, use_pbc, stt_index_in_scsmols);
                    if (do_hessian){
                        cout << "hessian calculation for PBC condition is not implemented." << endl;
                        exit(1);
                        //vector<vector<double > > hess;
                        //cal_hessian();
                    }
                    if (do_dp){
                        cout << "dipole moment calculation for PBC condition is not implemented." << endl;
                        exit(1);
                    }
                    double ewald_ene;
                    vector<vector<double > > ewald_grad;
                    vector<vector<double > > rec_vec = reciprical_lattice_vectors(lattice_vector);
                    ewald_ene_grad(unit_cell, super_cell, lattice_vector, rec_vec, ewald_ene,
                                   ewald_grad);
                    //cout << "Ewald_gradient = ";
                    //print_matrix(ewald_grad);
                    //cout << setw(18) << setprecision(9) << "ewald_ene = " << ewald_ene << endl;
                    tot_ene += ewald_ene;
                    for (int iatm=0; iatm<tot_grad.size(); iatm++) {
                        tot_grad[iatm][0] = tot_grad[iatm][0] + ewald_grad[iatm][0];
                        tot_grad[iatm][1] = tot_grad[iatm][1] + ewald_grad[iatm][1];
                        tot_grad[iatm][2] = tot_grad[iatm][2] + ewald_grad[iatm][2];
                    }
                } else {
                    cout << "Please give the total_state and the current_state." << endl;
                }
                // output energy and gradient.
                write_output(tot_ene, tot_grad, cal_grad);
                cout << "MBE calculation done. Program will exit with code 0\n";
            }
        } else if (order > 2 and order <= 3) {
            // build trimer
            build_trimers(trimers, unit_cell, super_cell, stt_index_in_scsmols, num_uc, dis_mat, rcut_qm, rcut_sr, map_trimer, reduce_order, ro_index);
            // do qm/dnn calculation with external packages.
            qm_sp_driver(qm_method, monomers, dimers, trimers, map_monomer, map_dimer, map_trimer, rcut_qm, reduce_order, ro_index, partial_nobgc,
                         mono_bgc, order, use_pbc);
            if (qm_method == 1) {
                cout << "MBE monomer & dimer sp files has been writen. Program will exit with code 0.\n";
            } else if (qm_method == 0) {
                // after sp calculation, calculate the energy and gradients of the system.
                double tot_ene;
                vector<vector<double> > tot_grad;
                if (argparser.exist("tot_state") and argparser.exist("current_state")) {
                    int tot_state = argparser.get<int>("tot_state");
                    int current_state = argparser.get<int>("current_state");
                    cal_eandg(tot_ene, tot_grad, monomers, dimers, trimers, has_freeze_mol, rcut_qm, rcut_sr, rcut_lr, map_monomer,
                              map_dimer, map_trimer,
                              directlyDelta, dd_index, reduce_order, ro_index, partial_nobgc, mono_bgc, cal_grad, order,
                              true, tot_state, current_state, use_pbc, stt_index_in_scsmols);
                    if (do_hessian){
                        cout << "hessian calculation for PBC condition is not implemented." << endl;
                        exit(1);
                        //vector<vector<double > > hess;
                        //cal_hessian();
                    }
                    if (do_dp){
                        cout << "dipole moment calculation for PBC condition is not implemented." << endl;
                        exit(1);
                    }
                    double ewald_ene;
                    vector<vector<double > > ewald_grad;
                    vector<vector<double > > rec_vec = reciprical_lattice_vectors(lattice_vector);
                    ewald_ene_grad(unit_cell, super_cell, lattice_vector, rec_vec, ewald_ene,
                                   ewald_grad);
                    //cout << "Ewald_gradient = ";
                    //print_matrix(ewald_grad);
                    //cout << setw(18) << setprecision(9) << "ewald_ene = " << ewald_ene << endl;
                    tot_ene += ewald_ene;
                    for (int iatm=0; iatm<tot_grad.size(); iatm++) {
                        tot_grad[iatm][0] = tot_grad[iatm][0] + ewald_grad[iatm][0];
                        tot_grad[iatm][1] = tot_grad[iatm][1] + ewald_grad[iatm][1];
                        tot_grad[iatm][2] = tot_grad[iatm][2] + ewald_grad[iatm][2];
                    }
                } else {
                    cout << "Please give the total_state and the current_state." << endl;
                }
                // output energy and gradient.
                write_output(tot_ene, tot_grad, cal_grad);
                cout << "MBE calculation done. Program will exit with code 0\n";
            }
        }
    }
    time(& end_time);
    cout<<"Running Time for MBE_Driver is "<< (double)(end_time-start_time) << " Second." << endl;
    return 0;
}
 
#include "../include/ewald_wk.h"

void get_crd_chg_mol(const Pbc_cell_smols & uc,
        vector<vector<double > > & coord,
        vector<double > & charge,
        vector<vector<int > > & sc_index,
        vector<int > & mol_index,
        vector<double > & masses){
    int nmol_uc = uc.smols.size();
    int natom = 0;

    for (int imol_uc=0; imol_uc<nmol_uc; imol_uc++) { 
        SingleMol smol_uc = uc.smols[imol_uc];
        for (int inatom_uc = 0; inatom_uc < smol_uc.ele.size(); inatom_uc++) {
            coord[natom] = 0.1 * smol_uc.coord[inatom_uc]; 
            sc_index[natom] = smol_uc.sc_index;
            mol_index[natom] = smol_uc.mol_index;
            charge[natom] = smol_uc.charge[inatom_uc];
            masses[natom] = ele2mass.at(smol_uc.ele[inatom_uc]);
            natom += 1;
        }
    }
}

void get_crd_chg_mol_uc(const Pbc_cell_smols & uc,
                       vector<vector<double > > & coord,
                       vector<double > & charge,
                       vector<int > & mol_index,
                       vector<double > & masses){
    int nmol_uc = uc.smols.size();
    int natom = 0;

    for (int imol_uc=0; imol_uc<nmol_uc; imol_uc++) { 
        SingleMol smol_uc = uc.smols[imol_uc];
        for (int inatom_uc = 0; inatom_uc < smol_uc.ele.size(); inatom_uc++) {
            coord[natom] = 0.1 * smol_uc.coord[inatom_uc]; 
            mol_index[natom] = smol_uc.mol_index;
            charge[natom] = smol_uc.charge[inatom_uc];
            masses[natom] = ele2mass.at(smol_uc.ele[inatom_uc]);
            natom += 1;
        }
    }
    const vector<vector<double > > & lv = uc.lattice_vector;

    put_in_cell(coord, lv);
}

void put_in_cell(vector<vector<double > > & coord,
        const vector<vector<double > > & lv ){

    vector<vector<double > > frac_crd;

    car2frac(coord,lv*0.1,frac_crd);
    int natom = coord.size();
    vector<vector<double > > frac_in_cell(natom, vector<double > (3,0.0));
    for (int i=0; i<natom; i++){
        for (int j=0; j<3; j++){
            frac_in_cell[i][j] = frac_crd[i][j]-(int)(frac_crd[i][j]/1.0)*1.0;
            if (frac_in_cell[i][j] < 0){
                frac_in_cell[i][j] = 1.0 + frac_in_cell[i][j];
            }
        }
    }

    frac2car(frac_in_cell,lv*0.1,coord);

}

void get_PBC(double &RX, double &RY, double &RZ,
             const double & half_L, const double & half_B, const double & half_H,
             const double & L, const double & B, const double & H)
{
    if (abs(RX) > half_L)
    {
        if(RX <0 )	RX = RX + L;
        else		RX = RX - L;
    }
    if (abs(RY) > half_B)
    {
        if(RY <0 )	RY = RY + B;
        else		RY = RY - B;
    }
    if (abs(RZ) > half_H)
    {
        if(RZ <0 )	RZ = RZ + H;
        else		RZ = RZ - H;
    }
}

void ewald_ene_grad(const Pbc_cell_smols & unit_cell,
                    const Pbc_cell_smols & super_cell,
                    const vector<vector<double > > & lattice_vec,
                    const vector<vector<double > > & rec_lac_vec,
                    double & ewald_ene,
                    vector<vector<double > > & ewald_grad) {
    double kjmol2hartree = 3.808798e-4;
    int nmol_uc = unit_cell.smols.size();
    int natom = 0;
    for (int imol_uc=0; imol_uc<nmol_uc; imol_uc++) { 
        SingleMol smol_uc = unit_cell.smols[imol_uc];
        for (int inatom_uc = 0; inatom_uc < smol_uc.ele.size(); inatom_uc++) {
            natom += 1;
        }
    }
    vector<vector<double > > coord(natom, vector<double >(3,0)); 
    vector<double > charge(natom,0); 
    vector<int > mol_index(natom, 0);
    vector<double > masses(natom, 0); 
    get_crd_chg_mol_uc(unit_cell,coord,charge,mol_index,masses);


    vector<double > pbc_info;
    cell2cellpar(lattice_vec,pbc_info);


    vector<vector<double > > R_force(natom, vector<double >(3,0));
    vector<vector<double > > K_force(natom, vector<double >(3,0));
    vector<vector<double > > Intra_force(natom, vector<double >(3,0));
    vector<vector<double > > DIP_force(natom, vector<double >(3,0));
    vector<vector<double > > FORCE(natom, vector<double >(3,0));
    vector<vector<double > > COM(natom, vector<double >(3,0));
    double L = pbc_info[0]*0.1; double B = pbc_info[1]*0.1; double H = pbc_info[2]*0.1;
    double box_VOL = L*B*H;
    int N=natom;
    double ONE_PI_EPS0 =138.9354859; 
    double U_CR, U_CF, U_CS, U_CLR, U_INTRA;
    double Vir_XX = 0.; double Vir_YY = 0.;double Vir_ZZ = 0.;
    double PXX, PYY, PZZ;
    
    double ALPHA;
    double half_L = L/2.;double half_B = B/2.;double half_H = H/2.;
    double max_box_length;
    max_box_length = L;
    if (max_box_length<B) max_box_length = B;
    if (max_box_length<H) max_box_length = H;
    ALPHA = 0.37; // not auto-determined
    int nmax = 20; // not auto-determined
    int nmax_x = nmax; int nmax_y = nmax; int nmax_z = nmax;
    int kmax = 13; // not auto-determined
    int kmax_x = kmax; int kmax_y = kmax; int kmax_z = kmax;
    // This part is for auto-determined the paramters of ewald summation.
    double accuracy_fac = 15; // This value is taken from HMBI packge of Gregory Beran group.
    ALPHA = sqrt(M_PI) * (1/L + 1/B + 1/H)/3.0;
    double rc = accuracy_fac/ALPHA;
    nmax_x = (int) ceil(rc/L);
    nmax_y = (int) ceil(rc/B);
    nmax_z = (int) ceil(rc/H);
    kmax_x = (int) ceil(accuracy_fac*L*ALPHA/M_PI);
    kmax_y = (int) ceil(accuracy_fac*B*ALPHA/M_PI);
    kmax_z = (int) ceil(accuracy_fac*H*ALPHA/M_PI);
    // ----------------- FIND COM ----------------------
    for(int i = 0; i < N; i ++) {
        COM[i] = unit_cell.smols[mol_index[i]].comass*0.1;
    }
    // -------------- FIND COM END ---------------------
    // ---------------- REAL SPACE ---------------------
    {
        double qi, qj, rsqr, r;
        double rx,ry,rz, tmp;
        double RX, RY, RZ;
        double GAMMA = 2*ALPHA/sqrt(M_PI);
        double ALP2 = ALPHA*ALPHA;

        U_CR = 0;
        //output<<"calculating real space values..."<<endl;
        for(int i = 0; i < N-1; i++)
        {
            qi = charge[i];
            for(int j = i+1; j < N; j++)
            {
                qj = charge[j];
                rx = coord[i][0] - coord[j][0];
                ry = coord[i][1] - coord[j][1];
                rz = coord[i][2] - coord[j][2];
                get_PBC(rx, ry, rz,half_L,half_B,half_H,L,B,H);
                for(int nx = -nmax_x; nx <= nmax_x ; nx++)
                {
                    for(int ny = -nmax_y; ny <= nmax_y; ny++)
                    {
                        for(int nz = -nmax_z; nz <= nmax_z; nz++)
                        {
                            if(nx == 0 && ny == 0 && nz == 0) //central cell
                            {
                                //ignore i == j
                                if(mol_index[i] != mol_index[j])
                                {
                                    rsqr = rx*rx + ry*ry + rz*rz;
                                    if(rsqr != 0)
                                    {
                                        r = sqrt(rsqr);
                                        U_CR += ONE_PI_EPS0 * qi * qj * erfc(ALPHA * r) / r;
                                        tmp = (GAMMA * exp(-ALP2*rsqr) + erfc(ALPHA*r)/r)/rsqr;
                                        tmp = tmp * ONE_PI_EPS0 * qi * qj;
                                        R_force[i][0] += tmp * rx;
                                        R_force[i][1] += tmp * ry;
                                        R_force[i][2] += tmp * rz;

                                        R_force[j][0] -= tmp * rx;
                                        R_force[j][1] -= tmp * ry;
                                        R_force[j][2] -= tmp * rz;

                                        /********* virial *********/
                                        Vir_XX += tmp * rx* rx;
                                        Vir_YY += tmp * ry* ry;
                                        Vir_ZZ += tmp * rz* rz;
                                    }
                                }
                                else
                                {
                                    //ignore, because same molecule
                                }
                            }
                            else //image cells
                            {
                                RX = rx + nx*L;
                                RY = ry + ny*B;
                                RZ = rz + nz*H;
                                rsqr = RX*RX + RY*RY + RZ*RZ;
                                if(rsqr != 0)
                                {
                                    r = sqrt(rsqr);
                                    U_CR += ONE_PI_EPS0 * qi * qj * erfc(ALPHA * r) / r;
                                    tmp = (GAMMA * exp(-ALP2*rsqr) + erfc(ALPHA*r)/r)/rsqr;
                                    tmp = tmp * ONE_PI_EPS0 * qi * qj;
                                    R_force[i][0] += tmp * RX;
                                    R_force[i][1] += tmp * RY;
                                    R_force[i][2] += tmp * RZ;

                                    R_force[j][0] -= tmp * RX;
                                    R_force[j][1] -= tmp * RY;
                                    R_force[j][2] -= tmp * RZ;

                                    /********* virial *********/
                                    Vir_XX += tmp * RX* RX;
                                    Vir_YY += tmp * RY* RY;
                                    Vir_ZZ += tmp * RZ* RZ;
                                }
                            }
                        }
                    }
                }
            }
        }
        //output<<"done with real space"<<endl;
    }
    // ---------------REAL SPACE END -------------------
    // --------------- FOURIER SPACE -------------------
    {
        complex<double> ak;
        double ksqr, kdotr, qi, tmp, tmp2;
        double a,b,akak, mx,my,mz, chak;
        double GAMMA = -1/(4*ALPHA*ALPHA);
        double recip = ONE_PI_EPS0 * 2* M_PI / box_VOL;
        //output.open("DEBUG_details.xls");
        //output<<"i\tkz\tmz\tqi\tkdotr\ttmp2\tforcez"<<endl;
        U_CF = 0;

        for(int kx = -kmax_x; kx <= kmax_x; kx++)
        {
            mx = 2*M_PI*kx/L;
            for(int ky = -kmax_y; ky <= kmax_y; ky++)
            {
                my = 2*M_PI*ky/B;
                for(int kz = -kmax_z; kz <= kmax_z; kz++)
                {
                    mz = 2*M_PI*kz/H;
                    ksqr = mx*mx + my*my + mz*mz;
                    if(ksqr != 0)
                    {
                        ak.real(0.);
                        ak.imag(0.);
                        //ak = (0.,0.);
                        for(int i = 0; i < N; i++)
                        {
                            kdotr = mx*coord[i][0] + my*coord[i][1] + mz*coord[i][2];
                            qi = charge[i];
                            ak.real(ak.real()+qi*cos(kdotr));
                            ak.imag(ak.imag()-qi*sin(kdotr));
                            //complex<double > tmp_ak = (qi*cos(kdotr),-1.*qi*sin(kdotr));
                            //ak += tmp_ak;
                        }
                        a = ak.real();
                        b = ak.imag();
                        akak = (a*a + b*b);
                        tmp = recip * exp(GAMMA * ksqr)/ksqr;
                        U_CF += tmp * akak;
                        chak = (2/ksqr) - 2*GAMMA;
                        /********* virial *********/
                        Vir_XX += (1 - chak * mx*mx) *tmp * akak;
                        Vir_YY += (1 - chak * my*my) *tmp * akak;
                        Vir_ZZ += (1 - chak * mz*mz) *tmp * akak;
                        for(int i = 0; i < N; i++)
                        {
                            kdotr = mx*coord[i][0] + my*coord[i][1] + mz*coord[i][2];
                            qi = charge[i];
                            tmp2 = 2*tmp*qi*(sin(kdotr) * a + cos(kdotr) * b);
                            K_force[i][0] += tmp2 * mx;
                            K_force[i][1] += tmp2 * my;
                            K_force[i][2] += tmp2 * mz;
                        }
                    }
                }
            }
        }
    }
    // --------------FOURIER SPACE END -----------------
    // ----------------- SELF --------------------------
    {
        double SR_PI = -ALPHA/sqrt(M_PI);

        U_CS = 0;
        for(int i = 0; i < N; i++)
        {
            U_CS += charge[i] * charge[i];
        }
        U_CS *= SR_PI * ONE_PI_EPS0;
    }
    // ----------------SELF END ------------------------
    // ---------------- DIPOLE -------------------------
    {
        double M_sqr, Mx, My, Mz;
        double qi ;
        U_CLR = 0;
        M_sqr = Mx = My = Mz = 0;
        for(int i = 0; i < N; i++)
        {
            qi = charge[i];
            Mx += qi * coord[i][0];
            My += qi * coord[i][1];
            Mz += qi * coord[i][2];
        }
        M_sqr = Mx*Mx + My*My + Mz*Mz;
        U_CLR = ONE_PI_EPS0 * 2 * M_PI * M_sqr / (box_VOL*3);
        for(int i = 0; i < N; i++)
        {
            qi = charge[i];
            DIP_force[i][0] -= ONE_PI_EPS0* qi * 4*M_PI * Mx / (3*box_VOL);
            DIP_force[i][1] -= ONE_PI_EPS0* qi * 4*M_PI * My / (3*box_VOL);
            DIP_force[i][2] -= ONE_PI_EPS0* qi * 4*M_PI * Mz / (3*box_VOL);
            /********* virial *********/
            Vir_XX += ONE_PI_EPS0 * 2*M_PI*(M_sqr - 2*Mx*Mx)/(box_VOL*box_VOL*3);
            Vir_YY += ONE_PI_EPS0 * 2*M_PI*(M_sqr - 2*My*My)/(box_VOL*box_VOL*3);
            Vir_ZZ += ONE_PI_EPS0 * 2*M_PI*(M_sqr - 2*Mz*Mz)/(box_VOL*box_VOL*3);
        }
    }
    // ---------------DIPOLE END -----------------------
    // ---------------- INTRA_MOL ----------------------
    {
        double rx, ry, rz, qi, qj;
        double rsqr, r, tmp;
        double KAPPA = 2*ALPHA/sqrt(M_PI);
        double ALPSQR = ALPHA * ALPHA;
        double FST;

        U_INTRA = 0;

        for(int i = 0; i < N-1; i++) // N is the number of atoms
        {
            for(int j = i+1; j < N; j++)
            {
                if(mol_index[i] == mol_index[j]) // same molecule
                {
                    qi = charge[i];
                    qj = charge[j];
                    rx = coord[i][0] - coord[j][0];
                    ry = coord[i][1] - coord[j][1];
                    rz = coord[i][2] - coord[j][2];
                    rsqr = rx*rx + ry*ry + rz*rz;
                    r = sqrt(rsqr);
                    FST =erf(ALPHA*r);
                    U_INTRA -= ONE_PI_EPS0*qi*qj*FST/r;
                    tmp = ONE_PI_EPS0*qi*qj*(KAPPA*exp(-ALPSQR*rsqr)/rsqr - FST/(r*rsqr));

                    Intra_force[i][0] += tmp * rx;
                    Intra_force[i][1] += tmp * ry;
                    Intra_force[i][2] += tmp * rz;

                    Intra_force[j][0] -= tmp * rx;
                    Intra_force[j][1] -= tmp * ry;
                    Intra_force[j][2] -= tmp * rz;
                    /********* virial *********/
                    Vir_XX += tmp * rx * rx;
                    Vir_YY += tmp * ry * ry;
                    Vir_ZZ += tmp * rz * rz;
                }
            }
        }
    }
    // ---------------- INTRA_MOL END ----------------------
    ewald_ene = (U_CF + U_CS + U_CR + U_INTRA) * kjmol2hartree;
    // ------------------ DISP_FORCE -----------------------
    {
        double nm2bohr = 18.89726;
        ewald_grad.resize(N);
        for(int i = 0; i < N; i++)
        {
            FORCE[i][0] = R_force[i][0] + K_force[i][0] + Intra_force[i][0];
            FORCE[i][1] = R_force[i][1] + K_force[i][1] + Intra_force[i][1];
            FORCE[i][2] = R_force[i][2] + K_force[i][2] + Intra_force[i][2];
            ewald_grad[i].resize(3);
            ewald_grad[i][0] = FORCE[i][0] * kjmol2hartree / nm2bohr;
            ewald_grad[i][1] = FORCE[i][1] * kjmol2hartree / nm2bohr;
            ewald_grad[i][2] = FORCE[i][2] * kjmol2hartree / nm2bohr;
        }
    }
    // ----------------- DISP_FORCE END --------------------

    // test: add ee ene and gradient

    double sc_ee_ene;
    vector<vector<double > > sc_ee_grad(natom, vector<double >(3,0));
    uc_ee_in_sc_ene_grad(unit_cell, super_cell, lattice_vec,sc_ee_ene,sc_ee_grad);
    //cout << "sc_ee_ene = " << sc_ee_ene;
    //print_matrix(sc_ee_grad);
    ewald_ene = ewald_ene - sc_ee_ene;
    ewald_grad = ewald_grad - sc_ee_grad;

    write_ewald_ene_grad(ewald_ene,ewald_grad);

}

void uc_ee_in_sc_ene_grad(const Pbc_cell_smols & unit_cell,
                          const Pbc_cell_smols & super_cell,
                          const vector<vector<double > > & lattice_vec,
                          double & sc_ee_ene,
                          vector<vector<double > > & sc_ee_grad){
    int nmol_uc = unit_cell.smols.size();
    int natom_uc = 0;
    for (int imol_uc=0; imol_uc<nmol_uc; imol_uc++) { 
        SingleMol smol_uc = unit_cell.smols[imol_uc];
        for (int inatom_uc = 0; inatom_uc < smol_uc.ele.size(); inatom_uc++) {
            
            natom_uc += 1;
        }
    }
    int nmol_sc = super_cell.smols.size();
    int natom_sc = 0;
    for (int imol_sc=0; imol_sc<nmol_sc; imol_sc++) { 
        SingleMol smol_sc = super_cell.smols[imol_sc];
        for (int inatom_sc = 0; inatom_sc < smol_sc.ele.size(); inatom_sc++) {
            natom_sc += 1;
        }
    }
    
    vector<vector<double > > coord_uc(natom_uc, vector<double >(3,0)); // unit: nm
    vector<double > charge_uc(natom_uc,0); // unit: a.u.
    vector<int > mol_index_uc(natom_uc,0);
    vector<vector<int > > sc_index_uc(natom_uc, vector<int >(3, 0));
    vector<double > masses_uc(natom_uc, 0); // unit: amu
    get_crd_chg_mol(unit_cell,coord_uc,charge_uc,sc_index_uc,mol_index_uc,masses_uc);
    
    vector<vector<double > > coord_sc(natom_sc, vector<double >(3,0)); // unit: nm
    vector<double > charge_sc(natom_sc,0); // unit: a.u.
    vector<int > mol_index_sc(natom_sc,0);
    vector<vector<int > > sc_index_sc(natom_sc, vector<int >(3, 0));
    vector<double > masses_sc(natom_sc, 0); // unit: amu
    get_crd_chg_mol(super_cell,coord_sc,charge_sc,sc_index_sc,mol_index_sc,masses_sc);

    vector<double > pbc_info;
    cell2cellpar(lattice_vec,pbc_info);
    
    double angstrom2bohr = 1.889726;
    pbc_info[0] = pbc_info[0] * angstrom2bohr;
    pbc_info[1] = pbc_info[1] * angstrom2bohr;
    pbc_info[2] = pbc_info[2] * angstrom2bohr;

    
    double nm2bohr = 18.89726;
    coord_uc = coord_uc * nm2bohr;
    coord_sc = coord_sc * nm2bohr;

    
    sc_ee_ene = 0.0;

    
    for(int i=0; i<natom_uc; i++){
        vector<double > crd_i = coord_uc[i];
        double qi = charge_uc[i];
        for (int j=0; j<natom_sc; j++){
            if (not vector_equal_or_not(sc_index_uc[i],sc_index_sc[j])) { 
                vector<double > crd_j = coord_sc[j];
                vector<double > r_ij_vec = crd_j - crd_i;
                double r_ij = sqrt(norm2(r_ij_vec));
                double qj = charge_sc[j];
                //energy
                sc_ee_ene += 0.5*qi*qj/r_ij;
                //gradient
                double r3 = r_ij * r_ij * r_ij;
                vector<double > grad = r_ij_vec * (-1. * qi * qj/r3);
                for (int xyz=0; xyz<3; xyz++){
                    sc_ee_grad[i][xyz] += grad[xyz];
                }
            }
        }
    }
    
    for(int i=0; i<natom_uc-1; i++){
        vector<double > crd_i = coord_uc[i];
        double qi = charge_uc[i];
        for (int j=i+1; j<natom_uc; j++){
            if(mol_index_uc[i] != mol_index_uc[j]) {//ignore i == j
                vector<double> crd_j = coord_uc[j];
                vector<double> r_ij_vec = crd_j - crd_i;
                double r_ij = sqrt(norm2(r_ij_vec));
                double qj = charge_uc[j];
                //energy
                sc_ee_ene += qi * qj / r_ij;
                //gradient
                double r3 = r_ij * r_ij * r_ij;
                vector<double> grad = r_ij_vec * (-1. * qi * qj / r3);
                for (int xyz = 0; xyz < 3; xyz++) {
                    sc_ee_grad[i][xyz] += grad[xyz];
                    sc_ee_grad[j][xyz] -= grad[xyz];
                }
            }
        }
    }
}

void write_ewald_ene_grad(
        const double & ewald_ene,
        const vector<vector<double > > & ewald_grad){
    stringstream enename;
    enename << "ewald_ene.ene";
    // energy file
    ofstream enefile;
    enefile.open(enename.str());
    enefile << fixed << right;
    enefile << setw(18) << setprecision(9) << ewald_ene << "\n";
    enefile.close();
    // gradient file
    stringstream gradname;
    gradname << "ewald_grad.grad";
    ofstream gradfile;
    gradfile.open(gradname.str());
    gradfile << fixed << right;
    for (int ii = 0; ii < ewald_grad.size(); ii++) {
        gradfile << setw(18) << setprecision(9) << ewald_grad[ii][0]
                 << setw(18) << setprecision(9) << ewald_grad[ii][1]
                 << setw(18) << setprecision(9) << ewald_grad[ii][2] << "\n";
    }
    gradfile.close();
    cout << "Done.\n";
}

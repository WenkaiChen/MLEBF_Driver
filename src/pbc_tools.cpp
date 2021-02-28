#include "../include/pbc_tools.h"

static const double PI = 3.1415926535898;

void
cell2cellpar(const vector<vector<double > > & lattice_vector,
             vector<double > & pbc_box_info){
    pbc_box_info.clear();
    pbc_box_info.resize(6);
    for (int i=0; i<3; i++){
        pbc_box_info[i] = sqrt(norm2(lattice_vector[i]));
    }
    for (int i=3; i<6; i++){
        int j = (i - 1)%3;
        int k = (i - 2)%3;
        double ll = pbc_box_info[j] * pbc_box_info[k];
        if (ll > 1e-16){
            double x = vector_multiply(lattice_vector[j], lattice_vector[k]) / ll;
            pbc_box_info[i] = 180.0 / PI * acos(x);
        }
    }
}

void
car2frac(const vector<vector<double > > & coord_cartesian,
         const vector<vector<double > > & lattice_vector,
         vector<vector<double > > & coord_fraction){
    int natoms = coord_cartesian.size();
    coord_fraction.clear();
    coord_fraction.resize(natoms);
    vector<vector<double > > M = mat_transpose(lattice_vector);
    vector<vector<double > > U;
    mat_inverse(M,U,3);
    for (int i=0; i<natoms; i++){
        coord_fraction[i] = mat_mul_vec(U,coord_cartesian[i]);
    }
}

void
frac2car(const vector<vector<double > > & coord_fraction,
         const vector<vector<double > > & lattice_vector,
         vector<vector<double > > & coord_cartesian){
    int natoms = coord_fraction.size();
    coord_cartesian.clear();
    coord_cartesian.resize(natoms);
    vector<vector<double > > M = mat_transpose(lattice_vector);
    for (int i=0; i<natoms; i++){
        coord_cartesian[i] = mat_mul_vec(M,coord_fraction[i]);
    }
}

vector<vector<double > >
reciprical_lattice_vectors(const vector<vector<double > > & lat_vec){
    vector<vector<double > > rec_vec;
    double temp;
    rec_vec.resize(3);
    rec_vec[0].resize(3);
    rec_vec[1].resize(3);
    rec_vec[2].resize(3);

    rec_vec[0][0] = lat_vec[1][1] * lat_vec[2][2] - lat_vec[1][2] * lat_vec[2][1];
    rec_vec[0][1] = lat_vec[1][2] * lat_vec[2][0] - lat_vec[1][0] * lat_vec[2][2];
    rec_vec[0][2] = lat_vec[1][0] * lat_vec[2][1] - lat_vec[1][1] * lat_vec[2][0];
    temp = lat_vec[0][0]*rec_vec[0][0] + lat_vec[0][1]*rec_vec[0][1] + lat_vec[0][2]*rec_vec[0][2];
    rec_vec[0][0] = rec_vec[0][0]/temp;
    rec_vec[0][1] = rec_vec[0][1]/temp;
    rec_vec[0][2] = rec_vec[0][2]/temp;

    rec_vec[1][0] = lat_vec[0][1] * lat_vec[2][2] - lat_vec[0][2] * lat_vec[2][1];
    rec_vec[1][1] = lat_vec[0][2] * lat_vec[2][0] - lat_vec[0][0] * lat_vec[2][2];
    rec_vec[1][2] = lat_vec[0][0] * lat_vec[2][1] - lat_vec[0][1] * lat_vec[2][0];
    temp = lat_vec[1][0]*rec_vec[1][0] + lat_vec[1][1]*rec_vec[1][1] + lat_vec[1][2]*rec_vec[1][2];
    rec_vec[1][0] = rec_vec[1][0]/temp;
    rec_vec[1][1] = rec_vec[1][1]/temp;
    rec_vec[1][2] = rec_vec[1][2]/temp;

    rec_vec[2][0] = lat_vec[0][1] * lat_vec[1][2] - lat_vec[0][2] * lat_vec[1][1];
    rec_vec[2][1] = lat_vec[0][2] * lat_vec[1][0] - lat_vec[0][0] * lat_vec[1][2];
    rec_vec[2][2] = lat_vec[0][0] * lat_vec[1][1] - lat_vec[0][1] * lat_vec[1][0];
    temp = lat_vec[2][0]*rec_vec[2][0] + lat_vec[2][1]*rec_vec[2][1] + lat_vec[2][2]*rec_vec[2][2];
    rec_vec[2][0] = rec_vec[2][0]/temp;
    rec_vec[2][1] = rec_vec[2][1]/temp;
    rec_vec[2][2] = rec_vec[2][2]/temp;

    return rec_vec;
}

double
compute_volume_cell(const vector<vector<double > > & lattice_vector){
    return mat_determinant(lattice_vector,3);
}

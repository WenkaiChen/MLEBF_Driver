// get if from https://github.com/kain88-de/pbc_distances/blob/master/pbc_distances/vector_helper.h
//

#ifndef VECTOR_HELPER_H
#define VECTOR_HELPER_H

#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

using namespace std;

template <typename T>
bool vector_equal_or_not(const std::vector<T> &a, const std::vector<T> &b){
    int size_a = a.size();
    int size_b = b.size();
    if (size_a != size_b) {
        return false;
    }
    for (int i=0; i<size_a; i++){
        if (a[i] != b[i]){
            return false;
        }
    }
    return true;
}

template <typename T>
bool matrix_equal_or_not(const std::vector<std::vector<T > > & a,
        const std::vector<std::vector<T > > & b){
    int size_a = a.size();
    int size_b = b.size();
    if (size_a != size_b) {
        return false;
    }
    for (int i=0; i<size_a; i++) {
        int size_a_i = a[i].size();
        int size_b_i = b[i].size();
        if (size_a_i != size_b_i) {
            return false;
        }
        for (int j=0; j<size_a_i; j++) {
            if (a[i][j] != b[i][j]){
                return false;
            }
        }
    }
    return true;
}

template <typename T>
std::vector<T> operator*(const std::vector<T> &a, const T b) {
    auto res = std::vector<T>(a.size());
    for (auto i = 0u; i < a.size(); ++i) {
        res[i] = a[i] * b;
    }
    return res;
}

template <typename T>
std::vector<T> operator*(const T b, const std::vector<T> &a) {
    auto res = std::vector<T>(a.size());
    for (auto i = 0u; i < a.size(); ++i) {
        res[i] = a[i] * b;
    }
    return res;
}

template <typename T>
std::vector<vector<T > > operator*(const std::vector<vector<T> > &a, const T b) {
    auto res = std::vector<vector<T > >(a.size());
    for (auto i = 0u; i < a.size(); ++i) {
        res[i].resize(a[i].size());
        for (auto j = 0u; j < a[i].size(); ++j){
            res[i][j] = a[i][j] * b;
        }
    }
    return res;
}

template <typename T>
std::vector<vector<T > > operator*(const T b, const std::vector<vector<T> > &a) {
    auto res = std::vector<vector<T > >(a.size());
    for (auto i = 0u; i < a.size(); ++i) {
        res[i].resize(a[i].size());
        for (auto j = 0u; j < a[i].size(); ++j){
            res[i][j] = a[i][j] * b;
        }
    }
    return res;
}

template <typename T>
std::vector<T> operator-(const std::vector<T> &a, const std::vector<T> &b) {
    auto res = std::vector<T>(a.size());
    for (auto i = 0u; i < a.size(); ++i) {
        res[i] = a[i] - b[i];
    }
    return res;
}

template <typename T>
std::vector<T> operator-(const T a, const std::vector<T> &b) {
    auto res = std::vector<T>(b.size());
    for (auto i = 0u; i < b.size(); ++i) {
        res[i] = a - b[i];
    }
    return res;
}

template <typename T>
std::vector<T> operator+(const std::vector<T> &a, const std::vector<T> &b) {
    auto res = std::vector<T>(a.size());
    for (auto i = 0u; i < a.size(); ++i) {
        res[i] = a[i] + b[i];
    }
    return res;
}

template <typename T>
std::vector<std::vector<T > > operator+(const std::vector<std::vector<T > > &a, const std::vector<T> &b) {
    std::vector<std::vector<T > > res;
    res.resize(a.size());
    for (auto i = 0u; i < a.size(); ++i) {
        res[i].resize(b.size());
        for (auto j = 0u; j < b.size(); ++j) {
            res[i][j] = a[i][j] + b[j];
        }
    }
    return res;
}

template <typename T> T norm2(const std::vector<T> &v) {
    auto n = T(0);
    for (const auto el : v) {
        n += el * el;
    }
    return n;
}

template <typename T> bool is_element_in_vector(std::vector<T> v,T element){
    std::vector<int>::iterator it;
    it=find(v.begin(),v.end(),element);
    if (it!=v.end()){
        return true;
    }
    else{
        return false;
    }
}

template <typename T>
void set_matrix_zero(vector<vector<T > > & matrix, int i1, int i2){
    matrix.clear();
    matrix.reserve(i1);
    matrix.resize(i1);
    for (int i=0;i<i1;i++){
        matrix[i].reserve(i2);
        matrix[i].resize(i2);
        for (int j=0;j<i2;j++){
            matrix[i][j]=(T) 0;
        };
    }
}

template <typename T>
void set_vector_zero(vector<T > & vec, int i1){
    vec.clear();
    vec.reserve(i1);
    vec.resize(i1);
    for (int i=0;i<i1;i++){
        vec[i]=(T) 0;
    }
}

template <typename T>
void set_part_matrix_zero(vector<vector<T > > & matrix, int a_stt, int a_final, int b_stt, int b_final){
    for (int i=a_stt;i<a_final;i++){
        for (int j=b_stt;j<b_final;j++){
            matrix[i][j] = (T) 0;
        }
    }
}

template <typename T>
vector<vector<T> > matrix_multiply(vector<vector<T> > arrA, vector<vector<T> > arrB)
{
    //矩阵arrA的行数
    int rowA = arrA.size();
    //矩阵arrA的列数
    int colA = arrA[0].size();
    //矩阵arrB的行数
    int rowB = arrB.size();
    //矩阵arrB的列数
    int colB = arrB[0].size();
    //相乘后的结果矩阵
    vector<vector<T>>  res;
    if (colA != rowB)//如果矩阵arrA的列数不等于矩阵arrB的行数。则返回空
    {
        return res;
    }
    else
    {
        //设置结果矩阵的大小，初始化为为0
        res.resize(rowA);
        for (int i = 0; i < rowA; ++i)
        {
            res[i].resize(colB);
        }

        //矩阵相乘
        for (int i = 0; i < rowA; ++i)
        {
            for (int j = 0; j < colB; ++j)
            {
                for (int k = 0; k < colA; ++k)
                {
                    res[i][j] += arrA[i][k] * arrB[k][j];
                }
            }
        }
    }
    return res;
}

template<typename T>
T vector_multiply(const vector<T > & a, const vector<T > & b){
    if (a.size() != b.size()){
        cout << "Error in Function vector_multiply. a.size() != b.size()\n";
    }
    T result = 0;
    for (int i = 0; i<a.size(); i++){
        result += a[i] * b[i];
    }
    return result;
}

template<typename T>
vector<T > mat_mul_vec(const vector<vector<T > > & mat, const vector<T > & vec){
    vector<T > result;
    int mat_row = mat.size();
    int mat_col = mat[0].size();
    if (mat_col != vec.size()){
        cout << "Error in Function mat_mul_vec. mat[0].size() != vec.size()\n";
    }
    result.clear();
    result.resize(mat_row);
    for (int i=0;i<mat_row;i++){
        for (int j=0;j<mat_col;j++){
            result[i] += mat[i][j] * vec[j];
        }
    }
    return result;
}

// 计算行列式
template<typename _Tp>
_Tp mat_determinant(const std::vector<std::vector<_Tp>>& mat, int N)
{
    if (mat.size() != N) {
        fprintf(stderr, "mat must be square matrix\n");
        return -1;
    }
    for (int i = 0; i < mat.size(); ++i) {
        if (mat[i].size() != N) {
            fprintf(stderr, "mat must be square matrix\n");
            return -1;
        }
    }

    _Tp ret{ 0 };

    if (N == 1) return mat[0][0];

    if (N == 2) {
        return (mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0]);
    }
    else {
        // first col
        for (int i = 0; i < N; ++i) {
            std::vector<std::vector<_Tp>> m(N - 1);
            std::vector<int> m_rows;
            for (int t = 0; t < N; ++t) {
                if (i != t) m_rows.push_back(t);
            }
            for (int x = 0; x < N - 1; ++x) {
                m[x].resize(N - 1);
                for (int y = 0; y < N - 1; ++y) {
                    m[x][y] = mat[m_rows[x]][y + 1];
                }
            }
            int sign = (int)pow(-1, 1 + i + 1);
            ret += mat[i][0] * sign * mat_determinant<_Tp>(m, N - 1);
        }
    }

    return ret;
}

// 计算伴随矩阵
template<typename _Tp>
int mat_adjoint(const std::vector<std::vector<_Tp>>& mat, std::vector<std::vector<_Tp>>& adj, int N)
{
    if (mat.size() != N) {
        fprintf(stderr, "mat must be square matrix\n");
        return -1;
    }
    for (int i = 0; i < mat.size(); ++i) {
        if (mat[i].size() != N) {
            fprintf(stderr, "mat must be square matrix\n");
            return -1;
        }
    }

    adj.resize(N);
    for (int i = 0; i < N; ++i) {
        adj[i].resize(N);
    }

    for (int y = 0; y < N; ++y) {
        std::vector<int> m_cols;
        for (int i = 0; i < N; ++i) {
            if (i != y) m_cols.push_back(i);
        }

        for (int x = 0; x < N; ++x) {
            std::vector<int> m_rows;
            for (int i = 0; i < N; ++i) {
                if (i != x) m_rows.push_back(i);
            }

            std::vector<std::vector<_Tp>> m(N - 1);
            for (int i = 0; i < N - 1; ++i) {
                m[i].resize(N - 1);
            }
            for (int j = 0; j < N - 1; ++j) {
                for (int i = 0; i < N - 1; ++i) {
                    m[j][i] = mat[m_rows[j]][m_cols[i]];
                }
            }

            int sign = (int)pow(-1, x + y);
            adj[y][x] = sign * mat_determinant<_Tp>(m, N-1);
        }
    }

    return 0;
}

template<typename _Tp>
void print_matrix(const std::vector<std::vector<_Tp>>& mat)
{
    int rows = mat.size();
    for (int y = 0; y < rows; ++y) {
        for (int x = 0; x < mat[y].size(); ++x) {
            fprintf(stderr, "  %f  ", mat[y][x]);
        }
        fprintf(stderr, "\n");
    }
    fprintf(stderr, "\n");
}

template<typename _Tp>
void print_matrix_to_file(const std::vector<std::vector<_Tp>>& mat, std::string filename)
{
    ofstream outfile;
    outfile.open(filename);
    outfile << fixed << right;
    int rows = mat.size();
    for (int y = 0; y < rows; ++y) {
        for (int x = 0; x < mat[y].size(); ++x) {
            outfile << setw(18) << setprecision(9) << mat[y][x];
        }
        outfile << "\n";
    }
    outfile.close();
}

// 求逆矩阵
template<typename _Tp>
int mat_inverse(const std::vector<std::vector<_Tp>>& mat, std::vector<std::vector<_Tp>>& inv, int N)
{
    if (mat.size() != N) {
        fprintf(stderr, "mat must be square matrix\n");
        return -1;
    }
    for (int i = 0; i < mat.size(); ++i) {
        if (mat[i].size() != N) {
            fprintf(stderr, "mat must be square matrix\n");
            return -1;
        }
    }

    _Tp det = mat_determinant(mat, N);
    if (fabs(det) < 1.0e-5) {
        fprintf(stderr, "mat's determinant don't equal 0\n");
        return -1;
    }

    inv.resize(N);
    for (int i = 0; i < N; ++i) {
        inv[i].resize(N);
    }

    double coef = 1.f / det;
    std::vector<std::vector<_Tp>> adj;
    if (mat_adjoint(mat, adj, N) != 0) return -1;

    for (int y = 0; y < N; ++y) {
        for (int x = 0; x < N; ++x) {
            inv[y][x] = (_Tp)(coef * adj[y][x]);
        }
    }

    return 0;
}

// transpose
template<typename T>
vector<vector<T > > mat_transpose(const vector<vector<T > > & mat){
    int row = mat.size();
    int col = mat[0].size();
    vector<vector<T > > trans_mat;
    trans_mat.resize(col);
    for (int i = 0; i<col; i++){
        trans_mat[i].resize(row);
        for (int j = 0; j<row; j++){
            trans_mat[i][j] = mat[j][i];
        }
    }
    return trans_mat;
}

#endif //VECTOR_HELPER_H

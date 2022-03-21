//
// Created by п on 24.01.2022.
//

#ifndef NEW_MATRIX_MATRIX2_0_H
#define NEW_MATRIX_MATRIX2_0_H
#include<utility>
#include <iostream>
#include <optional>
#include<cassert>
class matrix {
public:
    static double cin_item();
    matrix(unsigned int row, unsigned int column);
    matrix(unsigned row, unsigned  column,const double *mass);
    matrix(const matrix& matrix);
    matrix& operator= (const matrix& matrix);
    matrix& operator*= (const matrix &matrix2);
    matrix& operator+= (const matrix& matrix2);
    matrix& operator-= (const matrix& matrix2);
    matrix operator* (const matrix &matrix2);
    matrix operator+ (const matrix &matrix2);
    matrix operator- (const matrix &matrix2);
    double& operator() (unsigned int i, unsigned int j);
    double operator() (unsigned int i, unsigned int j) const;
    friend std::ostream& operator<< (std::ostream& os, const matrix& matrix2);
    friend std::istream& operator>> (std::istream &is, matrix& matrix2);
    void transponse();
    void find_det_and_rank();
    void find_rank ();
    void find_det ();
    void inverse();
    void subGauss0();
    matrix& slau_solution(matrix &b_column);
    ~matrix();
private:
    matrix& Gauss_Seidel(const matrix& b);
    void find_inverse();
    double* subGauss();
    unsigned int row;
    unsigned int column;
    double *data = nullptr;
    std::optional<double> det;
    std::optional<unsigned int> rank;
};


#endif //NEW_MATRIX_MATRIX2_0_H

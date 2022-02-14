//
// Created by п on 24.01.2022.
//

#ifndef NEW_MATRIX_MATRIX2_0_H
#define NEW_MATRIX_MATRIX2_0_H
#include<utility>
#include <iostream>
#include <experimental/optional>
#include<cassert>
class matrix {
public:
    static double cin_item();
    matrix(unsigned int row, unsigned int column);
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
    ~matrix();
private:
    unsigned int row;
    unsigned int column;
    double *data = nullptr;
    std::experimental::optional<double> det;
    std::experimental::optional<double> rank;
};


#endif //NEW_MATRIX_MATRIX2_0_H

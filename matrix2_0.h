//
// Created by п on 24.01.2022.
//

#ifndef NEW_MATRIX_MATRIX2_0_H
#define NEW_MATRIX_MATRIX2_0_H
#include<iostream>
#include<utility>
#include<cassert>
class matrix {
public:
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
    ~matrix(){
        delete [] data;
    };
private:
    unsigned int row;
    unsigned int column;
    double *data = nullptr;
    std::pair<bool, unsigned int> rank;
    std::pair<bool, double> det;

};


#endif //NEW_MATRIX_MATRIX2_0_H

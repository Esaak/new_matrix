#include "matrix2_0.h"
matrix:: matrix(unsigned int row, unsigned int column): row(row), column(column), rank(false, 0), det(false, 0){
    assert(column !=0 && row !=0);
    data = new double[row*column];
    for(unsigned int i=0; i<row; i++){
        for(unsigned int j=0; j<column; j++){
            std::cin>>data[j+i*column];
        }
    }
}
/*matrix:: matrix(unsigned int row, unsigned int column,const double *mass): matrix::matrix(row, column){
    for(unsigned int i=0; i<row; i++){
        for(unsigned int j=0; j<column; j++){
            data[j+i*column] = mass[j+i*column];
        }
    }
}*/
matrix::matrix(const matrix &matrix):row(matrix.row), column(matrix.column){
    assert(column !=0 && row !=0);
    data = new double[row*column];
    if(matrix.rank.first){
        rank = matrix.rank;
    }
    else{
        rank.second = false;
    }
    if(matrix.det.first){
        det = matrix.det;
    }
    else{
        det.second = false;
    }
    for (unsigned int i = 0; i < row; i++) {
        for (unsigned int j = 0; j < column; j++) {
            data[j+i*column] =matrix.data[j+i*column];
        }
    }
}

matrix& matrix::operator= (const matrix& matrix){
    if (this == &matrix) {
        return *this;
    }
    delete[] data;
    data = nullptr;
    column = matrix.column;
    row = matrix.row;
    data = new double [column*row];
    for (unsigned int i = 0; i <matrix.row; i++) {
        for (unsigned int j = 0; j < matrix.column; j++) {
            this->data[j+i*column] = matrix.data[j+i*column];
        }
    }
    if(matrix.det.first){
        this->det = matrix.det;
    }
    if(matrix.rank.first){
        this->rank = matrix.rank;
    }
    return *this;
}
matrix& matrix::operator*= (const matrix &matrix2) {
    assert(column==matrix2.row && row==matrix2.column);
    double *new_data = new double [row*matrix2.column];
    for(unsigned int i=0; i<row; i++){
        for(unsigned int j=0; j<matrix2.column; j++){
            new_data[j+i*matrix2.column]=0;
            for(unsigned int k=0; k<column; k++){
                new_data[j+i*matrix2.column] +=data[k+i*column]*matrix2.data[j+k*column];
            }
        }
    }
    delete [] data;
    data = new_data;
    return *this;
}
matrix& matrix::operator+= (const matrix &matrix2) {
    assert(column==matrix2.column && row==matrix2.row);
    for(unsigned int i=0; i<row; i++){
        for(unsigned int j=0; j<column; j++){
            data[j+i*column]+=matrix2.data[j+i*column];
        }
    }
    return *this;
}
matrix& matrix::operator-= (const matrix &matrix2) {
    assert(column==matrix2.column && row==matrix2.row);
    for(unsigned int i=0; i<row; i++){
        for(unsigned int j=0; j<column; j++){
            data[j+i*column]-=matrix2.data[j+i*column];
        }
    }
    return *this;
}
matrix matrix::operator+ (matrix const &matrix2){
    matrix new_matrix =*this;
    return new_matrix+=matrix2;
}


matrix matrix::operator- (matrix const &matrix2){
    matrix new_matrix = *this;
    return new_matrix-=matrix2;
}
matrix matrix::operator* (const matrix &matrix2){
    matrix new_matrix=*this;
    return new_matrix*=matrix2;
}
std:: ostream& operator<< (std::ostream& os ,const matrix &matrix2) {
    os<<std::fixed;
    os.precision(3);
    os<<" Your matrix : \n";
    for(unsigned int i=0; i<matrix2.row; i++){
        for(unsigned int j=0; j<matrix2.column; j++){
            std::cout<<matrix2.data[i*matrix2.column + j]<<" ";
        }
        std::cout<<'\n';
    }
    if(matrix2.det.first){
        os<<"Matrix determinant = "<< matrix2.det.second<<"\n";
    }
    if(matrix2.rank.first){
        os<<"Matrix rank = "<< matrix2.rank.second<<"\n";
    }
    os<< "Matrix rows = "<< matrix2.row<<"\n";
    os<<"Matrix columns = "<< matrix2.column<<"\n";
    return os;
}
std:: istream& operator>> (std::istream &is, matrix& matrix2){
    for(unsigned int i=0; i<matrix2.row; i++){
        for(unsigned int j=0; j<matrix2.column; j++){
            std::cin>>matrix2.data[i+j*matrix2.column];
        }
    }
    return is;
}
double& matrix:: operator() (unsigned int i, unsigned int j){
    return data[i+j*column];
}
double matrix:: operator() (unsigned int i, unsigned int j) const {
    return data[i+j*column];
}
#include "matrix2_0.h"
matrix::~matrix() {
    delete [] data;
}

double matrix:: cin_item(){
    double d;
    std::cin>>d;
    //std::cin>>d;
    return d;
}
matrix:: matrix(unsigned int row, unsigned int column):row(row), column(column) {
    assert(column !=0 && row !=0);
    data = new double[row*column];
    for(unsigned int i=0; i<row; i++){
        for(unsigned int j=0; j<column; j++){
            data[j+i*column] = cin_item();
        }
    }
}
matrix:: matrix(unsigned int row, unsigned  int column, const double *mass): row(row), column(column){
    data = new double [row*column];
    for(unsigned i=0; i<row; i++){
        for(unsigned j=0; j<column; j++) {
            data[i * column + j] = mass[i * column + j];
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
    rank = matrix.rank;
    det = matrix.det;
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
    this->det = matrix.det;
    this->rank = matrix.rank;
    return *this;
}
/*matrix& matrix::operator*= (const matrix &matrix2) {
    assert(column==matrix2.row);
    double *new_data = new double [row*matrix2.column];
    for(unsigned int i=0; i<row; i++){
        for(unsigned int j=0; j<matrix2.column; j++){
            new_data[j+i*matrix2.column]=0;//сильно не уверен, проверить
            for(unsigned int k=0; k<column; k++){
                new_data[j+i*matrix2.column] +=data[k+i*column]*matrix2.data[j+k*column];
            }
        }
    }
    delete [] data;
    data = new_data;
    return *this;
}*/
matrix& matrix::operator*= (const matrix &matrix2) {
    assert(column==matrix2.row);
    double *new_data = new double [row*matrix2.column]{};
    for(unsigned int i=0; i<row; i++){
        for(unsigned int k=0; k<matrix2.column; k++){//сильно не уверен, проверить
            for(unsigned int j=0; j<column; j++){
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
    std::streamsize ss = std::cout.precision();
    os<<std::fixed;
    os.precision(3);
    os<<" Your matrix : \n";
    for(unsigned int i=0; i<matrix2.row; i++){
        for(unsigned int j=0; j<matrix2.column; j++){
            std::cout<<matrix2.data[i*matrix2.column + j]<<" ";
        }
        std::cout<<'\n';
    }
    if(matrix2.det){
        os<<"Matrix determinant = "<< matrix2.det.value()<<"\n";
    }
    if(matrix2.rank){
        os<<"Matrix rank = "<< matrix2.rank.value()<<"\n";
    }
    os<< "Matrix rows = "<< matrix2.row<<"\n";
    os<<"Matrix columns = "<< matrix2.column<<"\n";
    os.precision(ss);
    return os;
}
std:: istream& operator>> (std::istream &is, matrix& matrix2){
    for(unsigned int i=0; i<matrix2.row; i++){
        for(unsigned int j=0; j<matrix2.column; j++){
            is>>matrix2.data[i+j*matrix2.column];
        }
    }
    return is;
}
double& matrix:: operator() (unsigned int i, unsigned int j){
    return data[j+i*column];
}
double matrix:: operator() (unsigned int i, unsigned int j) const {
    return data[j+i*column];
}

double* matrix:: subGauss() {
    unsigned int rank_temp = column;
    double *new_data = new double [column*row];
    // stupid operation
    for(unsigned i=0; i<row; i++){
        for(unsigned j=0; j<column; j++){
            new_data[i*column+ j]= data[i*column+ j];
        }
    }
    for (unsigned row_counter = 0; row_counter < rank_temp; row_counter++)
    {
        if (new_data[row_counter+row_counter*column]!=0){
            for (unsigned column_counter = row_counter+1; column_counter < row; column_counter++)
            {
                double mult = (double)new_data[column_counter*column+row_counter] /(double)new_data[row_counter*column+ row_counter];
                //std::cout<<mult<<'\n';
                for (unsigned i = row_counter; i < rank_temp; i++)//0->row_counter
                    new_data[column_counter*column +i] -= mult * new_data[row_counter*column+i];
            }
            /*for(unsigned i=0; i<row;i++){
                for(unsigned j=0; j<column; j++){
                    std::cout<<new_data[j+i*column]<<" ";
                }
                std::cout<<"\n";
            }*/
        }
        else
        {
            bool reduce = true;
            for (unsigned i = row_counter + 1; i < row; i++)
            {
                if (new_data[i*column+row_counter]!=0)
                {
                    //stupid swap
                    for (unsigned j = 0; j < rank_temp; j++)
                    {
                        double temp = new_data[row_counter*column+j];
                        new_data[row_counter*column+j] = new_data[i*column+j];
                        new_data[i*column+j] = temp;
                    }
                    reduce = false;
                    break ;
                }
            }
            if (reduce)
            {
                rank_temp--;
                for (unsigned i = 0; i < row; i ++)
                    new_data[i*column+row_counter] = new_data[i*column+rank_temp];
            }
            row_counter--;
        }
    }
    rank = rank_temp;
    return new_data;
}
void matrix::subGauss0() { //for rectangular
    unsigned rank_temp = column;
    bool *basis = new bool[column]{};
    for (unsigned row_counter = 0, row_counter_column=0; row_counter < rank_temp, row_counter_column<rank_temp; row_counter_column++, row_counter++)
    {

        if (data[row_counter_column+row_counter*column]!=0)
        {
            basis[row_counter_column] =true;
            for (unsigned column_counter =0; column_counter < row; column_counter++)
            {
                if(column_counter!= row_counter) {
                    //double mult = (double)new_data[column_counter*column+row] /
                    //             new_data[row_counter*column+ row_counter];te
                    double tempa = data[row_counter * column + row_counter_column];
                    double tempb = data[column_counter * column + row_counter_column];
                    for (unsigned i = 0; i < column; i++) {//may be i = 0, i<rank_temp
                        data[row_counter * column + i] /= tempa;
                        data[column_counter * column + i] -=
                                tempb * data[row_counter * column + i];
                    }
                }
            }
            /*for(unsigned i=0; i<row; i++){
                for(unsigned  j=0; j<column; j++){
                    //data[i*column+j] = new_data[i*column+j];
                    //std::cout<<"1\t";
                    std::cout<<data[i*column+j]<<' ';
                }
                std::cout<<"\n";
            }*/
            if(row_counter == row-1)
                break;
        }
        else
        {
            bool reduce = true;
            for (unsigned i = row_counter + 1; i < row; i++)
            {
                if (data[i*column+row_counter_column]!=0)
                {
                    //stupid swap
                    for (unsigned j = 0; j < column; j++)//j=0, j<rank_temp;
                    {
                        double temp1 = data[row_counter*column+j];
                        data[row_counter*column+j] = data[i*column+j];
                        data[i*column+j] = temp1;
                    }
                    reduce = false;
                    break ;
                }
            }
            if (reduce)
            {
                rank_temp--;
            }
            row_counter--;//row_counter --
        }

    }
    /*for(unsigned i=0; i<row; i++){
        for(unsigned  j=0; j<column; j++){
            data[i*column+j] = new_data[i*column+j];
            //std::cout<<"1\t";
            //std::cout<<data[i*column+j]<<' ';
        }
    }*/
}

/*
matrix& matrix::Gauss_Seidel(const matrix& b){
    //stupid operation
    double *a_big = new double[(column+1)*row];
    for(unsigned i=0; i<row; i++){
        for(unsigned j=0; j<column +1; j++) {
            if(j == column){
                a_big[j+i*(column+1)] = b(i,j);
            }
            else{
                a_big[j + i * (column + 1)] = b(i, j);
            }
        }
    }
    //stupid operation end;
    matrix solution(row, 1);//в конструкторе матриц есть cin_item




}

*/


void matrix:: transponse(){
    assert(column!=0 && row!=0);
    double* new_data = new double;
    if(column == row){
        for(unsigned int i=0; i<row; i++){
            for(unsigned int j=1+i; j<row; j++){
                double tmp = data[i*row + j];
                new_data[i*row+j] = data[j*row+i];
                new_data[j*row+i] = tmp;
            }
        }
        for(unsigned int i=0; i<column; i++){
            new_data[i*column + i]=data[i*column + i];
        }
        delete [] data;
        data = new_data;
        return;
    }
    for(unsigned int t=0; t<column; t++){
        for(unsigned int q=0;q<row; q++){
            new_data[t*row + q] =data[q*column + t];
        }
    }
    delete [] data;
    data = new_data;
}

void matrix:: find_det_and_rank(){
    double *new_data = subGauss();
    if(column == row){
        if(rank.value() == row) {
            double det_temp = 1;
            for (unsigned j = 0; j < row; j++) {
                det_temp *= new_data[j * column + j];
            }
            det = det_temp;
        } else{
            det =0;
        }
    }
    //for(unsigned i=0; i<row; i++){
      //  for(unsigned j=0; j<column; j++){
        //    std::cout<<new_data[row*i+j]<<" ";
        //}
        //std::cout<<"\n";
    //}
    delete [] new_data;
}
//I don't know, but it is code duplication;
void matrix:: find_rank() {
    find_det_and_rank();
    det.reset();
}

void matrix::find_det() {
    find_det_and_rank();
    rank.reset();
}

void matrix::find_inverse() {
    double *new_data = new double[column*row]{};
    //create identity matrix
    for(unsigned i=0; i<row; i++){
        new_data[row*i+i] = 1;
    }
    unsigned rank_temp = column;
    for (unsigned row_counter = 0, row_counter_column=0; row_counter < rank_temp, row_counter_column<rank_temp; row_counter_column++, row_counter++)
    {

        if (data[row_counter_column+row_counter*column]!=0)
        {
            for (unsigned column_counter =0; column_counter < row; column_counter++)
            {
                if(column_counter!= row_counter) {
                    //double mult = (double)new_data[column_counter*column+row] /
                    //             new_data[row_counter*column+ row_counter];te
                    double tempa = data[row_counter * column + row_counter_column];
                    double tempb = data[column_counter * column + row_counter_column];
                    for (unsigned i = 0; i < column; i++) {//may be i = 0, i<rank_temp
                        data[row_counter * column + i] /= tempa;
                        new_data[row_counter * column + i] /= tempa;
                        data[column_counter * column + i] -=
                                 tempb * data[row_counter * column + i];
                        new_data[column_counter * column + i] -=
                                tempb * new_data[row_counter * column + i];
                    }
                }
            }
            /*for(unsigned i=0; i<row; i++){
                for(unsigned  j=0; j<column; j++){
                    //data[i*column+j] = new_data[i*column+j];
                    //std::cout<<"1\t";
                    std::cout<<new_data[i*column+j]<<' ';
                }
                std::cout<<"\n";
            }*/
        }
        else
        {
            bool reduce = true;
            for (unsigned i = row_counter + 1; i < row; i++)
            {
                if (data[i*column+row_counter_column]!=0)
                {
                    //stupid swap
                    for (unsigned j = 0; j < column; j++)//j=0, j<rank_temp;
                    {
                        double temp = new_data[row_counter*column+j];
                        double temp1 = data[row_counter*column+j];
                        new_data[row_counter*column+j] = new_data[i*column+j];
                        new_data[i*column+j] = temp;
                        data[row_counter*column+j] = data[i*column+j];
                        data[i*column+j] = temp1;
                    }
                    reduce = false;
                    break ;
                }
            }
            if (reduce)
            {
                rank_temp--;
            }
            row_counter--;//row_counter --
        }

    }

    /*for(unsigned i=0; i<row; i++){
        for(unsigned  j=0; j<column; j++){
            data[i*column+j] = new_data[i*column+j];
            //std::cout<<"1\t";
            //std::cout<<data[i*column+j]<<' ';
        }
    }*/
    data = new_data;
    delete [] new_data;
}

void matrix:: inverse(){
    assert(column==row);
    if(det.has_value() ){
        if(det.value_or(0)){
            find_inverse();
        }
        else{
            std::cout<<"This matrix don't have inverse matrix, sorry gay";
        }
    }
    else{
        find_det();
        if(det.value_or(0)){
            find_inverse();
        }
        else{
            std::cout<<"This matrix don't have inverse matrix, sorry gay";
        }
        det.reset();
    }

}
matrix& matrix::slau_solution(matrix &b_column) {
    assert(b_column.column  == row && row<=  column );
    if(row== column) {
        //stupid operation
        matrix *new_matrix = new matrix(1,1);
        new_matrix = this;
        new_matrix->inverse();
        *new_matrix *=b_column;
        return *new_matrix;
    }
    else {
        //stupid operations
        double *new_data = new double[(column + 1) * row];
        for (unsigned i = 0; i < row; i++) {
            for (unsigned j = 0; j < column + 1; j++) {
                if (j == column) {
                    new_data[i * column + column] = b_column.data[i];
                } else {
                    new_data[i * column + j] = data[i * column + j];
                }
            }
        }
        matrix a_big(row, column+1, new_data);//такого конструктора нет и он тупой
        a_big.subGauss0();
        //unsigned i=0;
        matrix FSR(column, column+1-a_big.rank());//нулями заполненная матрица
        double *priv_solution = new double [column];
        for(unsigned i=0; i < row; i++){
          for(unsigned j=i; j < column+1; j++){// мысль такая, найти все бахисные столбцы в Subgauss.

          }
        }

    }
    //stupid operation end

}


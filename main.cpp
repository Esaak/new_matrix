#include "matrix2_0.h"
#include<cstdlib>
#include<ctime>
int main() {
    /*double a[9];
    double a1[9];
    int k=2;
    for(int i=0; i<3;i++){
        for(int j=0; j<3;j++){
            a[j+i*3] =k;
            k++;
            a1[j+i*3]=k++;
        }
    }
    //matrix b(1000,1000);
    //matrix c(1000,1000);
    matrix d=b+c;
    std::cout<<1;
    std::cout<<d;
    std::cout<<d(0,0);
    d(0,0) +=1;
    std::cout<<d;
    std::cin>>d;
    std::cin>>d(0,0);
    std::cout<<d;
    //b*=c;
    //std::cout<<b;
    //b.transponse();
    //std::cout<<"\n";
    //std::cout<<b;

    //d.show_matrix();
    matrix mat(3,3);
    //mat.find_det_and_rank();
    mat.inverse();
    std::cout<<mat;*/
    /*matrix a(3,4);
    a.subGauss0();
    std::cout<<a;*/
    matrix a(3,4);
    matrix b(3,1);
    std::pair<matrix, double*> wq =a.slau_solution(b) ;
    std::cout<<wq.first;
    std::cout<<wq.second;
    return 0;
}

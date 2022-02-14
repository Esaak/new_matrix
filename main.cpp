#include "matrix2_0.h"
#include<cstdlib>
#include<ctime>
int main() {
    double a[9];
    double a1[9];
    int k=2;
    for(int i=0; i<3;i++){
        for(int j=0; j<3;j++){
            a[j+i*3] =k;
            k++;
            a1[j+i*3]=k++;
        }
    }
    double time_begin = clock();
    matrix b(1000,1000);
    matrix c(1000,1000);
    /*matrix d=b+c;
    std::cout<<1;
    std::cout<<d;
    std::cout<<d(0,0);
    d(0,0) +=1;
    std::cout<<d;
    std::cin>>d;
    std::cin>>d(0,0);
    std::cout<<d;*/
    b*=c;
    double time_end = clock();
    std:: cout<<'\n'<< time_end-time_begin;
    //d.show_matrix();
    return 0;
}

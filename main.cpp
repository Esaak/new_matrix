#include <iostream>
#include "matrix2_0.h"
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
    matrix b(3,3);
    matrix c(3,3);
    matrix d=b+c;
    std::cout<<1;
    std::cout<<d;
    std::cout<<d(0,0);
    d(0,0) +=1;
    std::cout<<d;
    std::cin>>d;
    std::cin>>d(0,0);
    std::cout<<d;
    //d.show_matrix();
    return 0;
}

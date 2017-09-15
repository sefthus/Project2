#include <iostream>
#include "Armadillo"
#include "the_jacobi_method.hpp"
using namespace arma;
using namespace std;


int main(){
    int n = 100;
    double omega_r=0.25;
    bool eigtest = false;
    bool coloumb = true;
    mat A=zeros<mat>(n,n);
    mat R; R.eye(n,n);

    Jacobi_method(A,R,omega_r,n,eigtest,coloumb);
    //testEig();
    //testOffdiag();
}

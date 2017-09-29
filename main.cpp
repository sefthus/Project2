#include <iostream>
#include "armadillo"
#include "the_jacobi_method.hpp"
using namespace arma;
using namespace std;


int main(){
    int n = 170;

    // omega value for one non-interacting electron:
    double omega_r= 1.0;

    // different omega values for interacting case:
    //double omega_r= 0.01;
    //double omega_r= 0.5;
    //double omega_r= 1.0;
    //double omega_r= 5.0;

    //test omega to compare with article:
    //double omega_r= 0.25;

    // choose testing or interacting electrons:
    bool eigtest = false;
    bool coloumb = false;

    mat A=zeros<mat>(n,n);
    mat R; R.eye(n,n);

    Jacobi_method(A,R,omega_r,n,eigtest,coloumb);
    //testEig();
    //testOffdiag();
}

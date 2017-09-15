#include <iostream>
#include "Armadillo"
#include "the_jacobi_method.hpp"
using namespace arma;
using namespace std;


int main(){
    int n = 5;
    bool eigtest = false;
    mat A=zeros<mat>(n,n);
    mat R; R.eye(n,n);

    Jacobi_method(A,R,n,eigtest);
    //testEig();
    //testOffdiag();
}

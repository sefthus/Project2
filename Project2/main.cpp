#include <iostream>
#include "armadillo"


using namespace arma;
using namespace std;

int main()
{
    //setting up empty A matrix
    mat A = zeros<mat>(n-1,n-1);

    //Empty potential
    double V;
    double rho;

    //step size
    double h=(rho_max-rho_min)/(double) n;

    //matrix elements
    double d = 2./(h*h);
    double e = -1./(h*h);


    for (int i=0;i<(n-1);i++){
        rho = (i+1)*h;
        V = rho*rho;

        //setting diagonal elements
        A(i,i)=d+V;
        if(i,i+1){
            A(i,i+1)= A(i+1,i)=e;
        }
    }
    return A;



}

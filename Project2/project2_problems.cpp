#include <iostream>
#include "armadillo"

using namespace arma;
using namespace std;

mat makeAmatrix(double rho_max, double rho_min, int n)
{
    //setting up empty A matrix
    mat A = zeros<mat>(n,n);

    //Potential and varible rho
    double V;
    double rho;

    //step size
    double h=(rho_max-rho_min)/(double) (n+1);

    //matrix elements
    double d = 2./(h*h);
    double e = -1./(h*h);


    for (int i=0;i<(n);i++){
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
// the off diagonal elements, offdiag from Armadillo
double offdiag(mat A, int p, int q, int n){
    double max;
    for (int i=0; i<n;++i){
        for (int j = i+1; j<n;++j){
            double aij = fabs(A(i,j));
            if (aij > max){
                max = aij;
                p =i;
                q=j;
            }
        }

    }
    return max;
}

void Jacobi_rotate( mat A, mat R, int k, int l, int n){
    // s is sine, c is cosine, t is tangent, tau is cot2theta
    double s,c;
    if (A(k,l) !=0.){
        double t;
        double tau = (A(l,l)-A(k,k))/(2*A(k,l));
        if tau>=0{
            t = 1./(tau+sqrt(1. + tau*tau));
        }
        else{
            t= -1./(-tau+sqrt(1. + tau*tau));
        }
        c = 1./sqrt(1+t*t); //cosine
        s = c*t;            //sinus
        }
    else{
    c=1.;
    s=0.;
    }

    double a_kk = A(k,k);
    double a_ll = A(l,l);
    A(k,k) = a_kk*c*c - 2*c*s*A(k,l) + a_ll*s*s;
    A(l,l) = a_ll*c*c + 2*c*s*A(k,l) + a_kk*s*s;
    A(l,k) = 0.; //hardcoding non-diagonal elements
    A(k,l) = 0.;
    for (int i =0;i<n;i++){
        if (i !=k && i!=1){
            double a_ik = A(i,k);
            double a_il = A(i,l);
            A(i,k) = a_ik*c - a_il*s;
            A(k,i) = A(i,k);
            A(i,l) = a_il*c + a_ik*s;
            A(l,i) = A(i,l);
        }
        // new eigenvectors
        double r_ik = R(i,k);
        double r_il = R(i,l);
        R(i,k) = c*r_ik - s*r_il;
        R(i,l) = c*r_il + s*r_ik;
    }

}


















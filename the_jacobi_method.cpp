#include <iostream>
#include <iomanip>
#include <cmath>
#include "armadillo"
#include "time.h"
#include <chrono>
#include <fstream>

using namespace arma;
using namespace std;

void makeAmatrix(mat &A,double rho_min, double rho_max, double omega_r,int n,int coloumb)
{
    //setting up empty A matrix
    //mat A = zeros<mat>(n,n);

    //Potential and varible rho
    vec V(n);
    //double V;
    double rho;
    //step size
    double h=(rho_max-rho_min)/(n+1);

    //matrix elements
    double d = 2./(h*h);
    double e = -1./(h*h);


    for (int i=0;i<n;i++){
        rho = (i+1)*h;
        if (coloumb==1){
            V(i)= omega_r*omega_r*rho*rho + 1./rho;}
        else{
            V(i) = omega_r*omega_r*rho*rho;}




       /* //setting diagonal elements
        A(i,i)=d+V;
        // setting off diagonal elements
        if(i<n-1){
            A(i,i+1)= A(i+1,i)=e;
        }*/
    }
    //no if test here
   A(0,0)= d + V(0);
    A(0,1) = e;
    for (int i=1;i<n-1;i++){
        A(i,i-1) = A(i,i+1) =e;
        A(i,i) = d+V(i);
    }
    A(n-1,n-2)=e;
    A(n-1,n-1)= d + V(n-1);
    //return A;
}
// the off diagonal elements, offdiag from Armadillo
double offdiag(mat &A, int &k, int &l, int n){
    double max=0;
    for (int i=0; i<n;i++){
        for (int j = i+1; j<n;j++){
            double aij = fabs(A(i,j));
            if (aij > max){
                max = aij;
                k = i;
                l = j;
            }
        }

    }

    return max;
}

void Jacobi_rotate( mat &A, mat &R, int &k, int &l, int n){
    // s is sine, c is cosine, t is tangent, tau is cot2theta
    double s,c;
    if (A(k,l) !=0.){
        double t;
        double tau = (A(l,l)-A(k,k))/(2*A(k,l));
        if( tau>=0){
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
        if (i !=k && i!=l){
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
    //cout<<"AAAAAA: "<<A<<endl;
}
void output(double rho_min , double rho_max, double omega_r, mat &R, int n, vec &lambda,int ground_state)
{
    ofstream ofile;
    ofile.open("interacting_5.txt");
    ofile<<setw(15)<<" omega_r" <<setw(15)<<" rho_min "<<setw(15)<<" rho_max "<<setw(20)
        <<" Number of steps "<<setw(15)<<" ground state energy "<<endl;
    ofile<<setw(15) << setprecision(8) << omega_r <<setw(10) << setprecision(8)
        << rho_min<< setw(15) << setprecision(8) << rho_max << setw(15) << n <<setw(25) << lambda[0] <<endl;
    /*ofile <<"omega_r" << setw(15) << setprecision(8) << omega_r << endl;
    ofile<<"rho_min = " << setw(15) << setprecision(8) << rho_min << endl;
    ofile <<"rho_max = " << setw(15) << setprecision(8) << rho_max << endl;
    ofile <<"Number of steps = " << setw(15) << n << endl;
    ofile <<"ground state energy ="<<setw(15) << lambda[0] << endl;*/
    ofile <<"Eigenvector corresponding to lowest eigenvalue:"<<endl;
    for (int i=0;i<n;i++){
        ofile << R(i,ground_state)<<endl;
    }
    ofile<<endl;
    ofile.close();
}

void Jacobi_method(mat &A,mat &R,double omega_r, int n,int eigtest, int coloumb){

    double rho_min = 0.;
    double rho_max = 2.;
    // make A matrix
    if (eigtest==false){
        makeAmatrix(A,rho_min,rho_max,omega_r,n,coloumb);
    }
    // tolerance
    double eps = 1e-8;

    int k = 0;
    int l = 0;
    // get max off diagonal element
    double max_offdiag = offdiag(A,k,l,n);
    int counter=0;
    int maxcount=n*n*n;

    // find eigenvalues of A using armadillo (timed)
    clock_t start1, finish1;
    start1 = clock();

    vec eigval(n);
    eig_sym(eigval,A);
    eigval =sort(eigval);

    finish1=clock();
    double timeused1 = (double)(finish1-start1)/(CLOCKS_PER_SEC);
    cout <<"time used to run armadillo eigval: "<<timeused1<<endl;

    clock_t start2, finish2;
    start2 = clock();

    // find eigenvalues of A using Jacobi's method (timed)

    while (max_offdiag>eps && counter<=maxcount){
        Jacobi_rotate(A,R,k,l,n);
        max_offdiag=offdiag(A,k,l,n);
        counter++;

    }

    finish2=clock();
    double timeused2 = (double)(finish2-start2)/(CLOCKS_PER_SEC);
    cout <<"time used to run Jacobi: "<<timeused2<<endl;

    cout <<"max off diag "<<max_offdiag<<endl;
    //cout<<"A matrix: "<<A<<endl;
    cout <<"Number of similarity transformations needed:"<<counter<<endl;
    cout<<maxcount<<endl;

    // vector of eigenvalues, sorted
    vec lambda = A.diag();
    int ground_state= lambda.index_min();
    lambda=sort(lambda);

    // compare eigenvalues
    cout.precision(5);
    cout<<"Three lowest eigenvalues from armadillo: "<<endl;
    cout<<eigval[0]<<endl;
    cout<<eigval[1]<<endl;
    cout<<eigval[2]<<endl;
    cout<<"Three lowest eigenvalues from Jacobi: "<<endl;
    cout<<lambda[0]<<endl;
    cout<<lambda[1]<<endl;
    cout<<lambda[2]<<endl;
    cout<<"Three lowest known eigenvalues (non-interacting): "<<endl;
    cout<<"3.0000"<<endl;
    cout<<"7.0000"<<endl;
    cout<<"11.000"<<endl;

    output(rho_min,rho_max,omega_r,R,n,lambda,ground_state);


}

/* Test functions */

void testEig(){
    cout<<"Test of Jacobi mathod"<<endl;
    int n=3;
    double omega_r=0;
    bool coloumb=false;
    mat A;
    mat R; R.eye(n,n);
    bool eigtest = true;

    A<< 0 << 1 << 1 << endr
     << 1 << 0 << 2 << endr
     << 1 << 2 << 0 << endr;
    vec eigval(3);
    eig_sym(eigval,A);
    cout <<"Matrix A = "<<A<<endl;
    cout <<"Known eigenvalues: "<< eigval<<endl;
    Jacobi_method(A,R,omega_r, n,eigtest,coloumb);
}

void testOffdiag(){
    int k = 0;
    int l = 0;
    int n = 3;
    mat A;
    A<< 0 << 1 << 8 << endr
     << 1 << 0 << 4 << endr
     << 8 << 4 << 7 << endr;
    double max= offdiag(A,k,l,n);
    cout<<"Test of Offdiag function"<<endl;
    cout<<"Matrix A="<<endl;
    cout<<A<<endl;
    cout<<"Biggest element: "<<max<<endl;
}














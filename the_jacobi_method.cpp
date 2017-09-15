#include <iostream>
#include <iomanip>
#include "armadillo"

using namespace arma;
using namespace std;

void makeAmatrix(mat &A,double rho_min, double rho_max, int n)
{
    //setting up empty A matrix
    //mat A = zeros<mat>(n,n);

    //Potential and varible rho
    //vec V(n);
    double V;
    double rho;

    //step size
    double h=(rho_max-rho_min)/(n+1);

    //matrix elements
    double d = 2./(h*h);
    double e = -1./(h*h);


    for (int i=0;i<n;i++){
        rho = (i+1)*h;
        V = rho*rho;

        //setting diagonal elements
        A(i,i)=d+V;
        // setting off diagonal elements
        if(i<n-1){
            A(i,i+1)= A(i+1,i)=e;
        }
    }
    //no if test here
   /*A(0,0)= d + V(0);
    A(0,1) = e;
    for (int i=1;i<n-1;i++){
        A(i,i-1) = A(i,i+1) =e;
        A(i,i) = d+V(i);
    }
    A(n-1,n-2)=e;
    A(n-1,n-1)= d + V(n-1);*/
    cout<<"matrix A"<<A<<endl;
    //return A;
}
// the off diagonal elements, offdiag from Armadillo
double offdiag(mat &A, int &k, int &l, int n){
    double max;
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
void output(double rho_min , double rho_max, int n, vec& d)
{
    int i;
    cout << "RESULTS:" << endl;
    cout << setiosflags(ios::showpoint | ios::uppercase);
    cout <<"rho_min = " << setw(15) << setprecision(8) << rho_min << endl;
    cout <<"rho_max = " << setw(15) << setprecision(8) << rho_max << endl;
    cout <<"Number of steps = " << setw(15) << n << endl;
    cout << "Five lowest eigenvalues:" << endl;

    for(i = 0; i < 3; i++) {
        cout << setw(15) << setprecision(8) << d[i] << endl;
    }
}

void Jacobi_method(mat &A,mat &R,int n,int eigtest){

    double rho_min = 0.;
    double rho_max = 8.;
    // make A matrix
    if (eigtest==false){
        makeAmatrix(A,rho_min,rho_max,n);
    }
    // tolerance
    double eps = 1e-8;

    int k = 0;
    int l = 0;
    // get max off diagonal element
    double max_offdiag = offdiag(A,k,l,n);
    int counter=0;
    int maxcount=n*n*n;

    while (max_offdiag>eps && counter<=maxcount){
        Jacobi_rotate(A,R,k,l,n);
        max_offdiag=offdiag(A,k,l,n);
        if (counter<17){
            cout<<A<<endl;
            cout<<"BB:  "<<max_offdiag<<" k,l "<<k<<l<<endl;
        }

        counter++;

    }
    cout <<"max off diag "<<max_offdiag<<endl;
    cout<<"A matrix: "<<A<<endl;
    cout <<"Number of similarity transformations needed:"<<counter<<endl;
    cout<<maxcount<<endl;
    // vector of eigenvalues, sorted
    vec lambda = A.diag();
    lambda=sort(lambda);

    output(rho_min,rho_max,n,lambda);


}

/* Test functions */

void testEig(){
    cout<<"Test of Jacobi mathod"<<endl;
    int n=3;
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
    Jacobi_method(A,R,n,eigtest);
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














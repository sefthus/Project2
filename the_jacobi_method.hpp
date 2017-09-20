#ifndef THE_JACOBI_METHOD_HPP
#define THE_JACOBI_METHOD_HPP
#include "armadillo"
using namespace arma;
void testEig();
void testOffdiag();
void Jacobi_method(mat &A,mat &R,double omega_r,int n,int eigtest,int coloumb);

#endif // THE_JACOBI_METHOD_HPP

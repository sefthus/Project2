#ifndef THE_JACOBI_METHOD_HPP
#define THE_JACOBI_METHOD_HPP
#include "Armadillo"
using namespace arma;
void testEig();
void testOffdiag();
void Jacobi_method(mat &A,mat &R,int n,int eigtest);

#endif // THE_JACOBI_METHOD_HPP

#include <iostream>
#include <cmath>
#include <typeinfo>
#include <time.h>
#include <stdio.h>
#include <string>
#include "armadillo"


using namespace std;
using namespace arma;


void lu_decomp(int n, double h, double x){

    vec y, w, x_lu;
    mat L, U, P;
    mat A = zeross<mat>(n+1,n+1);



    //Filling the matrix
    A.diag() = 2.;
    A.diag(1) = -1.;
    A.diag(-1) = -1.;



    for(int i=1; i<=n; i++){
        w[i] = h*h*100*exp(-10*x[i]);
    }

    lu(L,U,P,A);

    y = solve(L,w);
    x_lu = solve(U,y);

}

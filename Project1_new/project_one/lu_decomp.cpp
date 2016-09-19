#include <iostream>
#include <cmath>
#include <typeinfo>
#include <time.h>
#include <stdio.h>
#include <string>
#include "armadillo"
#include "lu_decomp.h"
#include "numerical.h"


using namespace std;
using namespace arma;


mat lu_decomp(int n, double h, mat x, mat x_lu){


    mat L, U, P, f, y;

    clock_t start, finish ;

    mat A = zeros<mat>(n+1,n+1);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //Filling the diagonal of the matrix
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    A.diag() += 2.;
    A.diag(1) += -1.;
    A.diag(-1) += -1;

    mat w = zeros<mat>(n+1,1);

    for(int i=1; i<=n; i++){
        w(i) = h*h*100*exp(-10*x(i));
    }

    f = w;

    start = clock();

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //Using the Armadillo-functions to solve our equation:
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    lu(L,U,P,A);

    y = solve(L,w);
    x_lu = solve(U,y);

    finish = clock();

    //~~~~~~~~~~~~~~~~~~~~~~
    //Reshaping the matrix:
    //~~~~~~~~~~~~~~~~~~~~~~

    x_lu.reshape(n+3,1);
    for (int i=n; i>=0; i--){
        x_lu(i+1) = x_lu(i);
    }

    x_lu(0) = 0.;

   printf("Time taken for LU-decomposition to run: %5.9f seconds. \n", ((float)(finish - start)/CLOCKS_PER_SEC ));

   cout << "We did it! Yes we did" << endl;

    return x_lu;
}

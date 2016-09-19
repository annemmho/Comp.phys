
#include <iostream>
#include <time.h>
#include "armadillo"
#include "lu_arma.h"

using namespace std;
using namespace arma;

////////////////////////////////////////////////
// Solving the same problem by LU-decomposition
////////////////////////////////////////////////

mat lu_arma(int n, double h, mat x_mat, mat x3){

    mat w, y, L_inv, U_inv;

    clock_t start , finish ; // declare start and final time start = clock () ;

    // Initializing matrices for the calculations
    mat A = zeros<mat>(n+1,n+1);

    // Filling matrix A
    A.diag() += 2.;
    A.diag(1) += -1.;
    A.diag(-1) += -1.;

    mat f2 = zeros<mat>(n+1,1);

    for (int i=1; i<=n; i++){
        f2(i) = h*h*100*exp(-10*x_mat(i));
    }

    w = f2;

    mat L,U,P;

    // Doing the LU-decomposition and finding x
    // as described in the lecture notes.

    start = clock();

    lu(L,U,P,A);

    y = solve(L, f2);
    x3 = solve(U, y);

    finish = clock();

    // Shifting the array one element
    x3.reshape(n+3,1);
    for (int i=n; i>=0; i--){
        x3(i+1) = x3(i);
    }

    x3(0) = 0.;

    printf ("Elapsed time LU-decomposition: %5.9f seconds.\n", ( (float)( finish - start ) / CLOCKS_PER_SEC ));

    cout << "LU_decomp. success" << endl;
    return x3;
}

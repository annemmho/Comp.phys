#include <iostream>
#include "numerical.h"
#include <time.h>
#include "lu_decomp.h"
#include <cmath>
#include <typeinfo>
#include <stdio.h>
#include <string>
#include "armadillo"




using namespace std;

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Algorithm for solving tridiagonal systems
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


void numerical(int n, double *a, double *b, double *c, double *v, double *b_tilde){
    int i;

    //Time:
    clock_t start, finish ;

    double t[n+1];
    double b_temp = b[1];

    // Starting the timer
    start = clock();



    v[1] = b_tilde[1]/b_temp;

    //~~~~~~~~~~~~~~~~~~~~~
    //Forward substitution
    //~~~~~~~~~~~~~~~~~~~~~

    for(i=2; i<= n; i++){
        t[i] = c[i-1]/b_temp;
        b_temp = b[i] - a[i]*t[i];
        v[i] = (b_tilde[i] - a[i]*v[i-1])/b_temp;
    }



    //~~~~~~~~~~~~~~~~~~~~~~
    //Backward substitution
    //~~~~~~~~~~~~~~~~~~~~~~

    for(i=n-1; i>=1; i--){
        v[i] -= t[i+1]*v[i+1];
    }

    //Ending timer
    finish = clock();


    printf ("Time taken for Substitution to run: %5.9f seconds.\n", ( (float)( finish - start ) / CLOCKS_PER_SEC ));

    cout << "Numerical Substition works! Yay!" << endl;
}


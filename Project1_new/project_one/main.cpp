// Project 1 main program
#include <iostream>
#include <cmath>
#include <typeinfo>
#include <time.h>
#include <stdio.h>
#include <string>
#include "armadillo"
#include "numerical.h"
#include "lu_decomp.h"





using namespace std;
using namespace arma;


//~~~~~~~~~~~~~~~~
//Functions used:
//~~~~~~~~~~~~~~~~

void SaveToFile(string Filename, int n, double *x, double *v, double *v_cf, mat x_lu, double *error);

void ClosedForm(int , double * , double * );

void Error(int n, double *v, double *v_cf, double *error);




int main(){


    int n = 1000;
    double a[n+1], b[n+1], c[n+1], v[n+1], h, x[n+1], b_tilde[n+1], v_cf[n+1], error[n+1];

    string FileName = "LU_numerical_analytic_n1000.m";


    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //Setting boundary conditions and initializing:
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    v[0] = v[n+1] = 0;
    a[1] = c[n] = 0;
    x[0] = 0;
    x[n+1] = 1;
    v_cf[n] = 0;

    h = 1/(float(n)+1);



    for(int i=1; i<=n; i++){
        x[i] = i*h;
        b_tilde[i] = h*h*100*exp(-10*x[i]);
    }


    mat A = zeros<mat>(n,n);

    for(int i=1; i<=n; i++){
        a[i] = -1.;
        b[i] = 2.;
        c[i] = -1.;
    }


    cout << "The following are the results for n = " << n << endl;

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Solving the problem with the substitution method:
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    numerical(n, a, b, c, v, b_tilde);


    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Analytic solution (closed form solution)
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ClosedForm(n, x, v_cf);


    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //Constructing and reshaping matrices:
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    mat x_matrix = zeros<mat>(n+1,1);
    mat x_lu = zeros<mat>(n+1,1);

    for (int i=1; i<=n; i++){
        x_matrix(i) = i*h;
    }
    x_matrix.reshape(n+2,1);
    x_matrix(n+1) = 1.;


    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //Solving using LU-decomposistion:
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    x_lu = lu_decomp(n, h, x_matrix, x_lu);


    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //Finding the maximum error for each value of n:
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Error(n, v, v_cf, error);


    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Saving our results to a MATLAB-file:
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


    SaveToFile(FileName, n, x, v, v_cf, x_lu, error);



    return 0;


}

//~~~~~~~~~~~~~~~~~~~
//Analytic Solution:
//~~~~~~~~~~~~~~~~~~~

void ClosedForm(int n, double *x, double *v_cf){

    for(int i=1; i<=n; i++){
        v_cf[i] = 1 - (1-exp(-10))*x[i] - exp(-10*x[i]);
    }
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Saving results to a MATLAB-fil for plotting:
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


void SaveToFile(string Filename, int n, double *x, double *v, double *v_cf, mat x_lu, double *error){

    ofstream myfile;
    myfile.open(Filename);
    myfile << "x " << "= [";
    for(int i=0; i<n; i++){
        myfile << x[i] << ", ";

    }
    myfile << "];" << endl;

    myfile << "v " << "= [";
    for(int i=0; i<n; i++){
        myfile << v[i] << ", ";

    }
    myfile << "];" << endl;

    myfile << "v_cf " << "= [";
    for(int i=0; i<n; i++){
        myfile << v_cf[i] << ", ";

    }
    myfile << "];" << endl;


    myfile << "x_lu " << "= [";
    for(int i=0; i<n; i++){
        myfile << x_lu[i] << ", ";
    }
    myfile << "];" << endl;


    myfile << "Error " << "= [";
    for(int i=0; i<n; i++){
        myfile << error[i] << ", ";
    }
    myfile << "];" << endl;

    myfile << "plot(x, v)" << endl;
    myfile << "hold('on')" << endl;
    myfile << "plot(x, v_cf)" << endl;
    myfile << "hold('on')" << endl;
    myfile << "plot(x, x_lu)" << endl;
    myfile << "title('Numerical vs. Analytic Solution for n = " <<  n << "');"  << endl;
    myfile << "xlabel('x');" << endl;
    myfile << "ylabel('v(x)');" << endl;
    myfile << "legend('Numerical Substitution', 'Analytic', 'LU-decomposition');" << endl;

    myfile.close();

}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//Function to find the relative error and its maximum:
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void Error(int n, double *v, double *v_cf, double *error){

    double maximum;

    for(int i=1; i<n; i++){
        error[i] = log10(fabs( (v[i] - v_cf[i])/v_cf[i]  ));
        if (fabs(error[i]>error[i-1])){
            maximum = error[i];
        }


    }
    cout << "The maximum error is: " << maximum << endl;

}



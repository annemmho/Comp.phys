// Project 1 main program
#include <iostream>
#include <cmath>
#include <typeinfo>
#include <time.h>
#include <stdio.h>
#include <string>
#include "armadillo"
//#include "arma_solve.h"
#include "numerical.h"




using namespace std;
using namespace arma;

void SaveToFile(string Filename, int n, double *x, double *v, double *v_cf, vec x_lu);

void ClosedForm(int , double * , double * );

void Error(int n, double *v, double *v_cf, double *h);

mat LU(int n, double h, double *x);


int main(){

    int n = 1000;
    double a[n+1], b[n+1], c[n+1], v[n+1], h, x[n+1], b_tilde[n+1], v_cf[n+1];

    string FileName = "numerical_analytic_n1000.m";


    //vec x_lu(n+1);


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

    numerical(n, a, b, c, v, b_tilde);

    ClosedForm(n, x, v_cf);

    vec x_lu;
    //for(int i=1; i<=n; i++){
      //  x_lu(i) = 0;
    //}


    x_lu = LU(n,h,x);




    SaveToFile(FileName, n, x, v, v_cf, x_lu);






    return 0;




}


void ClosedForm(int n, double *x, double *v_cf){

    for(int i=1; i<=n; i++){
        v_cf[i] = 1 - (1-exp(-10))*x[i] - exp(-10*x[i]);
    }
}


mat LU(int n, double h, double *x){
    mat A2 = zeros<mat>(n+1,n+1);

    A2.diag() += 2.;
    A2.diag(1) += -1.;
    A2.diag(-1) += -1.;

    mat L, U, P;
    vec y, x_lu;
    //double b_tilde;
    mat b_tilde = zeros<mat>(n+1,1);


    for(int i=1; i<=n; i++){
        x[i] = i*h;
        b_tilde(i) = h*h*100*exp(-10*x[i]);
    }


    lu(L,U,P,A2);

    y = solve(L,b_tilde);
    x_lu = solve(U,y);

    return x_lu;
}

void SaveToFile(string Filename, int n, double *x, double *v, double *v_cf, vec x_lu){

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


    myfile << "plot(x, v)" << endl;
    myfile << "hold('on')" << endl;
    myfile << "plot(x, v_cf)" << endl;
    myfile << "hold('on')" << endl;
    //myfile << "plot(x, x_lu)" << endl;
    myfile << "title('Numerical vs. Analytic Solution for n = " <<  n << "');"  << endl;
    myfile << "xlabel('x');" << endl;
    myfile << "ylabel('v(x)');" << endl;
    myfile << "legend('Numerical', 'Analytic');" << endl;

    myfile.close();

}




void Error(int n, double *v, double *v_cf, double *h){

    double *epsilon;

    for(int i=1; i>=n; i++){
        epsilon[i] = log10(fabs( (v[i] - v_cf[i])/v_cf[i]  ));

    }

}

#include <iostream>
#include <armadillo>
#include <string>

using namespace std;
using namespace arma;

void findMaximumElementOnNonDiagonal(mat &A, int &k, int &l, double &max_A, int n);
void findSinAndCos(mat &A, int k, int l, double &c, double &s);
void rotateMatrix(mat &A, int k, int l, double c, double s, int n, mat &R);
void printMatlabMatrix(string name, mat &A);
void JacobiRotations(int n, double rho_max, mat &A, uvec, mat &);
void SavingResults(string, string, string, vec, vec);
void AnalyticSolution(int n, double h, vec &Psi_anal);

int main()
{
    //============
    //Defining variable type
    //============
    int n =200;
    mat A(n,n);
    vec rho(n+2);
    vec v(n+2);

    //======================
    // Defining parameters
    //======================
    double rho_max = 8. ;
    double h = rho_max / (n+2);
    double e = -1 / (h*h);
    double h_temp = 2 / (h*h);

    //===================
    // Producing our rho
    //===================
    for(int i=0; i<n+2; i++){
        rho(i) += h*i;
    }

    rho(0) = 0;

    //===============================
    // Vector w/discretized potential
    //===============================
    v = rho % rho;


    //====================================
    // Filling the matrix A with elements:
    //====================================

    for(int i=0; i<n; i++){
        A(i,i) = h_temp + v(i+1);
    }
    A.diag(1) += e;
    A.diag(-1) += e;


    mat A_Arm = A;
    uvec diagonalA_1;
    mat R_1 = eye<mat>(n,n);
    //cout << "The single electron case: " << endl;
    //JacobiRotations(n, rho_max, A, diagonalA_1, R_1);


    //====================================================
    //Finding the analytic solution and saving the result
    //====================================================
    vec Psi_anal;
    AnalyticSolution(n, h, Psi_anal);
    //cout << "Analytic Solution: " << AnalyticSolution << endl;
    SavingResults("anal_sol.m", "rho", "Psi", rho, Psi_anal);




    //=====================================
    //Using Armadillo to solve the problem:
    //=====================================
    mat eigvec1;
    vec eigval1;
    //eig_sym(eigval1, eigvec1, A_Arm);


    //=================================
    //Printing the Armadillo solutions
    //=================================
//    cout << "Armadillo Calculations " << endl;
//    cout << "lambda0: " << eigval1(0) << "   lambda1: " << eigval1(1) << "   lambda2: " << eigval1(2) << endl;


    //=========================================
    // The potential for the two-electron case
    //=========================================

    double omega_r = 1.;
    vec TwoElectronPotential = omega_r*(rho % rho) + 1 /rho;

    //New matrix
    mat Psi(n,n);
    for(int i=0; i<n; i++){
        Psi(i,i) = h_temp + TwoElectronPotential(i+1);
    }
    Psi.diag(1) += e;
    Psi.diag(-1) += e;
    mat Psi_Arma = Psi;
    uvec diagonalPsi_1;
    mat R_2 = eye<mat>(n,n);
//    cout << "The two-electron case" << endl;
//    JacobiRotations(n, rho_max, Psi, diagonalPsi_1, R_2);

//    //cout << "her er jeg" << endl;
//    SavingResults("noen.m", "rho", "Psi", rho, eigvec1.col(2));



}

void SavingResults(string filename, string name_Vec1, string name_Vec2, vec Vec1, vec Vec2){
    //cout << "dritt" << endl;
    ofstream myfile;
    myfile.open(filename);
    myfile << name_Vec1 << "= [";
    for(int i=0; i<int(Vec1.n_rows); i++){
        myfile << Vec1(i) << ", ";
    }
    myfile << "];" << endl;

    myfile << name_Vec2 << "= [";
    for(int i=0; i<int(Vec2.n_rows); i++){
        myfile << Vec2(i) << ", ";
    }

    myfile << "];" << endl;
    myfile << "plot(" << name_Vec1 << " ," << name_Vec2 << ");" << endl;
    myfile.close();
}

//======================
// Rotating the matrix
//======================

void rotateMatrix(mat &A, int k, int l, double c, double s, int n, mat &R) {
    for(int i=0; i<n; i++){
        if(i != k && i != l){
            A(i,i) = A(i,i);
            double a_il = A(i,l);
            double a_ik = A(i,k);
            A(i,k) = c*a_ik - s*a_il;
            A(i,l) = c*a_il + s*a_ik;
            A(k,i) = A(i,k);
            A(l,i) = A(i,l);

        }

        //Finding
        double r_ik = R(i,k);
        double r_il = R(i,l);
        R(i,k) = c*r_ik - s*r_il;
        R(i,l) = c*r_il + s*r_ik;
    }

    double a_kk = A(k,k);
    double a_ll = A(l,l);

    A(k,k) = c*c*a_kk - 2.*c*s*A(k,l) + s*s*a_ll;
    A(l,l) = s*s*a_kk + 2.*c*s*A(k,l) + c*c*a_ll;
    A(k,l) = 0;
    A(l,k) = 0;


}


void JacobiRotations(int n, double rho_max, mat &A, uvec diagonalA, mat &R){

    double epsilon = pow(10,-8);
    //cout << epsilon << endl;

    // Starting w/A(1,2)
    //cout << A << endl;
    int k = 0; int l = 0;
    double max_A = 0.0;
    findMaximumElementOnNonDiagonal(A, k, l, max_A, n);
    //cout << max_A << endl;

    int numberOfIterations = 0;
    int maxIterations = pow(10, 6);

    //==================================================================================
    //While-loop to change the non-diagonal elements to be approximately zero.
    //==================================================================================


    while( fabs(max_A) > epsilon && (double) numberOfIterations < maxIterations){
        double c = 0.0;
        double s = 0.0;
        findSinAndCos(A, k, l, c, s);
        // Rotate
        rotateMatrix(A, k, l, c, s, n, R);
        //Find max off-diagonal matrix element
        findMaximumElementOnNonDiagonal(A, k, l, max_A, n);
        numberOfIterations++;
    }

    //Preparing results
    diagonalA = sort_index(A.diag());

    //Printing

    //printMatlabMatrix("A", A);
    cout << "Number of iterations: " << numberOfIterations << endl;
    //cout << rho_max << endl;
    cout << "lambda0:  " << A.diag()(diagonalA(0)) << "  lambda1:  " << A.diag()(diagonalA(1)) << "  lambda2:  " << A.diag()(diagonalA(2)) << endl;

    //cout << "hallo" << endl;


}

//=================================================================
//Function to find the non-diagonal element with the maximum value.
//=================================================================

void findMaximumElementOnNonDiagonal(mat &A, int &k, int &l, double &max_A, int n) {
    // Find max off-diagonal matrix element
    max_A = 0;
    for(int i=0; i<n; i++){
        //cout << max_A << endl;
        for(int j=i+1; j<n; j++){
            if( fabs(A(i,j)) > max_A){
                k = i;
                l = j;
                max_A = fabs(A(i,j));
            }
        }
    }
}


//==============================================
//Function to find the sinus and cosinus values
//==============================================

void findSinAndCos(mat &A, int k, int l, double &c, double &s) {
    double tau = ( A(l,l) - A(k,k) ) / ( 2 * A(k,l) );
    double t;
    if(tau<0) {
        t = -tau - sqrt(1.0 + tau*tau);
    }
    else {
        t = -tau + sqrt(1.0 + tau*tau);
    }

    c = 1 / sqrt( 1 + t*t);
    s =  c*t;
}

//=========================================================
// Function calculating the analytic solution, to compare.
//=========================================================

void AnalyticSolution(int n, double h, vec &Psi_anal){
    vec r(n+2);
    int l=0;
    for(int i=0; i<n+2; i++){
        r(i) += i*h;
    }
    r(0) = 0;
    Psi_anal = r % exp(- pow(r, 2) / 8 ) % (1 + r / 2);
    //Psi_anal = pow(r, 1) % exp(- pow(r, 2) / (8*(l+1)) ) % ( 1 + r / (2*(l+1) ) );

}


//=====================================================================================
//Function to print the matrix such that it is possible to check the answere in MATLAB
//=====================================================================================

void printMatlabMatrix(string name, mat &A) {
    cout << name << " = [";
    for(int i=0; i<A.n_rows; i++) {
        for(int j=0; j<A.n_cols; j++) {
            cout << A(i,j) << " ";
        }
        cout << "; ";
    }
    cout << "];" << endl;
    cout << "eig("<<name<<")" << endl;
}

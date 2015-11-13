#include <iostream>
#include <cmath>
#include "lib.h"
#include "gauss-laguerre.cpp"

using namespace std;

double integral_function(double theta1, double theta2, double phi1, double phi2, double r1, double r2);

int main()
{
    int N = 20;
    double a = -2.3;
    double b = 2.3;
    double pi = 3.14159;
    double *x = new double [N];
    double *w = new double [N];
    double exact = 5*pi*pi/(16*16);

    gauleg(a, b, x, w, N);
    double integral = 0.;
//   Brute force integration with six for-loops.
    for (int i=0;i<N;i++){
            for (int j = 0;j<N;j++){
               for (int k = 0;k<N;k++){
                   for (int l = 0;l<N;l++){
                       for (int m = 0;m<N;m++){
                           for (int n = 0;n<N;n++){
                               integral += w[i]*w[j]*w[k]*w[l]*w[m]*w[n]*integration_function(x[i],x[j],x[k],x[l],x[m],x[n]);

                           }
                       }
                   }
               }
            }
       }
    cout << "Gauss Legendre: " << "     " << integral << endl;
    cout << "Exact value: " << "     " << exact << endl;

}


#include <iostream>
#include "lib.h"
#include <cmath>

using namespace std;

//  this function defines the function to integrate
double integration_function(double x1, double y1, double z1, double x2, double y2, double z2);
double Spherical_Coord_function(double x1, double x2, double y1, double y2, double z1, double z2, double theta1, double theta2, double phi1, double phi2);


int main()
{
     //  Defining variables, mesh points and weights
         int N = 20;
         double a = -2.3;
         double b = 2.3;
         double pi = 3.14159;
         double *x = new double [N];
         double *w = new double [N];
         double exact = 5*pi*pi/(16*16);

         gauleg(a, b, x, w, N);
         double integral = 0.;
    //   Brute forces integration with six for-loops.
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

double integration_function(double x1, double y1, double z1, double x2, double y2, double z2){
    double alpha = 2.;

 // Computing the different exponentials and the denominator.
    double exp1 = -2*alpha*sqrt(x1*x1 + y1*y1 + z1*z1);
    double exp2 = -2*alpha*sqrt(x2*x2 + y2*y2 + z2*z2);
    double denom = sqrt(pow((x1-x2), 2) + pow((y1-y2), 2) + pow((z1-z2), 2));

 // Checking the denominator to avoid dividing by zero.
    if(denom <pow(10.,-6.)) { return 0;}
        else return exp(exp1 + exp2)/denom;
 }





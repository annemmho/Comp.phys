#include <iostream>


using namespace std;

int main()
{
         int N = 100;
         double a, b;
    // array for integration points and weights using Legendre polynomials
         double *x = new double [N];
         double *w = new double [N];
    //   set up the mesh points and weights
         gauleg(a,b,x,w, N);
    //   evaluate the integral with the Gauss-Legendre method
    //   Note that we initialize the sum
         double int_gauss = 0.;
    //   six-double loops
         for (int i=0;i<N;i++){
                 for (int j = 0;j<N;j++){
                 for (int k = 0;k<N;k++){
                 for (int l = 0;l<N;l++){
                 for (int m = 0;m<N;m++){
                 for (int n = 0;n<N;n++){
            int_gauss+=w[i]*w[j]*w[k]*w[l]*w[m]*w[n]
           *int_function(x[i],x[j],x[k],x[l],x[m],x[n]);
                    }}}}}
            }

}


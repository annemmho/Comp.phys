#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include "lib.h"
#include "time.h"



using namespace std;


//=========================================
// To get the periodic boundary condition.
//=========================================
inline int periodic(int i, int limit, int add){
    return (i+limit+add) % (limit);
}

void Metropolis(int **s_matrix, int spin, double &E, double &M, double *w, long &idum, int &accepted, int **other_spin);

void Initialize(int **, double, double &, double &, long &, bool, int **);

void Expectation_Values(double *average, double norm, double &total_number_accepted, double &Heat_Capacity,
                        double &Susceptibility, double &Mean_M, int L, double &Average_E, double temp);

void Energy_Probability(int **s_matrix, int spin, double &E, int **other_spin);

void Save_Results_Matlab(string Filename, double *Heat_Capacity, double *Susceptibility, double *Mean_M, int N,
                         double *number_mcs, double *total_number_accepted, double *Probability, double *Average_E);


int main()
{
// The parameters:
    int mcs = 10000000;
    int mc = 100000;
    int N = mcs/mc;
    double initial_temp = 1.0;
    int L = 2;
    bool random = true;
    string Filename = "met_is_mc_den.m";

    double average[5];
    double E = 0;
    double M = 0;
    double w[17];
    long idum = -1;


    int **other_spin;
    other_spin = (int **) matrix(L, 2, sizeof(L));
    for(int i=0; i<L; i++){
        other_spin[i][0] = periodic(i, L, -1);
        other_spin[i][1] = periodic(i, L, 1);
    }

    double Probability[N];
    double total_number_accepted[N];
    int accepted;
    double Heat_Capacity[N];
    double Mean_M[N];
    double Susceptibility[N];
    double Average_E[N];

    //To initalize w:
    for(int dE=-8; dE<=8; dE++) w[dE+8] = exp(-dE/initial_temp);
    for(int i=0; i<N; i++) total_number_accepted[i]=0, Probability[i]=0;


    double number_mcs[N];
    number_mcs[0] = mc;
    for(int i=1; i<N; i++) number_mcs[i] = mc + i*mc;


    for(int i=0; i<5; i++) average[i] = 0;

    int**s_matrix;

    s_matrix = (int**) matrix(L, L, sizeof(L));

    Initialize(s_matrix, L, E, M, idum, random, other_spin);

    for(int k=0; k<N; k++){
        for(int l=0; l<mc; l++){
            accepted = 0;
            Metropolis(s_matrix, L, E, M, w, idum, accepted, other_spin);
            average[0] += E;
            average[1] += E*E;
            average[2] += M;
            average[3] += M*M;
            average[4] += fabs(M);
            total_number_accepted[k] += accepted;
        }
        cout << average[0]/number_mcs[k] << endl;
        Expectation_Values(average, number_mcs[k], total_number_accepted[k], Heat_Capacity[k], Susceptibility[k], Mean_M[k], L, Average_E[k], initial_temp);

        if(number_mcs[k] > 100000) Energy_Probability(s_matrix, L, Probability[k], other_spin);


    }
    free_matrix(((void **) s_matrix)); //freeing memory
    Save_Results_Matlab(Filename, Heat_Capacity, Susceptibility, Mean_M, N, number_mcs, total_number_accepted, Probability, Average_E);
    return 0;

}

//=========================
//The Metropolis algorithm
//=========================
void Metropolis(int **s_matrix, int spin, double &E, double &M, double *w, long &idum, int &accepted, int **other_spin){
    for(int r=0; r<spin; r++){
        for(int c=0; c<spin; c++){
            int mcr = (int) (ran2(&idum)*(double)spin);
            int mcc = (int) (ran2(&idum)*(double)spin);

            int dE = 2*s_matrix[mcr][mcc]*
                    (s_matrix[mcr][other_spin[mcc][0]] +
                    s_matrix[other_spin[mcr][0]][mcc] +
                    s_matrix[mcr][other_spin[mcc][1]] +
                    s_matrix[other_spin[mcr][1]][mcc]);
            if(dE <= 0 || ran2(&idum) <= w[dE+8]){
                accepted += 1;
                s_matrix[mcr][mcc] *= -1;
                E += (double) dE;
                M += (double) 2*s_matrix[mcr][mcc];
            }
        }
    }
}


//=======================
//Function to initialize
//=======================
void Initialize(int **s_matrix, double spin, double &E, double &M, long &idum, bool random, int **other_spin){
    if(random){
        for(int r=0; r<spin; r++){
            for(int c=0; c<spin; c++){
                if(ran2(&idum) <= 0.5) s_matrix[r][c] = 1.;
                else s_matrix[r][c] = -1;
            }
        }
    }

    else{
        for(int r=0; r<spin; r++){
            for(int c=0; c<spin; c++){
                s_matrix[r][c] = 1.;
            }
        }
    }
    for(int r=0; r<spin; r++){
        for(int c=0; c<spin; c++){
            E -= (double) s_matrix[r][c]*(s_matrix[other_spin[r][0]][c] +
                    s_matrix[r][other_spin[c][0]]);
            M += (double) s_matrix[r][c];
        }
    }
}

//==================================================================================
// Function to save my results in a matlab file, easy to further analyse the result
//==================================================================================
void Save_Results_Matlab(string Filename, double *Heat_Capacity, double *Susceptibility, double *Mean_M, int N, double *mcs, double *total_number_accepted, double *Probability, double *Average_E){

    ofstream myfile;
    myfile.open(Filename);
    myfile << "mc " << "= [";
    for (int i=0; i<N; i++){
        myfile << mcs[i] << ", ";
    }
    myfile << "];" << endl;

    myfile << "M_var " << "= [";
    for (int i=0; i<N; i++){
        myfile << Susceptibility[i] << ", ";
    }
    myfile << "];" << endl;

    myfile << "E_var " << "= [";
    for (int i=0; i<N; i++){
        myfile << Heat_Capacity[i] << ", ";
    }
    myfile << "];" << endl;

    myfile << "MeanM " << "= [";
    for (int i=0; i<N; i++){
        myfile << Mean_M[i] << ", ";
    }
    myfile << "];" << endl;

    myfile << "Probability " << "= [";
    for (int i=0; i<N; i++){
        myfile << Probability[i] << ", ";
    }
    myfile << "];" << endl;

    myfile << "E_avg " << "= [";
    for (int i=0; i<N; i++){
        myfile << Average_E[i] << ", ";
    }
    myfile << "];" << endl;

    myfile << "accepted_tot " << "= [";
    for (int i=0; i<N; i++){
        myfile << total_number_accepted[i] << ", ";
    }
    myfile << "];" << endl;

    myfile << "plot(mc, E_avg)" << endl;

    myfile.close();

}

//============================================
//Function to calculate the energy probability
//============================================
void Energy_Probability(int **s_matrix, int spin, double &E, int **other_spin){
    E = 0;
    for(int r=0; r<spin; r++){
        for(int c=0; c<spin; c++){
            E -= (double) s_matrix[r][c]*(s_matrix[other_spin[r][0]][c] +
                    s_matrix[r][other_spin[c][0]]);
        }
    }
}


//============================================
//Function to calculate the expectation values.
//=============================================
void Expectation_Values(double *average, double norm, double &total_number_accepted, double &Heat_Capacity, double &Susceptibility, double &Mean_M, int L, double &Average_E, double temp){
    total_number_accepted = total_number_accepted/norm;
    Heat_Capacity = (average[1]/norm - average[0]*average[0]/(norm*norm))/(L*L*temp*temp);
    Susceptibility = (average[3]/norm - average[4]*average[4]/(norm*norm))/(L*L*temp);
    Mean_M = average[4]/(norm*L*L);
    Average_E = average[0]/(norm*L*L);
}

#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include "lib.h"
#include "time.h"

using namespace std;

// Periodic boundary conditions.

inline int periodic(int i, int limit, int add) {
    return (i+limit+add) % (limit);
}


void Metropolis(int **s_matrix, int spin, double &E, double &M, double *w, long &idum, int &accepted, int **Spins);

void SaveResults(string FileName, double *Heat_Capacity, double *Susceptibility, double *Mean_M, double *MCS, int N, double *average_total, double *, double *Average_E);

void Initialize(int **, double, double & , double & , long &, bool random, int **);

void EnergyProbability(int **s_matrix, int spin, double &E, int **Spins);

void Expectation_Values(double *values, double normalization, double &accepted_total, double &Heat_Capacity, double &Susceptibility, double &Mean_M, int L, double &Average_E, double temp);

int main(){
    //================
    // Main parameters
    //================
    int MCS = 10000;
    int MC = 10;
    int N = MCS/MC;
    double initial_temp = 1;
    int L = 2; // Lattice
    bool random = true; // To start with random state
    string FileName = "this_file_MC.m";
    int temp = 1;
    double ran2(long *);

    double average[5]; // Results array
    double E=0, M=0;
    double w[17]; // Energies array
    long idum  = -1;

    //Spins
    int **Spins;
    Spins = (int **) matrix(L,2,sizeof(L));
    for(int i=0; i<L; i++){
        Spins[i][0] = periodic(i,L,-1);
        Spins[i][1] = periodic(i,L,1);
    }


    //====================
    // Monte Carlo cycles:
    //====================
    double number_MCS[N];
    number_MCS[0] = MC;
    for(int i=1; i<N; i++) number_MCS[i] = MC + i*MC; //cout << number_MCS[i] << endl;



    //    double number_MCS[N];
    //    number_MCS[N] = MC;


    //cout << N << endl;

    double Heat_Capacity[N];
    double Susceptibility[N];
    double Mean_M[N];
    double Probability[N];
    double Average_E[N];
    double average_total[N];
    int accepted;
    //double t[N];


    // Initialising w
    //for(int delta_E=-8;delta_E<=8;delta_E++) w[delta_E+8] = exp(-delta_E/start_temperature);
    //for(int i=0;i<N;i++) accepted_total[i]=0, probability[i]=0;


    //for(int dE =-8; dE <= 8; dE++) w[dE+8] = 0;
    for(int dE =-8; dE <= 8; dE++) w[dE+8] = exp(-dE/temp); //hva skjer her?
    for(int i=0; i<N; i++) average_total[i]=0, Probability[i]=0;

    for(int i = 0; i < 5; i++) average[i] = 0.;
    //initialize(L, temp, s_matrix, E, M);




    clock_t start, finish;
    int **s_matrix;


    s_matrix = (int **) matrix(L, L, sizeof(L));
    Initialize(s_matrix, L, E, M, idum, random, Spins);
    //cout << average[0] << endl;
    //cout << E << endl;
    int bolle =0;
    start = clock();
    for(int i=0; i<N; i++){
        for(int j=0; j<MC; j++){
            accepted = 0;
            Metropolis(s_matrix, L, E, M, w, idum, accepted, Spins);
            //cout<< E << endl;
            average[0] += E;
            average[1] += E*E;
            average[2] += M;
            average[3] += M*M;
            average[4] += fabs(M);
            average_total[i] += accepted;
            //Expectation_Values(average, number_MCS[j], average_total[j], Heat_Capacity[j], Susceptibility[j], Mean_M[j], L, Average_E[j], initial_temp);
            bolle +=1;
        }
        //cout << accepted << endl;
        //cout << average[0]/bolle << endl;
        //cout << E/number_MCS[i] << endl;
        Expectation_Values(average, number_MCS[i], average_total[i], Heat_Capacity[i], Susceptibility[i], Mean_M[i], L, Average_E[i], initial_temp);
        cout << average[0]/bolle << endl;
        if( number_MCS[i] > 100000 ) EnergyProbability(s_matrix, L, Probability[i], Spins);
        finish = clock();
        //t[i] = ((double) (finish-start)/CLOCKS_PER_SEC);

    }

    free_matrix(((void **) s_matrix)); // free memory
    SaveResults(FileName, Heat_Capacity,Susceptibility, Mean_M, number_MCS, N, average_total, Probability, Average_E);
    return 0;
}

//===========================
// The Metropolis algorithm
//===========================

void Metropolis(int **s_matrix, int spin, double &Energy, double &M, double *w, long &idum, int &accepted, int **Spins){
    for(int j=0; j<spin; j++){
        for (int i=0; i<spin; i++){
            int x = (int) (ran2(&idum)*(double)spin);
            int y = (int) (ran2(&idum)*(double)spin);
            int dE = 2*s_matrix[x][y]*(s_matrix[x][Spins[y][0]] +
                    s_matrix[Spins[x][0]][y] +
                    s_matrix[x][Spins[y][1]] +
                    s_matrix[Spins[x][1]][y]);
            if (dE <= 0 || ran2(&idum) <= w[dE+8] ){
                accepted += 1;
                s_matrix[y][x] *= -1;
                Energy += (double) dE;
                M += (double) 2*s_matrix[y][x];
            }
        }
    }
    cout << Energy << endl;
}

//=========================
// Initialzing the system
//=========================
void Initialize(int **s_matrix, double spin, double &E, double &M, long &idum, bool random, int **Spins){
    if(random){
        for(int r=0; r<spin; r++){
            for(int c=0; c<spin; c++){
                if(ran2(&idum) <=0.5) s_matrix[r][c] = 1.;
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
        for(int c=0;c<spin; c++){
            E -= (double) s_matrix[r][c]*(s_matrix[Spins[r][0]][c] +
                    s_matrix[r][Spins[c][0]]);
            M += (double) s_matrix[r][c];
        }
    }




    //    for(int i=0; i < spin; i++) {
    //        for (int j= 0; j < spin; j++){
    //            s_matrix[i][j] = 1; // spin orientation for the ground state
    //            M +=  (double) s_matrix[i][j];
    //        }
    //    }
    //    // setup initial energy
    //    for(int i=0; i < spin; i++) {
    //        for (int j=0; j < spin; j++){
    //            E -=  (double) s_matrix[i][j]*
    //                    (s_matrix[periodic(i,spin,-1)][j] +
    //                    s_matrix[i][periodic(j,spin,-1)]);
    //        }

    //    }
    //    //cout << E << endl;

}// end function initialise


//==============================================
// Function to save my results in a matlab-file
//==============================================
void SaveResults(string FileName, double *Heat_Capacity, double *Susceptibility, double *Mean_M, double *MCS, int N, double *average_total, double *Probability, double *Average_E){

    ofstream myfile;
    myfile.open(FileName);
    myfile << "MCS " << "= [";
    for (int i=0; i<N; i++){
        myfile << MCS[i] << ", ";
    }
    myfile << "];" << endl;

    myfile << "E_var " << "= [";
    for (int i=0; i<N; i++){
        myfile << Heat_Capacity[i] << ", ";
    }
    myfile << "];" << endl;

    myfile << "M_var " << "= [";
    for (int i=0; i<N; i++){
        myfile << Susceptibility[i] << ", ";
    }
    myfile << "];" << endl;

    myfile << "mean_M " << "= [";
    for (int i=0; i<N; i++){
        myfile << Mean_M[i] << ", ";
    }
    myfile << "];" << endl;

    //    myfile << "time_used " << "= [";
    //    for (int i=0; i<N; i++){
    //        myfile << t[i] << ", ";
    //    }
    //    myfile << "];" << endl;
    myfile << "accepted" << "= [";
    for (int i=0; i<N; i++){
        myfile << average_total[i] << ", ";
    }
    myfile << "];" << endl;

    myfile << "Probability " << "= [";
    for (int i=0; i<N; i++){
        myfile << Probability[i] << ", ";
    }
    myfile << "];" << endl;

    myfile << "av_E " << "= [";
    for (int i=0; i<N; i++){
        myfile << Average_E[i] << ", ";
    }
    myfile << "];" << endl;

    myfile << "plot(MCS, E_var)" << endl;

    myfile.close();
}


//============================================
// Function to compute the energy probability
//============================================
void EnergyProbability(int **s_matrix, int spin, double &E, int **Spins){
    E = 0;
    for(int i=0; i<spin; i++){
        for(int j=0; j<spin; j++){
            E -= (double) s_matrix[i][j]*(s_matrix[Spins[i][0]][j] +
                    s_matrix[i][Spins[j][0]]);
        }
    }
}



//=============================================
// Function to compute the expectation values
//=============================================
void Expectation_Values(double *average, double norm, double &accepted_total, double &Heat_Capacity, double &Susceptibility, double &Mean_M, int L, double &Average_E, double temp){

    accepted_total = accepted_total/norm;
    Heat_Capacity = (average[1]/norm - average[0]*average[0]/(norm*norm))/(L*L*temp*temp);
    Susceptibility = (average[3]/norm - average[4]*average[4]/(norm*norm))/(L*L*temp);

    Mean_M = average[4]/(norm*L*L);
    Average_E = average[0]/(norm*L*L);
    //cout << Susceptibility << endl;


    //    double E2average = average[1]*norma;
    //    double Maverage = average[2]*norma;
    //    double M2average = average[3]*norma;
    //    double Mabsaverage = average[4]*norma;
    //    // all expectation values are per spin, divide by 1/n_spins/n_spins
    //    double Evariance = (E2average- Eaverage*Eaverage)/L/L;
    //    double Mvariance = (M2average - Maverage*Maverage)/L/L;
    //    double M2variance = (M2average - Mabsaverage*Mabsaverage)/L/L;
    //    double Mvariance = (M2average - Mabsaverage*Mabsaverage)/n_spins/n_spins;
}

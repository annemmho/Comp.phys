#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include "lib.h"
#include "time.h"
#include "mpi.h"


using namespace std;

ofstream ofile;

inline int periodic(int i, int limit, int add){
    return (i+limit+add) % (limit);
}

void Metropolis(int **s_matrix, int spin, double &E, double &M, double *w, long &idum, int &accepted, int **other_spin);

void Initialize(int **, double, double &, double &, long &, bool, int **);

void Expectation_Values(double *average, double norm, double &total_number_accepted, double &Heat_Capacity,
                        double &Susceptibility, double &Mean_M, int L, double &Average_E);

void Energy_Probability(int **s_matrix, int spin, double &E, int **other_spin);

void Saving(string Filename, double *E_var, double *M_var, double *Mean_M, double *temp, int steps, double *Mean_E);


int main(int argc, char* argv[])
{
    // The parameters
    int mcs = 10000000;
    double initial_temp = 2.1, final_temp = 2.6, temp_step = 0.05;
    int N = (final_temp - initial_temp)/temp_step + 1;
    int L = 40;
    bool random = true;
    string Filename = "met_is_mc_para.m";


    int  my_rank, numprocs, **s_matrix, **other_spin;



    //    double number_mcs[N];
    //    number_mcs[0] = mc;
    //    for(int i=1; i<N; i++) number_mcs[i] = mc + i*mc;

    //MPI initialization
    MPI_Init (&argc, &argv);
    MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
    //    if (my_rank == 0 && argc <= 1) {
    //      cout << "Bad Usage: " << argv[0] <<
    //        " read output file" << endl;
    //      exit(1);
    //    }


    int no_int = N/numprocs;
    int myloop_start = my_rank*no_int - 1;
    int myloop_end = (my_rank+1)*no_int;
    if((my_rank == numprocs-1) &&(myloop_end < N) ) myloop_end = N-1;

    MPI_Bcast (&L, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast (&mcs, 1, MPI_INT, 0, MPI_COMM_WORLD);
    //MPI_Bcast (&final_temp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    //MPI_Bcast (&temp_step, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    double average[5];
    double E = 0;
    double M = 0;
    double w[17];


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
    double Mean_E[N];
    double temp[N];
    double accepted_tot[N];
    for(int i=0;i<N;i++) temp[i] = initial_temp + i*temp_step;


    double total_HC[N];
    double total_S[N];
    double total_t[N];
    double total_ME[N];
    double total_MM[N];

    for(int i=0; i<N; i++) total_HC[i]=0, total_S[i]=0, total_t[i]=0,
            total_ME[i]=0, total_MM[i]=0;

    for(int i=0; i<N; i++) Heat_Capacity[i]=0, Susceptibility[i]=0,
            Mean_M[i]=0, Mean_E[i]=0, accepted_tot[i]=0;

    long idum = -1- my_rank;

    s_matrix = (int **) matrix(L,L, sizeof(L));

    Initialize(s_matrix, L, E=0, M=0, idum, random, other_spin);

    int a = 0;

    for(int cycle = myloop_start; cycle <= myloop_end; cycle++){
        a = cycle;

        for(int i=0; i<5; i++) average[i] = 0;

        for(int dE =-8; dE <= 8; dE++) w[dE+8] = 0;
        for(int dE = -8; dE <= 8; dE+=4) w[dE+8] = exp(-dE/temp[a]);
        for(int i=0; i<mcs; i++){
            int accepted = 0;

            Metropolis(s_matrix, L, E, M, w, idum, accepted, other_spin);

            average[0] += E;
            average[1] += E*E;
            average[2] += M;
            average[3] += M*M;
            average[4] += fabs(M);
            total_number_accepted[cycle] += accepted;
        }
        Expectation_Values(average, mcs, total_number_accepted[a], Heat_Capacity[a], Susceptibility[a], Mean_M[a], L, Mean_E[a]);
        if(my_rank == 0){
            cout << "Something: " << temp << endl;
        }

    }







    free_matrix(((void **) s_matrix)); //freeing memory

    MPI_Reduce(&Heat_Capacity, &total_HC, N, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&Susceptibility, &total_S, N, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&Mean_M, &total_MM, N, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&temp, &total_t, N, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&Mean_E, &total_ME, N, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);


    MPI_Finalize ();

    if ( my_rank == 0) {
        for(int i=0; i<N; i++) cout << total_HC[i] << "  ";
        Saving(Filename, total_HC, total_S, total_MM, temp, N, total_ME);
    }

    return 0;

}

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


void Saving(string Filename, double *E_var, double *M_var, double *Mean_M, double *temp, int steps, double *Mean_E){
    ofstream myfile;
    myfile.open(Filename);
    myfile << "E_var " << "= [";
    for (int i=0; i<steps; i++){
        myfile << E_var[i] << ", ";
    }
    myfile << "];" << endl;

    myfile << "Mean_E " << "= [";
    for (int i=0; i<steps; i++){
        myfile << Mean_E[i] << ", ";
    }
    myfile << "];" << endl;

    myfile << "M_var " << "= [";
    for (int i=0; i<steps; i++){
        myfile << M_var[i] << ", ";
    }
    myfile << "];" << endl;

    myfile << "Mean_M " << "= [";
    for (int i=0; i<steps; i++){
        myfile << Mean_M[i] << ", ";
    }
    myfile << "];" << endl;

    myfile << "Temperature " << "= [";
    for (int i=0; i<steps; i++){
        myfile << temp[i] << ", ";
    }
    myfile << "];" << endl;


    myfile.close();

}

void Energy_Probability(int **s_matrix, int spin, double &E, int **other_spin){
    E = 0;
    for(int r=0; r<spin; r++){
        for(int c=0; c<spin; c++){
            E -= (double) s_matrix[r][c]*(s_matrix[other_spin[r][0]][c] +
                    s_matrix[r][other_spin[c][0]]);
        }
    }
}

void Expectation_Values(double *average, double norm, double &total_number_accepted, double &Heat_Capacity, double &Susceptibility, double &Mean_M, int L, double &Mean_E){
    total_number_accepted = total_number_accepted/norm;
    Heat_Capacity = (average[1]/norm - average[0]*average[0]/(norm*norm))/(L*L);
    Susceptibility = (average[3]/norm - average[4]*average[4]/(norm*norm))/(L*L);
    Mean_M = average[4]/(norm*L*L);
    Mean_E = average[0]/(norm*L*L);
}












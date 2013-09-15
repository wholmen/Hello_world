#include <iostream>
#include <armadillo>
#include <cmath>
#include <fstream>
#include <stdio.h>
#include <time.h>

using namespace std;
using namespace arma;

inline float analytical(float x);
inline float logerror(float v, float u);

int main()
{
    // This program will solve a general tri diagonal linear matrix equation Av = b.

    int N = 10;
    double x = 1; double x0 = 0; double h = (x-x0)/(N+1);
    double a1[N]; double a2[N]; double a3[N]; double v[N]; double b[N];

    // Filling in vectors
    int j;
    for (j=0; j<N; j++){
        a1[j] = -1.0; a3[j] = -1.0; a2[j] = 2.0; b[j] = 100*exp(-10*h*j)*pow(h,2); v[j] = 0.0;
    }

    // Recording time
    clock_t start, finish;
    start = clock();

    // Forward Substitution. Row reducing the matrix equation
    int i;
    for (i = 1; i <= N-1; i++){
        float factor = a1[i]/a2[i-1];
        a2[i] = a2[i] - a3[i-1]*factor;
        b[i] = b[i] - b[i-1]*factor;
    }

    // Backward Substitution. Solving the equation for vector v.
    int k;
    v[N-1] = b[N-1]/a2[N-1];
    for (k = N-1; k >= 0; k--){
        v[k] = (b[k] - a3[k] * v[k+1])/a2[k];
    }

    // Computing the analytical solution, u. Will be compared with numerical solution v for variying step length, h.
    double u[N]; double epsilon[N]; //Declaring u and errorfunction epsilon.
    for (i=0;i<N-1;i++){
        u[i] = analytical(h*i);
        epsilon[i] = logerror(v[i],u[i]);  // Logaritmic numerical error.
    }

    finish = clock(); // Stopping time recording
    float algotime = ( (finish - start)/CLOCKS_PER_SEC );


    // The following code will implement armadillos LU-decomposition and compare computational speed with my algorithm

    // Implementing LU-decomposition to solve the equation.
    mat A = zeros<mat>(N,N);
    // Filling in elements in A
    A(N-1,N-1) = 2;
    for (i=0;i<N-1;i++){
        A(i,i) = 2;
        A(i,i+1) = -1; A(i+1,i) = -1;
    }

    // Recording time
    clock_t startLU, finishLU;
    startLU = clock();


    // Algorithm for LU-decomposition using armadillo library.
    mat L,U;
    lu(L, U, A);

    vec z[N]; vec v2[N]; //initializing vector v2, solution vector and z = U*v2

    for (i=0;i<N;i++){ // Computing Lz = b
        z[i] = b[i];
        for (j=0;j<i;j++){
            z[i] = z[i] - z[j];
        }
    }

    for (i=N-1;i>=0;i--){
        v2[i] = z[i];
        for (j=N-1;j>i;j--)
            v2[i] = v2[i] - U(i,j)*v[j];
        v2[i] = v2[i]/U(i,i);
    }


    finishLU = clock(); // Stopping time recording
    float algotimeLU = ( (finishLU - startLU)/CLOCKS_PER_SEC );

    // Writing the output to file. PS: Change name of file if N changes.
    ofstream myfile;
    myfile.open("b_n10.txt");

    myfile << "Sorted order: x, numerical v, analytical u, error epsilon, time own algo, time LU algo, LU solution." << endl;
    for (i=0; i<N; i++){
        myfile << i*h <<" "<< v[i] <<" "<< u[i] <<" "<< epsilon[i] <<" "<< algotime <<" "<< algotimeLU <<" "<< v2[i];
    }
    myfile.close();

}

inline float analytical(float x){
    return 1 - (1-exp(-10))*x-exp(-10*x);
}

inline float logerror(float v, float u){
    return log(abs((v-u)/u));
}

#include <iostream>
#include <armadillo>
#include <cmath>
#include <fstream>

using namespace std;
using namespace arma;

float u(float x){
    return 1 - (1-exp(-10))*x-exp(-10*x);
}

int main()
{
    int N = 10; float n = 10;
    double x = 1; double x0 = 0;
    double h = (x-x0)/(n+1);

    // This Program will be divided into 3 parts. One part will solve generally for a matrix
    // One part will solve as a tri diagonal matrix, and one will be specific for this situation.

    // Part 2.
    double a1[N]; double a2[N]; double a3[N]; double v[N]; double b[N];

    // Filling in vectors
    int j;
    for (j=0; j<N; j++){
        a1[j] = -1.0; a3[j] = -1.0; a2[j] = 2.0; b[j] = 100*exp(-10*h*j); v[j] = 0.0;
    }

    // Forward Substitution
    int i;
    for (i = 1; i <= N-1; i++){
        float factor = a1[i]/a2[i-1];
        a2[i] = a2[i] - a3[i-1]*factor;
        b[i] = b[i] - b[i-1]*factor;
    }

    // Backward Substitution
    int k;
    v[N-1] = b[N-1]/a2[N-1];
    for (k = N-1; k >= 0; k--){
        v[k] = (b[k] - a3[k] * v[k+1])/a2[k];
    }

    //float analytical[N];
    //for (i=0; i<N-1; i++){
    //    analytical[i] = u(i*h);
    //}

    ofstream myfile;
    myfile.open("b_n10.txt");

    for (i=0; i<N; i++){
        cout << v[i] << endl;
    }

    for (i=0; i<N; i++){
        myfile << i*h << " " << v[i] << endl;
    }
    myfile.close();

}

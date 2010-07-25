#include<iostream>
using namespace std;

void allocateMem( double **&u, int a, int b ) {
    u = new double*[a];
    for( int i = 0; i < a; i++)
    u[i]= new double[b];
}

void freeMem( double **&u, int a ) {
    for( int i = 0; i < a; i++)
    delete [] u[i];
    delete [] u;
}
   
int main() {
    int a=2024;
    double **u;
    allocateMem( u, a, a );
    u[a-1][a-1] = 20;
    cout<<u[a-1][a-1]<<endl;
    freeMem( u, a );
}

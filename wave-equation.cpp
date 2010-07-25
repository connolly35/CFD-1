#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include "gnuplot_i.hpp"

using namespace std;

void ftcs( double *u, int n, double cfl ) {
    double *v = new double[n];
    for( int i = 0; i < n ; i++ ) 
        v[i] = u[i];
    for( int i = 1; i < n - 1 ; i++ ) 
        u[i] = v[i] - ( cfl / 2 ) * ( v[i+1] - v[i-1] );
    u[n-1] = u[n-2];
    delete v;
    
}

void ftfs( double *u, int n, double cfl ) {
    double *v = new double[n];
    for( int i = 0; i < n ; i++ ) 
        v[i] = u[i];
    for( int i = 1; i < n-1 ; i++ ) 
        u[i] = v[i] - cfl * ( v[i+1] - v[i] );
    u[n-1] = u[n-2];
    delete v;
    
    
}

void ftbs( double *u, int n, double cfl ) {
    double *v = new double[n];
    for( int i = 0; i < n ; i++ ) 
        v[i] = u[i];
    for( int i = 1; i < n; i++ ) 
        u[i] = v[i] - cfl * ( v[i] - v[i-1] );
    delete v;
    
}

int main() {
    int n, c, s, ts = 1;
    double cfl;
    char a='n';
    double *u;
	Gnuplot g1("lines");
    g1.cmd("set term wxt persist");
    g1.cmd("set zero 1e-300");

    cout<<"Enter the scheme to be used \n\n1 for FTCS\n2 for FTFS\n3 for FTBS\nScheme (default FTBS) : ";
    cin>>s;
    cout<<endl; 
    cout<<"Enter the number of grid points to be used\n";
    cin>>n;
    cout<<endl; 
    cout<<"Enter the CFL number to be used\n";
    cin>>cfl;
    
    /* Initialize the BC 
        */
    u = new double[n];
    u[0] = 1;
    for( int i = 1; i < n; i++ ) u[i] = 0;

    vector <double> x(n);
    vector <double> U(n);
    for( int i = 0; i < n; i++ ) {
        x[i] = double(i)/double(n);
    }

    while (a!='e') {
        
        switch ( a ) {
            case 'm': {
                cout<<"\n\n\n";
                cout<<" Runtime Menu:\n 1. Set time step ( current value is "<<ts<<" ) \n 2. Change CFL number ( current value is "<<cfl<<" )\n 3. Change the scheme\n 4. Return\n ";
                cin>>c;
                switch(c) {
                    cout<<endl;
                    case 1: cout<<"Enter the time step : "; cin>>ts; break;
                    case 2: cout<<"Enter the new cfl number : "; cin>>cfl; break;
                    case 3: cout<<"Enter the scheme to be used \n\n1 for FTCS\n2 for FTFS\n3 for FTBS\nScheme (default FTBS) : "; cin>>s; break;
                    case 4: break;
                    default: cout<<"Invalid choice! \n ";
                }
                break; 
            }

            case 'n': {
                for( int i = 0; i< ts; i++)
                    switch(s) {
                        case 1: ftcs( u, n, cfl ); break;
                        case 2: ftfs( u, n, cfl ); break;
                        case 3:
                        default:
                            ftbs( u, n, cfl ); break;
                    }
                    
                for( int i = 0; i < n; i++ ) {
                    U[i] = u[i];
                }
   	       	    g1.set_style("lines").plot_xy(x, U, "wave", " pt 40");
                g1.reset_plot();
                break;
            }

            case 'e': {
                delete u;
                break;
            }
        }
        if( a != 'e' ) {
            cout<<"\n\n\n";
            cout<<"Press \n n for next iteration \n m for menu \n e for exit \n ";
            cin>>a;
        }
    }
} 

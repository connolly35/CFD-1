/* Program to solve the Euler equation for a 1-D constant area nozzle, with pressure boundary conditions specified.
   - M.Prakash
     AE07B014
    */
#include<iostream>
#include<math.h>
#include <vector>
#include <string>
#include <fstream>
#include "gnuplot_i.hpp"

using namespace std;

//Defining global constants
const double G = 1.4, Patm = 1e5, R = 287, T = 300;

//Function to allocate memory to a 2D array
void malloc2d( double **&array, unsigned int nrow, unsigned int ncol ) {
    array = new double*[nrow];
    for( int i = 0; i < nrow; i++)
        array[i]= new double[ncol];
}

//Function to release the memory allocated to a 2D array
void freeMem( double **&array, unsigned int nrow ) {
    for( int i = 0; i < nrow; i++)
        delete [] array[i];
    delete [] array;
}

//Function to obtain rho u P and T from Q
void q2vars(double q[3],double array[], double gam = 1.4) {
    array[0] = q[0];
    array[1] = q[1]/q[0];
    array[2] = (gam - 1) * (q[2] - (array[0]*array[1]*array[1]) *0.5);
    array[3] = array[2]/(array[0] * 287.026);
    return;
}

//Function to obtain the q matrix from primitive vars
void vars2q(double rho, double u, double P, double array[], double gam = 1.4) {
    array[0] = rho;
    array[1] = rho * u;
    array[2] = P/(gam-1) + (rho * u * u * 0.5);
    return;
}

//Function to get E given Q
void EofQ( double q[3], double e[3] ) {
    double p;
    e[0] = q[1];
    p = ( q[2] - 0.5 * pow( q[1], 2 ) / q[0] ) * (G-1);
    e[1] = p + pow(q[1],2)/q[0];
    e[2] = (q[2] + p)*q[1]/q[0]; 
}

main() {

    double **Q,**E,**q; 
    double delT, delX = 0.1, dTbydX;
    int ts,N, choice;
    double var1temp[4], var2temp[4], M;
    Gnuplot g1("lines");
    Gnuplot g2("lines");
    Gnuplot g3("lines");

    g1.cmd("set term wxt persist");
    g2.cmd("set term wxt persist");
    g3.cmd("set term wxt persist");
    //User Inputs

    cout<<"\nEnter the number of grid points: ";
    cin>>N;
    delX = 1.0/double(N);
    cout<<"\nEnter CFL(working range is around 1e-4 for 10, 1e-6 for 100, please try lower CFL if code fails):  ";
    cin >> dTbydX;
    cout<<"\nEnter Number of time steps:   ";
    cin>>ts;

    //Initializing Q,E

    malloc2d( Q, N, 3);
    malloc2d( q, N, 3);
    malloc2d( E, N, 3);
    vector <double> rho(N),pressure(N),temp(N),vel(N),x(N);

    Q[0][0]= q[0][0]= 1.2*Patm/(R*T);
    Q[0][1]= q[0][1]= 0;
    Q[0][2]= q[0][2]= 1.2*Patm/(G -1);
    E[0][0]= 0;
    E[0][1]= 1.2*Patm;
    E[0][2]= 0;

    for (int i=1; i<N; i++) {
        Q[i][0] = q[i][0] = Patm/(R*T);
        Q[i][1] = q[i][1] = 0;
        Q[i][2] = q[i][2] = Patm/(G -1);
        E[i][0] = 0;
        E[i][1] = Patm;
        E[i][2] = 0;
    }

    // Forward march in time 
iterate:

    while ( ts ) {
        
        for (int j=1; j<N-1; j++) {
            for (int k=0; k<3; k++) {
                q[j][k] = Q[j][k] - (dTbydX/2)*(E[j+1][k] - E[j-1][k]) + (Q[j+1][k] + Q[j-1][k] - 2*Q[j][k])*0.1*dTbydX/delX;
                if( j > 1 && j < N-2 )
                q[j][k] -= ( ( Q[j+2][k] + Q[j-2][k] ) + 6*Q[j][k] - 4*( Q[j+1][k] + Q[j-1][k] ) ) * 0.01 * dTbydX/pow(delX,3);
            }

        }

        //Transfer velocity of the second point to the initial point
		q2vars(q[0], var1temp, 1.4);                      //obtain primitive vars at Q[0]
		var1temp[1] = q[1][1] / q[1][0];                                    //u(0) = u(1)
		var1temp[3] = 300 - (var1temp[1] * var1temp[1] * 0.5)/1004.6; //change in T due to nonzero u
		M = var1temp[1]/sqrt(1.4*287.026*var1temp[3]);                 //define a Mach number for faster reference
		var1temp[2] = 1.2*Patm/pow((1 + 0.2 * M * M),3.5);                    //new value of P
		var1temp[0] = var1temp[2] / ( 287.026 * var1temp[3] );          //new rho
		vars2q(var1temp[0], var1temp[1], var1temp[2], q[0], 1.4);   //convert back to q variables

    
        //Iteration for the last point

		//set of transformations at the end point, since we transfer u and T to end point
		q2vars(q[N-1], var1temp, 1.4);                //obtain prim vars at Q[N-1]
        q2vars(q[N-2], var2temp, 1.4);                //obtain prim vars at Q[N-2]
        var1temp[1] = var2temp[1];                                       //equating velocities, i.e u[N-1] = u[N-2]
        var1temp[3] = var2temp[3];                                       //eqauting Temperatures, i.e. T[N-1] = T[N-2]
        M = (var1temp[1]/sqrt(1.4*287.026*var1temp[3]));                 //define Mach number for faster reference
        var1temp[0] = var1temp[2] / (287.026 * var1temp[3]);          //new rho
        vars2q(var1temp[0], var1temp[1], var1temp[2], q[N-1], 1.4);   //convert back to q variables

        //Update the Q and E matrices
        for (int i=0; i<N; i++) {
            for (int j=0; j<3; j++) 
                Q[i][j]=q[i][j];
            EofQ( Q[i], E[i] ); //Update E

        }
          
        ts--;
    }

    for (int i=0; i<N; i++) {
        q2vars(Q[i],var1temp);
        x[i] = double(i)/double(N);
        rho[i] = var1temp[0];
        vel[i] = var1temp[1];
        pressure[i] = var1temp[2];
    }
    
    //Plot the results
    g1.set_style("lines").plot_xy(x, vel, "Velocity", " lt 2 pt 40");
    g1.reset_plot();
    g2.set_style("lines").plot_xy(x, rho, "Density", " lt 3 pt 40");
    g2.reset_plot();
    g3.set_style("lines").plot_xy(x, pressure, "Pressure", " lt 4 pt 40");
    g3.reset_plot();
    
    //End user menu
    cout<<"\n\nRuntime Menu:\n"<<"\t1. Iterate for more timesteps (default)\n"<<"\t2. Exit\n\n"<<"Enter Choice: ";
    cin>>choice;
    if(choice != 2) {
        cout<<"\n\nEnter the number of time steps :";
        cin>>ts;
        goto iterate;
    }
       
    //Cleaning up
    freeMem(Q, N);
    freeMem(E, N);
    freeMem(q, N);
}


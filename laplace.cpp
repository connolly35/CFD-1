/* Program to solve 2D laplace equation with the bcs read from a file and the other parameters taken in as user input.

   - M.Prakash, AE07B014
*/
#include<iostream>
#include<fstream>
#include<stdlib.h>
#include<stdio.h>
#include<math.h>
using namespace std;

//Function to allocate memory to a 2D array
void allocateMem( double **&array, unsigned int nrow, unsigned int ncol ) {
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

//Function to initialize the u matrix
void init( double **&u, int *nn ) {
    int n;
    fstream f;
    f.open( "bc.dat", ios::in );
    f>>n;
    *nn=n;
    allocateMem( u, n, n );
    for( int i = 0; i < n; i++ )
        f>>u[i][0];
    for( int j = 1; j < n; j++ )
        f>>u[n-1][j];
    for( int i = n-2; i >= 0; i-- )
        f>>u[i][n-1];
    for( int j = n - 2; j > 0; j-- )
        f>>u[0][j];

    for(int i = 1; i < n-1; i++ )
        for(int j = 1; j < n-1; j++ )
            u[i][j]=0;
}

//Function for 2D interpolation
double value( double **u, double x, double y, int nx, int ny ){

    double delx = 1.0/(nx-1.0);
    double dely = 1.0/(ny-1.0);
    int i = int(x/delx);
    int j = int(y/dely);
    
    if( i * delx < x )
        i++;
    if ( ( i-1 ) *delx > x )
        i--;
    if( j* dely < y )
        j++;
    if ( (j-1)*dely > y )
        j--;
    if( i > 0 )
        i--;
    if( j > 0 )
        j--;
    if( x == 0 )
        i = 0;
    if( x == 1 )
        i = nx-1;
    if( y == 0 )
        j = 0;
    if( y == 1 )
        j = ny-1;

    if( ( x == 1 && y == 1 )  || ( x == 0 && y == 0 )  || ( x == 1 && y == 0 ) || ( x == 0 && y == 1 ) )
        return u[i][j];

    if( x == 1 || x == 0 )
        return( ( u[i][j] * ( (j+1)*dely - y ) +  u[i][j+1] *  ( y - j*dely )  ) / dely );
    
    if( y == 1 || y == 0 )
        return( ( u[i][j] * ( (i+1)*delx - x ) +  u[i+1][j] *  ( x - i*delx )  ) / delx );

        return( ( u[i][j] * ( (i+1)*delx - x ) * ( (j+1)*dely - y ) + u[i+1][j] * ( x - i*delx ) * ( (j+1)*dely - y) + u[i][j+1] * ( (i+1)*delx - x ) * ( y - j*dely ) + u[i+1][j+1] * ( x - i*delx ) * ( y - j*dely ) ) / ( delx * dely ) );

}

int main() {
    int i, j, n, iter;
    double **u, w, e, norm, temp, prevNorm, x0,y0;
    int t = 15, ocount = 0, dcount = 0, level = 0;
    
    // For checking if the code is entering into oscillations or divergence
    bool trigger = false;
    
    // Inputs
    cout<<"Enter the value of w for SOR: ";
    cin>>w;
    cout<<"Enter the value of error tolerance: ";
    cin>>e;
    cout<<"Enter coordinates for which the output has to be given \n ( x followed by y, both in the range 0 to 1 ) : ";
    cin>>x0>>y0;

    // Variable init
    norm = 2*e;
    iter = 0;
    init( u, &n );

    //Iterate
    while( norm > e ) {
        prevNorm = norm;
        norm = 0;
        iter++;
        if( dcount > t ) {
            cout<<"The code seems to be diverging! Wrong boundary conditions?";
            break;
        }
        if( ocount > t ) {
            cout<<"\n \n Code went into oscillation and has been interrupted \n This is due to a higher value of precision requested than the noise that arises due to subtracting two very close numbers \n The error limits are variable, according to the given combination of boundary conditions, grid size, w and error tolerance. \n The iterations have been stopped at this stage and the accuracy has been reduced to a norm value of: "<<prevNorm<<"\n";
            break;
        }
        for( i = 1; i < n  - 1; i++ ) {
            for( j = 1; j < n - 1 ; j++ ) {
                temp = u[i][j];
                u[i][j] = w * ( ( u[i+1][j] + u[i-1][j] + u[i][j+1] + u[i][j-1] ) * 0.25 ) + ( 1.0 - w ) * u[i][j];
                norm += pow( ( u[i][j] - temp ), 2 );
            }
        }
        norm = pow( norm, 0.5 );
        if(trigger) {
            if( norm < prevNorm ) {
                trigger = false;
                ocount++;
                dcount = 0;
                level = 0;
            } else {
                level++;
                dcount++;
                if( level > 1 ) {
                    trigger = false;
                    ocount = 0;
                    level = 0;
                }
            }
        } else {
            if( prevNorm < norm ) {
                ocount++;
                dcount++;
                trigger = true;
                level = 0;
            } else {
                level++;
                if( level > 1 ) {
                    ocount = 0;
                    dcount = 0;
                    level=0;
                }
            }
        }
    }
    cout<<"\n\n\nThe final value of u(x0,y0), with linear interpolation between the grid points is: "<< value(u,x0,y0,n,n)<<endl;

// Free the allocated memory
    freeMem( u, n );
}


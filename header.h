// header for toyModel.cpp methods

#include <iostream>
#include <fstream>
#include <tuple>
#include <math.h>
#include <string.h>
#define _USE_MATH_DEFINES            // for access to M_PI constant

using namespace std;

//   toyModel.cpp methods (SEE DESCRIPTIONS IN toyModel.cpp)

void genData( int, double*, double, 
              double, double );

tuple< double, double > toyEpsilonInit( double, double, 
                                        double, int ); 
                                              
double toyEpsilon( double, double, 
                   double, int );

double toyPrior( double, double, double, 
                 double, double, double );

double toyLikelihood( double*, double, double, 
                      double, int );

double toyPosterior( double*, double, double, 
                     double, double, double, 
                     double, double, int );

double toyIntProb( double*, double, double, double, 
                   double, double, double, double, 
                   double, double, int, int, int, string fn = string("") );
                   
void toyPriorGrad( double*, double, double, double,
                   double, double, double );
                   
void toyLikelihoodGrad( double*, double*, double,
                        double, double, int, int );                    
                   
void sgldUpdate( double*, double*, double*, double*,
                 double, double, double, double, double,  
                 double, double, double, int, int, int );
                 
void mhldUpdate( double*, double*, double*, double*,
                 double, double, double, double, double,  
                 double, int*, int*, int, int);
                 
void randData( double*, double *, int, int );

double estStat( double*, double, double, double,
                double, double, double, double,
                double, int ); 
              
double estStatAdj( double*, double, double, double,
                   double, double, double, double,
                   double, double, double, double, int ); 
                 
// auto-correlation related functions

double kuboTime(    double, double * );      // calculates the autocorrelation time
                                                // used in the Kubo formula
void   aCovEst( double *, double *, double,
                double, double, double, int, int ); // return estimated autocovariance, C(t),
                                                // for t = 1,2,...,n (specify n)
                   
//   Random number procedures, in file rng.c , written in C

extern "C" {void   sgenrand (unsigned long seed);} // Set the seed
extern "C" {double genrand ();                   } // Make a standard uniform
extern "C" {double InvNormDist (double p);       } // Inverse cumulative normal
                                                   // Z = InvNormDist( genrand() )
                                                   // produces a standard normal.

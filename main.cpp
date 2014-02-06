// This code does... some stuff to test the bias of the SGLD sampling method

#include "header.h"

int main(){

  // Model: SGLD with multimodal tied means gaussian

  // Parameters
  
  // calculate true statistic (true)
  bool bolTrueStat = true;
  
  // data vector
  int N = 1e+2;      // size of data
  double * X;        // data vector
  X = new double[N]; // set size of data vector
  
  // gradients
  double * priGrad;
  priGrad = new double[2];
  
  double * likhGrad;
  likhGrad = new double[2];
  
  // size of batch
  int * n;
  n = new int[6];
  
  // sigmas - prior and likelihood
  double sig1 = sqrt( 10. ); // sigma on theta1 prior
  double sig2 = 1.;  // sigma on theta2 prior
  double sigX = sqrt( 2. );  // sigma on likelihood mixture gaussians
  
  // theta - prior means
  double theta1pri = 0.; // mean of theta1 prior
  double theta2pri = 0.; // mean of theta2 prior
  
  // thetas - data generator
  double theta1dat = 0.; // theta1 in data generation process
  double theta2dat = 1.; // theta2 in data generation process
  
  // annealing schedule params
  double * epsMin;  // minimum step size
  epsMin = new double[4];
  double * epsMax;  // maximum step size
  epsMax = new double[4];
  double gam    = 0.55;  // neg. exponent on polynomial annealing schedule
  int T         = 1e+6;  // number of steps
  
  // MH-ld step size
  double * eps;
  eps = new double[4];
  
  // region 1, 3 and 2
  double l11 = -0.5;//0.9;//0.;//  //lower bound on theta 1
  double u11 = 0.25;//1.65;//1.;//  // upper bound on theta 1
  double l21 = 0.8;//-1.9;//-0.75;//   // lower bound on theta 2
  double u21 = 1.8;//-0.9;//0.65;//   // upper bound on theta 2
  
  // reference region
  double l12 = -10.0;  // lower bound on theta 1
  double u12 = 10.0;   // upper bound on theta 1
  double l22 = -10.0;  // lower bound on theta 2
  double u22 = 10.0;   //  upper bound on theta 2
  
  // number of grid evaluation points along theta 1 and 
  // theta 2 dimensions (for calculating test statistics)
  int n1     =  1e+3;
  int n2     =  1e+3;
  
  // Open results file
  ofstream results;
  results.open( "output/results_seed19_reg1_n100.txt" );

  // Initialize random number generator
  static unsigned long seed = 19;
  sgenrand( seed );
  
  // Generate data from bimodal distribution
  genData( N, X, theta1dat, theta2dat, sigX );
  cout << N << " data points generated " << endl;
  
  
  // evaluate test stat ( modal ratio -- f(mode1)/f(mode2) )
  double testStat;
  testStat = toyPosterior( X, 0., 1., theta1pri, theta2pri, sig1, sig2, sigX, N ) /
             toyPosterior( X, 1., -1., theta1pri, theta2pri, sig1, sig2, sigX, N );
             
  cout << "test statistic (modal ratio) = " << testStat << endl;
  
  // evaluate "true" test statistic ( Pr[ region1 ] / Pr[ region 2 ]
  if( bolTrueStat ) {
    double trueStat;
    trueStat = toyIntProb( X, theta1pri, theta2pri, sig1, sig2, 
                           sigX, l11, u11, l21, u21, n1, n2, N ) /
               toyIntProb( X, theta1pri, theta2pri, sig1, sig2, 
                           sigX, l12, u12, l22, u22, n1, n2, N, "output/PostGrid_seed19.txt" );

    cout << "true statistic = " << trueStat << endl;
    results << trueStat << ", nan, nan, nan, nan, nan, nan, nan" << endl;
  }
  
  // Run SGLD (or MHLD)
  
  // Initialize thetas
  double * thet_t;
  thet_t = new double[2];
  
  // Create array to hold theta draws
  double * thet_draws;
  thet_draws = new double[T*2];
  
  // Create array to hold random data draw
  double * Xi;
  
  // Gradients
  double * lGr;
  lGr = new double[2];
  double * pGr;
  pGr = new double[2];
  
  // Initialize polynomial annealing schedule for epsilon
  tuple< double, double > ab;
  double a;
  double b;
  
  // Autocovariance and autocorrelation time variables
  int lags = 1e+3;
  double * ac;
  ac = new double[lags];
  int M = 4;
  
  // parameters to loop through
  
  // batch sizes
  n[0] = 1; n[1] = 2; n[2] = 4; n[3] = 8; n[4] = 16; n[5] = 100;
  
  // epsMin and max (SGLD)
  epsMax[0] = 1e-2; epsMax[1] = 1e-2; epsMax[2] = 1e-1; epsMax[3] = 1e-1;
  epsMin[0] = 1e-4; epsMin[1] = 1e-5; epsMin[2] = 1e-3; epsMin[3] = 1e-4;
  
  // eps (MHLD)
  eps[0] = 1e-3; eps[1] = 5e-3; eps[2] = 1e-2; eps[3] = 5e-2;
  
  // loop through all specifications of SGLD and MHLD (LMC)
  for( int k = 0; k < 6; k++ ) { // loop through n
    Xi = new double[n[k]];
    for( int j = 0; j < 4; j++ ) { // loop through epsilon
      
      results << "nan" << ", " << n[k] << ", ";
      
      if( k != 5 ) // use annealing schedule for SGLD
        results << epsMin[j] << ", " << epsMax[j] << ", ";
      else // use fixed epsilon for MHLD (LMC)
        results << eps[j] << ", " << 0 << ", ";
      
      // get implied a and b from eps_t = a(b+t)^-gam
      ab = toyEpsilonInit( epsMin[j], epsMax[j], gam, T );
      a = get<0>(ab);
      b = get<1>(ab);
      
      // track rejections
      int rej_ct =  0;
      int all_ct = 0;
  
      // initialize thetas
      thet_t[0] = theta1pri;
      thet_t[1] = theta2pri;
  
      // Sampler
  
      for( int t = 1; t <= T; t++ ) { // generate T draws
  
        if( n[k] < N ) { // SGLD step
          randData( Xi, X, N, n[k] );
          sgldUpdate( thet_t, Xi, lGr, pGr, sig1, sig2, sigX, 
                      theta1pri, theta2pri,  gam, a, b, N, n[k], t );
        }
        else if( n[k] == N ) { // MHLD (LMC) step
          mhldUpdate( thet_t, X, lGr, pGr, sig1, sig2, sigX, 
                      theta1pri, theta2pri,  eps[j], &rej_ct, &all_ct, N, t );
        }
        
        // save draws           
        thet_draws[(t-1)*2]   = thet_t[0];
        thet_draws[(t-1)*2+1] = thet_t[1];
      }
      // output estimate
      results << estStat( thet_draws, l11, u11, l21, u21, l12, u12, l22, u22, T ) << ", ";
      if( n[k] != N ) // if SGLD output weighted estimate
        results << estStatAdj( thet_draws, l11, u11, l21, u21, 
                               l12, u12, l22, u22, a, b, gam, T ) << ", " << "nan" << ", ";  
      else // if MHLD (LMC) output rejection percentage
        results << "nan, " << double( rej_ct ) / double( all_ct ) << ", ";
      
      // calculate autocovariance and autocorrelation time  
      aCovEst( ac, thet_draws, l11, u11, l21, u21, T, lags );
      results << kuboTime( M, ac ) << endl;
    }
  }
  results.close();
}

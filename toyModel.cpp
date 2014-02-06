// toyModel.cpp defines all methods necessary to run SGLD test suite

#include "header.h"

// Generate data from bimodal gaussian with negatively correlated means
void genData( int N,        // Size of dataset
              double X[],   // dataset array
              double th1,   // theta_1 for generating distribution
              double th2,   // theta_2 for generating distribution
              double sX ) { // standard deviation for generating process
  
  double Z;
  double U;
  for( int i = 0; i < N; i++ ) {
    Z = InvNormDist( genrand() );
    U = genrand();
    
    if( U < 0.5 )  // choose which gaussian to sample from
      X[i] = th1 + sX * Z;
    else
      X[i] = th1 + th2 + sX * Z;
  }           
}

// Given max, min epsilon and time steps determine a,b in 
// polynomial annealing schedule: eps_t = a(b+t)^-gam
tuple< double, double > toyEpsilonInit( double epsMin, // max epsilon of annealing schedule
                                        double epsMax, // min epsilon of annealing schedule
                                        double gam,    // polynomial degree
                                        int T ) {      // number of timesteps

  double a;
  double b;
  double epsRatio = epsMin / epsMax;
  
  b = ( T * pow( epsRatio, 1. / gam ) - 1. ) / ( 1. - pow( epsRatio, 1. / gam ) );
  a = epsMax * pow( ( b + 1. ), gam );
  
  return make_tuple( a, b );
}

// return time t epsilon determined by polynomial annealing schedule: a(b+t)^-gam
double toyEpsilon( double a,
                   double b,
                   double gam,
                   int t ) {

  double eps;
  
  eps = a * pow( ( b + t ), -gam );
  return eps;
}
// calculate prior (without normalization)
double toyPrior( double th1,    // theta_1
                 double th2,    // theta_2
                 double th1pri, // theta_1 prior mean
                 double th2pri, // theta_2 prior mean
                 double s1,     // theta_1 prior std dev
                 double s2 ) {  // theta_2 prior std_dev

  double pri;

  pri = //( 1. / ( 2. * M_PI * s1 * s2 ) ) * 
        exp( -pow( ( th1 - th1pri ), 2 ) / ( 2. * pow( s1, 2 ) ) ) *
        exp( -pow( ( th2 - th2pri ), 2 ) / ( 2. * pow( s2, 2 ) ) );

  return pri;
}

// calculate likelihood (without normalization)
double toyLikelihood( double X[], // data
                      double th1, // theta_1 draw
                      double th2, // theta_2 draw
                      double sX,  // data std dev
                      int N ) {   // size of data
                      
  double likh = 1.;
  
  double C = 1. / ( pow( ( 2. * M_PI ), 0.5 ) * sX ); // Don't use constants since
                                                      // it leads to numerical underflow

  for( int i = 0; i < N; i++ ) {
    likh *= //( C / 2. ) * 
            ( exp( - pow( ( X[i] - th1 ), 2 ) / ( 2 * pow( sX, 2 ) ) ) +
              exp( - pow( ( X[i] - th1 - th2 ), 2 ) / ( 2 * pow( sX, 2 ) ) ) );
  }

  return likh;
}
// calculate posterior (without normalization)
double toyPosterior( double X[],    // data
                     double th1,    // theta_1 draw
                     double th2,    // theta_1 draw
                     double th1pri, // theta_1 prior mean
                     double th2pri, // theta_2 prior mean
                     double s1,     // theta_1 prior std dev
                     double s2,     // theta_2 prior std dev
                     double sX,     // data std dev
                     int N ) {      // size of data

  double pri; // calculate prior
  pri = toyPrior( th1, th2, th1pri, th2pri, s1, s2 );

  double likh; // calculate likelihood
  likh = toyLikelihood( X, th1, th2, sX, N );

  double post;  // calculate posterior
  post = pri * likh;

  return post;  
}

// integrate posterior (un-normalized) over some region
double toyIntProb( double X[],  // data
                   double th1p, // theta_1 prior mean
                   double th2p, // theta_2 prior mean
                   double s1,   // theta_1 prior std dev
                   double s2,   // theta_2 prior std dev
                   double sX,   // data std dev
                   double l1,   // theta_1 LB of region
                   double u1,   // theta_1 UB of region
                   double l2,   // theta_2 LB of region
                   double u2,   // theta_2 UB of region
                   int n1,      // number of theta_1 grid points
                   int n2,      // number of theta_2 grid points
                   int N,       // size of data
                   string fn ) { // file name to which posterior mesh is output

  double h1;
  h1 = ( u1 - l1 ) / ( n1 - 1 );
  
  double h2;
  h2 = ( u2 - l2 ) / ( n2 - 1 );
  
  double val = 0;
  
  ofstream writeFile;
  
  if( fn.length() != 0 )
    writeFile.open( fn ); // open posterior mesh save file

  for( int i = 1; i < n1; i++ ) {   // loop over theta_1 grid points
    for( int j = 1; j < n2; j++ ) { // loop over theta_2 grid points
      val += ( h1 * h2 / 4. ) * ( toyPosterior( X, l1+i*h1, l2+j*h2, th1p, th2p, s1, s2, sX, N )
             + toyPosterior( X, l1+i*h1, l2+(j-1)*h2, th1p, th2p, s1, s2, sX, N )
             + toyPosterior( X, l1+(i-1)*h1, l2+j*h2, th1p, th2p, s1, s2, sX, N )
             + toyPosterior( X, l1+(i-1)*h1, l2+(j-1)*h2, th1p, th2p, s1, s2, sX, N ) );
      if( fn.length() != 0 ) {
        // save point to posterior mesh file
        writeFile << l1+i*h1 << "," << l2+j*h2 << "," <<
          toyPosterior( X, l1+i*h1, l2+j*h2, th1p, th2p, s1, s2, sX, N ) << endl;
      }    
    }
  }
  if( fn.length() != 0 )  
    writeFile.close();
  
  return val;
}

// Calculate gradient of prior
void toyPriorGrad( double pGr[], // prior gradient vector (theta_1, theta_2)
                   double th1,   // theta draws
                   double th2,
                   double th1pri, // theta prior means
                   double th2pri,
                   double s1,     // theta prior std dev
                   double s2 ) {

  pGr[0] = 0;
  pGr[1] = 0;
  
  pGr[0] = -( th1 - th1pri ) / pow( s1, 2 );
  pGr[1] = -( th2 - th2pri ) / pow( s2, 2 );
}
// Calculate gradient of likelihood
void toyLikelihoodGrad( double X[],   // data (batch for SGLD, full data for MHLD)
                        double lGr[], // likelihood gradient vector (theta_1, theta_2)
                        double th1,   // theta draws
                        double th2,  
                        double sX,    // data std dev
                        int n,        // size of batch 
                        int N ) {     // size of data
  
  lGr[0] = 0;
  lGr[1] = 0;

  for( int i = 0; i < n; i++ ) { // messy algebra
    lGr[0] += ( ( X[i] - th1 ) / pow( sX, 2 ) *
              exp( - pow( ( X[i] - th1 ), 2 ) / ( 2 * pow( sX, 2 ) ) ) +
              ( X[i] - th1 - th2 ) / pow( sX, 2 ) *
              exp( - pow( ( X[i] - th1 - th2 ), 2 ) / ( 2 * pow( sX, 2 ) ) ) ) /
              ( exp( - pow( ( X[i] - th1 ), 2 ) / ( 2 * pow( sX, 2 ) ) ) +
                exp( - pow( ( X[i] - th1 - th2 ), 2 ) / ( 2 * pow( sX, 2 ) ) ) );

    
    lGr[1] += ( X[i] - th1 - th2 ) / pow( sX, 2 ) *
              exp( - pow( ( X[i] - th1 - th2 ), 2 ) / ( 2 * pow( sX, 2 ) ) ) /
              ( exp( - pow( ( X[i] - th1 ), 2 ) / ( 2 * pow( sX, 2 ) ) ) +
                exp( - pow( ( X[i] - th1 - th2 ), 2 ) / ( 2 * pow( sX, 2 ) ) ) );
  }
  
  lGr[0] *= double( N ) / double( n );
  lGr[1] *= double( N ) / double( n );
  
}
// SGLD proposal step
void sgldUpdate( double th[],   // theta draw
                 double X[],    // data or batch of data
                 double lGr[],  // lik. gradient
                 double pGr[],  // pri. gradient
                 double s1,     // theta prior std dev
                 double s2,
                 double sX,     // data std dev
                 double th1pri, // theta prior means
                 double th2pri,  
                 double gam,    // annealing parameters
                 double a,
                 double b,
                 int N,         // size of data
                 int n,         // size of batch
                 int t ) {      // current step (time)
  
  // Update Likelihood Gradient              
  toyLikelihoodGrad( X, lGr, th[0], th[1], sX, n, N );
  // Update Prior Gradient
  toyPriorGrad( pGr, th[0], th[1], th1pri, th2pri, s1, s2 );
  
  // Injected noise
  double eps; 
  eps = toyEpsilon( a, b, gam, t ); 

  double eta1, eta2; // draw noise
  eta1 = InvNormDist( genrand() ) * pow( eps, 0.5 );
  eta2 = InvNormDist( genrand() ) * pow( eps, 0.5 );
  
  // SGLD step
  th[0] += eps * 0.5 * ( pGr[0] + lGr[0] ) + eta1;
  th[1] += eps * 0.5 * ( pGr[1] + lGr[1] ) + eta2;     
}  

// MHLD (LMC) proposal step
void mhldUpdate( double th[],   // theta draws
                 double X[],    // data
                 double lGr[],  // lik. gradient
                 double pGr[],  // pri. gradient
                 double s1,     // theta prior std dev
                 double s2,
                 double sX,     // data std dev 
                 double th1pri, // theta prior means
                 double th2pri,  
                 double eps,    // step size
                 int *r_ct,     // pointer to rejection count
                 int *a_ct,     // pointer to proposal count
                 int N,         // data size
                 int t ) {      // current step (time)
  
  double eta1, eta2;
  double pr_a;
  
  double * tmp_th;
  tmp_th = new double[2];
  
  double * tmp_lGr;
  tmp_lGr = new double[2];
  
  double * tmp_pGr;
  tmp_pGr = new double[2];
  

  // Update Likelihood Gradient              
  toyLikelihoodGrad( X, lGr, th[0], th[1], sX, N, N );
  // Update Prior Gradient
  toyPriorGrad( pGr, th[0], th[1], th1pri, th2pri, s1, s2 );
    
  // Draw noise
  eta1 = InvNormDist( genrand() ) * pow( eps, 0.5 );
  eta2 = InvNormDist( genrand() ) * pow( eps, 0.5 );
    
  // save current thetas and gradients (in case of rejection)
  tmp_th[0] = th[0];
  tmp_th[1] = th[1];
    
  tmp_lGr[0] = lGr[0];
  tmp_lGr[1] = lGr[1];
    
  tmp_pGr[0] = pGr[0];
  tmp_pGr[1] = pGr[1];
    
  // SGLD step
  th[0] += eps * 0.5 * ( pGr[0] + lGr[0] ) + eta1;
  th[1] += eps * 0.5 * ( pGr[1] + lGr[1] ) + eta2;
    
  // acceptance probability
  pr_a = toyPosterior( X, th[0], th[1], th1pri, th2pri, s1, s2, sX, N ) /
         toyPosterior( X, tmp_th[0], tmp_th[1], th1pri, th2pri, s1, s2, sX, N );
    
  // Update Likelihood Gradient              
  toyLikelihoodGrad( X, lGr, th[0], th[1], sX, N, N );
  // Update Prior Gradient
  toyPriorGrad( pGr, th[0], th[1], th1pri, th2pri, s1, s2 );
    
  // update acceptance probability (since proposal is not symmetric)
  pr_a *= exp( -( pow( eta1 + eps * 0.5 * 
                      ( tmp_lGr[0] + tmp_pGr[0] + lGr[0] + pGr[0] ), 2) +
                  pow( eta2 + eps * 0.5 * 
                       ( tmp_lGr[1] + tmp_pGr[1] + lGr[1] + pGr[1] ), 2) -
                  pow( eta1, 2) - pow( eta2, 2 ) ) / ( 2 * eps ) );
                                                 
  if( genrand() > min( pr_a, 1. ) ) { // go back to previous theta
    th[0] = tmp_th[0];
    th[1] = tmp_th[1];
    *r_ct += 1; // count rejections
  }
  *a_ct += 1; // count total proposals
}

// generate random data without replacement (Better results in practice, SOURCE: Niu11)
void randData( double Xi[], // batch data array to be populated
               double X[],  // full data array
               int N,
               int n ) {
  int max = N;
  int ind;
  double tmp;

  for( int i = 0; i < n; i++ ) {
    ind = int( genrand() * max ); // random draw
    Xi[i] = X[ind]; // add to batch
    // switch sampled point with the last point in X
    tmp = X[max-1]; 
    X[max-1] = X[ind];
    X[ind] = tmp;
    // decrement max so as not to resample previously sampled point
    max -= 1;
  }
  
}
// calculate estimated statistic ( Pr[Region] )
double estStat( double th_draws[], // theta draws
              double l11, // these are the bounds of the region of interest
              double u11,
              double l21,
              double u21,
              double l12, // these are the bounds of the large region used as an estimate
              double u12, // for the normalizing constant of the distribution
              double l22,
              double u22,
              int T ) {
   
  int ct1 = 0;
  int ct2 = 0;
   
   for( int t = 0; t < T; t++ ) {
    if( th_draws[t*2] <= u11 && th_draws[t*2] >= l11 
        && th_draws[t*2+1] <= u21 && th_draws[t*2+1] >= l21 )
        ct1 += 1; // count draws that fall in region of interest
     if( th_draws[t*2] <= u12 && th_draws[t*2] >= l12 
         && th_draws[t*2+1] <= u22 && th_draws[t*2+1] >= l22 )
        ct2 += 1; // count draws that fall in our approximate support for distribution
   }
   cout << "count 1 = " << ct1 << endl;
   cout << "count 2 = " << ct2 << endl;
   // estimated statistic is essentially a conditional probability
   cout << "estimated statistic = " << double( ct1 ) / double( ct2 ) << endl; 
   return double( ct1 ) / double( ct2 ); 
}   

// Calculates estimated statistic adjusted using epsilon_t weights (same structure as above)
double estStatAdj( double th_draws[],
              double l11,
              double u11,
              double l21,
              double u21,
              double l12,
              double u12,
              double l22,
              double u22,
              double a,
              double b,
              double gam,
              int T ) {
   
  double sum1 = 0;
  double sum2 = 0;
   
   for( int t = 0; t < T; t++ ) {
    if( th_draws[t*2] <= u11 && th_draws[t*2] >= l11 
        && th_draws[t*2+1] <= u21 && th_draws[t*2+1] >= l21 )
        sum1 += toyEpsilon( a, b, gam, t+1 );
     if( th_draws[t*2] <= u12 && th_draws[t*2] >= l12 
         && th_draws[t*2+1] <= u22 && th_draws[t*2+1] >= l22 )
        sum2 += toyEpsilon( a, b, gam, t+1 );
   }
   cout << sum1 << endl;
   cout << sum2 << endl;
   cout << "adjusted estimated statistic = " << sum1 / sum2 << endl; 
   return sum1 / sum2;            
}

double kuboTime(    // calculate self-consistent auto-correlation time
  double M,         // M = 4
  double ac[] ) {   // vector of calculated autocovariances, C(t)

  double kubo = 1.;
  int W = 1;
  
  while( W < kubo * M ) {  // loop until self consistent
    kubo += 2 * ac[W++] / ac[0]; // add 2 * autocorrelation at time W
  }
  
  if ( W < kubo * M ) // if we run out of autocorrelations
  	cout << "WARNING: tau is inaccurate... not enough autocovariances " << endl;
  
  return kubo;

}


// aCovEst(...) calculates the estimated auto-covariance for 
// the series Y for lags t = 0,1,...,n

void aCovEst( double ac[], // pass in empty auto-cov array
	          double th_draws[],  // pass in series of interest
	          double l11,
              double u11,
              double l21,
              double u21,
	          int T,       // length of series
	          int n ) {    // max number of lags for which to
                           // compute auto-cov
  
  double * I;
  I = new double[T];
  
  for( int t = 0; t < T; t++ ) {
    if( th_draws[t*2] <= u11 && th_draws[t*2] >= l11 
        && th_draws[t*2+1] <= u21 && th_draws[t*2+1] >= l21 )
        I[t] = 1;
    else
        I[t] = 0;
   }

  double Ibar = 0.;
  int i, j; 

  for( i = 0; i < T ; i ++ ) {  // calculate mean of series, Ybar
    Ibar += I[i];
  }

  Ibar = Ibar / T;

  for( i = 0; i < n; i++ ) {    // loop through lags
    ac[i] = 0;

    for( j = 0; j < T - i; j++ ) {  // calculate auto-cov for lag i
      ac[i] += ( I[j] - Ibar ) * ( I[j+i] - Ibar );
    }

    ac[i] = ac[i] / ( T - i - 1 );  // account for degree of freedom lost
                                    // from using Ybar instead of true mean
  } 
}
  
  

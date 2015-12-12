//#include <stdlib.h>
//#include <math.h>
//#include <stdio.h>
//#include <sys\timeb.h>

#include "randgen.h"
#include "birthdeath_strauss.h"

extern double M;

/*
 * kiss
 */
/*{{{*/
double kiss(unsigned long *seed) {

  /* Generator proposed by Marsaglia and Zaman, 1993. See Robert and Casella
   * (1999, pages 41-43) for details.  
   * Watch out: the last line
   *      x = ((double) (*i+*j+*k)*exp(-32*log(2.0)));
   * must be calibrated, depending on the precision of the computer. */

  double x;

  seed[1] = seed[1] ^ (seed[1]<<17);
  seed[2] = (seed[2] ^ (seed[2]<<18)) & 0x7FFFFFFF;
  seed[0] = 69069*seed[0]+23606797;
  seed[1] ^= (seed[1]>>15);
  seed[2] ^= (seed[2]>>13);
  x = ((double) (seed[0]+seed[1]+seed[2])*M);

  return(x);

}
/*}}}*/




/*
 * runiform_n
 */
/*{{{*/
int runiform_n(int n, unsigned long *seed) {

  int i;
  double *prob;

  prob = (double *)calloc(n, sizeof(double));

  for (i=0;i<n;i++) {
    prob[i] = (i+1)/(double)n;
  }

  i = (int)rmultinomial(prob, n, seed);
  free(prob);

  return i;

}
/*}}}*/




/* 
 * rmultinomial
 */
/*{{{*/
long rmultinomial(double *prob, long len, unsigned long *seed) {

  long i;
  double y;

  y = kiss(seed);

  for (i=0;i<len;i++) {
    if (y <= prob[i]) {
      return i;
    }
  }

  return len-1;

}
/*}}}*/




double runif_atob(unsigned long *seed, double a, double b) {
  return ((b-a)*kiss(seed) + a);
}

double rexp(double beta, unsigned long *seed) {
  return -log(kiss(seed))/beta;
}




/*
 * sgamma
 */
/*{{{*/
double sgamma(double a, unsigned long* seed) {
  /*
   *
   * (STANDARD-)  GAMMA  DISTRIBUTION
   *
   * -----------------------------------------------------------------------------
   *
   * PARAMETER  A >= 1.0  !
   * 
   * For details see:
   * 
   *           AHRENS, J.H. AND DIETER, U.
   *           GENERATING GAMMA VARIATES BY A
   *           MODIFIED REJECTION TECHNIQUE.
   *           COMM. ACM, 25,1 (JAN. 1982), 47 - 54.
   * 
   * Step numbers correspond to algorithm 'GD' in the above paper
   *                             (straightforward implementation)
   * 
   * Modified by Barry W. Brown, Feb 3, 1988 to use RANF instead of SUNIF.  The
   * argument IR thus goes away.
   * 
   * -----------------------------------------------------------------------------
   *
   * PARAMETER  0.0 < A < 1.0  !
   *
   * For details see:
   *
   *           AHRENS, J.H. AND DIETER, U.
   *           COMPUTER METHODS FOR SAMPLING FROM GAMMA,
   *           BETA, POISSON AND BINOMIAL DISTRIBUTIONS.
   *           COMPUTING, 12 (1974), 223 - 246.
   *
   * (Adapted implementation of algorithm 'GS' in the above paper)
   *
   * -----------------------------------------------------------------------------
   * 
   * Input: a =parameter (mean) of the standard gamma distribution
   * Output: sgamma = sample from the gamma-(a)-distribution
   * Coefficients q(k) - for q0 = sum(q(k)*a**(-k))
   * Coefficients a(k) - for q = q0+(t*t/2)*sum(a(k)*v**k)
   * Coefficients e(k) - for exp(q)-1 = sum(e(k)*q**k)
   * Previous a pre-set to zero - aa is a', aaa is a"
   * Sqrt32 is the squareroot of 32 = 5.656854249492380
   *
   */

  extern double fsign( double num, double sign );
  static double q1 = 4.166669E-2;
  static double q2 = 2.083148E-2;
  static double q3 = 8.01191E-3;
  static double q4 = 1.44121E-3;
  static double q5 = -7.388E-5;
  static double q6 = 2.4511E-4;
  static double q7 = 2.424E-4;
  static double a1 = 0.3333333;
  static double a2 = -0.250003;
  static double a3 = 0.2000062;
  static double a4 = -0.1662921;
  static double a5 = 0.1423657;
  static double a6 = -0.1367177;
  static double a7 = 0.1233795;
  static double e1 = 1.0;
  static double e2 = 0.4999897;
  static double e3 = 0.166829;
  static double e4 = 4.07753E-2;
  static double e5 = 1.0293E-2;
  static double aa = 0.0;
  static double aaa = 0.0;
  static double sqrt32 = 5.656854;
  static double sgamma,s2,s,d,t,x,u,r,q0,b,si,c,v,q,e,w,p;
  double snorm(unsigned long *);
  double kiss(unsigned long *);
  double rexp(double,unsigned long *);

  if(a == aa) goto S10;
  if(a < 1.0) goto S120;

  /*
   * Step  1:  recalculations of s2,s,d if a has changed
   */
  aa = a;
  s2 = a-0.5;
  s = sqrt(s2);
  d = sqrt32-12.0*s;

S10:
  /*
   * Step  2:  t=standard normal deviate,
   *           x=(s,1/2)-normal deviate.
   *           immediate acceptance (i)
   */
  t = snorm(seed);
  x = s+0.5*t;
  sgamma = x*x;
  if(t >= 0.0) return sgamma;

  /*
   * Step  3:  u= 0,1 -uniform sample. squeeze acceptance (s)
   */
  u = kiss(seed);
  if(d*u <= t*t*t) return sgamma;

  /*
   *  Step  4:  recalculations of q0,b,si,c if necessary
   */
  if(a == aaa) goto S40;
  aaa = a;
  r = 1.0/ a;
  q0 = ((((((q7*r+q6)*r+q5)*r+q4)*r+q3)*r+q2)*r+q1)*r;

  /*
   * Approximation depending on size of parameter a the constants in the
   * expressions for b, si and c were established by numerical experiments
   */
  if(a <= 3.686) goto S30;
  if(a <= 13.022) goto S20;

  /*
   * Case 3:  a .gt. 13.022
   */
  b = 1.77;
  si = 0.75;
  c = 0.1515/s;
  goto S40;

S20:
  /*
   * Case 2:  3.686 .lt. a .le. 13.022
   */
  b = 1.654+7.6E-3*s2;
  si = 1.68/s+0.275;
  c = 6.2E-2/s+2.4E-2;
  goto S40;

S30:
  /*
   * Case 1:  a .le. 3.686
   */
  b = 0.463+s+0.178*s2;
  si = 1.235;
  c = 0.195/s-7.9E-2+1.6E-1*s;

S40:
  /*
   * Step  5:  No quotient test if x not positive
   */
  if(x <= 0.0) goto S70;

  /*
   * Step  6:  Calculation of v and quotient q
   */
  v = t/(s+s);
  if(fabs(v) <= 0.25) goto S50;
  q = q0-s*t+0.25*t*t+(s2+s2)*log(1.0+v);
  goto S60;

S50:
  q = q0+0.5*t*t*((((((a7*v+a6)*v+a5)*v+a4)*v+a3)*v+a2)*v+a1)*v;

S60:
  /*
   * Step  7:  Quotient acceptance (q)
   */
  if(log(1.0-u) <= q) return sgamma;

S70:
  /*
   * Step  8:  e = standard exponential deviate
   *           u =  0,1 - uniform deviate
   *           t = (b,si) - double exponential (laplace) sample
   */
  e  = rexp(1.0, seed);
  u  = kiss(seed);
  u += (u - 1.0);
  t  = b + fsign(si * e, u);

  /*
   * Step  9:  Rejection if t .lt. tau(1) = -.71874483771719
   */
  if(t < -0.7187449) goto S70;

  /*
   * Step 10:  Calculation of v and quotient q
   */
  v = t/(s+s);
  if(fabs(v) <= 0.25) goto S80;
  q = q0-s*t+0.25*t*t+(s2+s2)*log(1.0+v);
  goto S90;

S80:
  q = q0+0.5*t*t*((((((a7*v+a6)*v+a5)*v+a4)*v+a3)*v+a2)*v+a1)*v;

S90:
  /*
   * Step 11:  Hat acceptance (h) (if q not positive go to step 8)
   */
  if(q <= 0.0) goto S70;
  if(q <= 0.5) goto S100;
  w = exp(q)-1.0;
  goto S110;

S100:
  w = ((((e5*q+e4)*q+e3)*q+e2)*q+e1)*q;

S110:
  /*
   * If t is rejected, sample again at step 8
   */
  if(c*fabs(u) > w*exp(e-0.5*t*t)) goto S70;
  x = s+0.5*t;
  sgamma = x*x;
  return sgamma;

S120:
  /*
   * Alternate method for parameters a below 1  (.3678794=exp(-1.))
   */
  aa = 0.0;
  b = 1.0+0.3678794*a;

S130:
  p = b*kiss(seed);
  if(p >= 1.0) goto S140;
  sgamma = exp(log(p)/ a);
  if(rexp(1.0,seed) < sgamma) goto S130;
  return sgamma;

S140:
  sgamma = -log((b-p)/ a);
  if(rexp(1.0,seed) < (1.0-a)*log(sgamma)) goto S130;
  return sgamma;

}
/*}}}*/




/* 
 * fsign
 */
double fsign( double num, double sign ) {
  /* Transfers sign of argument sign to argument num */
  if ((sign > 0.0f && num < 0.0f) || (sign < 0.0f && num > 0.0f))
    return -num;
  else return num;
}




/*
 * snorm
 */
/*{{{*/
double snorm(unsigned long *seed)
  /*
   **********************************************************************


   (STANDARD-)  N O R M A L  DISTRIBUTION


   **********************************************************************
   **********************************************************************

   FOR DETAILS SEE:

   AHRENS, J.H. AND DIETER, U.
   EXTENSIONS OF FORSYTHE'S METHOD FOR RANDOM
   SAMPLING FROM THE NORMAL DISTRIBUTION.
   MATH. COMPUT., 27,124 (OCT. 1973), 927 - 937.

   ALL STATEMENT NUMBERS CORRESPOND TO THE STEPS OF ALGORITHM 'FL'
   (M=5) IN THE ABOVE PAPER     (SLIGHTLY MODIFIED IMPLEMENTATION)

   Modified by Barry W. Brown, Feb 3, 1988 to use RANF instead of
   SUNIF.  The argument IR thus goes away.

   **********************************************************************
   THE DEFINITIONS OF THE CONSTANTS A(K), D(K), T(K) AND
   H(K) ARE ACCORDING TO THE ABOVEMENTIONED ARTICLE
   */
{
  static double a[32] = {
    0.0,3.917609E-2,7.841241E-2,0.11777,0.1573107,0.1970991,0.2372021,0.2776904,
    0.3186394,0.36013,0.4022501,0.4450965,0.4887764,0.5334097,0.5791322,
    0.626099,0.6744898,0.7245144,0.7764218,0.8305109,0.8871466,0.9467818,
    1.00999,1.077516,1.150349,1.229859,1.318011,1.417797,1.534121,1.67594,
    1.862732,2.153875
  };
  static double d[31] = {
    0.0,0.0,0.0,0.0,0.0,0.2636843,0.2425085,0.2255674,0.2116342,0.1999243,
    0.1899108,0.1812252,0.1736014,0.1668419,0.1607967,0.1553497,0.1504094,
    0.1459026,0.14177,0.1379632,0.1344418,0.1311722,0.128126,0.1252791,
    0.1226109,0.1201036,0.1177417,0.1155119,0.1134023,0.1114027,0.1095039
  };
  static double t[31] = {
    7.673828E-4,2.30687E-3,3.860618E-3,5.438454E-3,7.0507E-3,8.708396E-3,
    1.042357E-2,1.220953E-2,1.408125E-2,1.605579E-2,1.81529E-2,2.039573E-2,
    2.281177E-2,2.543407E-2,2.830296E-2,3.146822E-2,3.499233E-2,3.895483E-2,
    4.345878E-2,4.864035E-2,5.468334E-2,6.184222E-2,7.047983E-2,8.113195E-2,
    9.462444E-2,0.1123001,0.136498,0.1716886,0.2276241,0.330498,0.5847031
  };
  static double h[31] = {
    3.920617E-2,3.932705E-2,3.951E-2,3.975703E-2,4.007093E-2,4.045533E-2,
    4.091481E-2,4.145507E-2,4.208311E-2,4.280748E-2,4.363863E-2,4.458932E-2,
    4.567523E-2,4.691571E-2,4.833487E-2,4.996298E-2,5.183859E-2,5.401138E-2,
    5.654656E-2,5.95313E-2,6.308489E-2,6.737503E-2,7.264544E-2,7.926471E-2,
    8.781922E-2,9.930398E-2,0.11556,0.1404344,0.1836142,0.2790016,0.7010474
  };
  static long i;
  static double snorm,u,s,ustar,aa,w,y,tt;
  u = kiss(seed);
  s = 0.0;
  if(u > 0.5) s = 1.0;
  u += (u-s);
  u = 32.0*u;
  i = (long) (u);
  if(i == 32) i = 31;
  if(i == 0) goto S100;
  /*
     START CENTER
     */
  ustar = u-(double)i;
  aa = *(a+i-1);
S40:
  if(ustar <= *(t+i-1)) goto S60;
  w = (ustar-*(t+i-1))**(h+i-1);
S50:
  /*
     EXIT   (BOTH CASES)
     */
  y = aa+w;
  snorm = y;
  if(s == 1.0) snorm = -y;
  return snorm;
S60:
  /*
     CENTER CONTINUED
     */
  u = kiss(seed);
  w = u*(*(a+i)-aa);
  tt = (0.5*w+aa)*w;
  goto S80;
S70:
  tt = u;
  ustar = kiss(seed);
S80:
  if(ustar > tt) goto S50;
  u = kiss(seed);
  if(ustar >= u) goto S70;
  ustar = kiss(seed);
  goto S40;
S100:
  /*
     START TAIL
     */
  i = 6;
  aa = *(a+31);
  goto S120;
S110:
  aa += *(d+i-1);
  i += 1;
S120:
  u += u;
  if(u < 1.0) goto S110;
  u -= 1.0;
S140:
  w = u**(d+i-1);
  tt = (0.5*w+aa)*w;
  goto S160;
S150:
  tt = u;
S160:
  ustar = kiss(seed);
  if(ustar > tt) goto S50;
  u = kiss(seed);
  if(ustar >= u) goto S150;
  u = kiss(seed);
  goto S140;
}/*}}}*/





double rnorm(double mean, double stdev, unsigned long *seed) {
  double snorm(unsigned long *);
  return mean + stdev*snorm(seed);
}

double rgamma(double alpha, double beta, unsigned long *seed) {
  double sgamma(double, unsigned long *);
  return sgamma(alpha, seed)/beta;
}

double inverse_gamma(double alpha, double beta, unsigned long *seed) {
  double sgamma(double, unsigned long *);
  return beta/sgamma(alpha, seed);
}

double fact_ln(int k) {
  int i;
  double fact;
  if (k == 0) return 0;
  fact = 0.0;
  for (i=1;i<=k;i++) fact += log((double)i);
  return fact;
}

int rpois(double lambda, unsigned long *seed) {
  int i;
  double fact_ln(int);
  double kiss(unsigned long *);
  double U, P;

  U = kiss(seed);
  i = 0;

  P = exp(-lambda);
  while (P <= U) {
    i++;
    P += exp(-lambda + i * log(lambda) - fact_ln(i));
  }
  return i;
}

double *rdirichlet(double *alpha, int len, unsigned long *seed) {
  int i;
  double *theta, denom=0.0;
  double sgamma(double, unsigned long *);

  theta = (double *)calloc(len, sizeof(double));
  for (i=0;i<len;i++) {
    theta[i] = sgamma(alpha[i], seed);
    denom += theta[i];
  }
  for (i=0;i<len;i++) theta[i] /= denom;

  return theta;
}

double rbeta(double alpha, double beta, unsigned long *seed) {
  double a, b;
  double sgamma(double, unsigned long *);

  a = sgamma(alpha, seed);
  b = sgamma(beta, seed);

  return a/(a+b);
}


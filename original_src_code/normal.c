#include "deconvolution_main.h"

#define EPS 1.0e-13

static double C = 0.1591549;
//static double C0 = 0.7071068;
//static double C1 = 0.999999;
//static double C2 = 6.5;
//static double C3 = 1.193e-07;
static double EXPOV = 174.673;


/*
 * MDTNF
 */
/*{{{*/
double MDTNF(double y, double z, double eps, int *nnn) {

  int    i, j, f, g, g1;
  double t=0.0, ta, a, b, hsqb, ber, ter, sum, aeps;
  double bexp, asq, a4, b4, a4b4, ahsqb, ab4, d2, d;
  double prob_std_norm();

  *nnn = 0;

  if (eps == 0) eps = EPS;

  b = fabs(y);
  a = fabs(z);

  if (a >= 0.0) {

    ta = atan(a);

    if (a*b > 4.0) {

      t = prob_std_norm(b,&d);
      t = C*(ta + atan(1.0/a)) - 0.5*(t-0.5);

    } else {

      hsqb = 0.5*b*b;
      if (hsqb > EXPOV) return t;

      bexp = exp(-hsqb);
      asq = a*a;
      a4 = asq*asq;
      b4 = hsqb*hsqb;
      a4b4 = a4*b4;
      ahsqb = a*hsqb;
      ab4 = a*b4*0.5;
      f = 1;
      sum = 0;
      g = 3;

      for (j=0;j<100;j++) {
        *nnn = j;
        g1 = g;
        ber = ter  = ab4;
        for (i=0;;i++) {
          if (ter > ber*eps) {
            ter *= hsqb/(double)g1;
            g1++;
            ber += ter;
          } else {
            break;
          }
        }

        sum += (ber + ahsqb)/(double)f - ber*asq/(double)(f+2);

        t = ta - sum*bexp;

        aeps = eps*t;

        ahsqb *= a4b4/(((double)g-1.0)*(double)g);
        ab4   *= a4b4/(((double)g+1.0)*(double)g);

        d2 = ber*asq/(double)(f+2);

        f += 4;
        g += 2;

        if (d2*bexp < aeps) {	  
          t = t*C;
          break;
        }
      }
    }
  } else {
    return t;
  }

  t = (z < 0.0) ? -t:t;

  return t;

}
/*}}}*/




///*
// * MDBNOR
// */
///*{{{*/
//double MDBNOR(double x, double y, double rho) {
//
//    int ind, iax, iay, i, nnn, ier;
//    double p=0, d, ay, ax, yyy, xxx, qx, qy, tx, ty, f1;
//    void UERTST();
//    double prob_std_norm();
//
//    if (fabs(rho) >= C1) {
//        UERTST(129, "6HMDBNOR");
//        return p;
//    }
//
//    if ((x > C2) && (y > C2)) {
//        p = 1;
//        return p;
//    } else if ((x < -C2) && (y < -C2)) {
//        p = 0;
//        return p;
//    } else if ((x > C2) && (C2 < y)) {
//        p = prob_std_norm(y, &d);
//        if (p < C3) p = 0.0;
//        if ((1-p) < C3) p = 1.0;
//        return p;
//    } else if ((y > C2) && (C2 < x)) {
//        p = prob_std_norm(x, &d);
//        if (p < C3) p = 0.0;
//        if ((1-p) < C3) p = 1.0;
//        return p;
//    }
//    
//    f1  = 1.0 / sqrt(1.0 - rho * rho);
//
//    ind = 0;
//    iax = 0;
//    iay = 0;
//
//    if (x*y == 0) {
//        if (x != 0) {
//            ty = 0.25;
//            if (x < 0) ty = -ty;
//            ax = -f1*rho;
//            ind = 1;
//        } else if  (y != 0) { 
//            tx = 0.25;
//            if (y < 0) tx = -tx;
//            ay = -f1*rho;
//        } else {
//            ay = ax = f1*(1.0 - rho);
//        }
//    } else {
//        ax = f1*(y/x - rho);
//        ay = f1*(x/y - rho);
//    }
//
//    if ((x*y != 0) || (x !=0) || ((x == 0) && (y == 0))) {
//        if (ax < 0) {
//            iax = 1;
//            ax = -ax;
//        }
//        
//        xxx = x;
//
//        for (i=0; i<30; i++) {
//            tx = MDTNF(xxx, ax, 0.0, &nnn);
//            if (nnn >= 99) {
//                Rprintf(".........xxx = %lf\n", xxx);
//                if (xxx > 0) xxx -= 0.000001;
//                else         xxx += 0.000001;
//            } else {
//                break;
//            }
//        }
//
//        if (i == 30) {
//            ier = 1000;
//            return 0.0;
//        }
//
//        if (iax > 0) tx = -tx;
//    }
//
//    if (ind == 0) {
//
//        if (ay < 0) {
//            iay = 1;
//            ay = -ay;
//        }
//        yyy = y;
//
//        for (i=0;i<30;i++) {
//            ty = MDTNF(yyy,ay,0.0,&nnn);
//            if (nnn >= 99) {
//                Rprintf(".........yyy = %lf\n",y);
//                if (yyy > 0) yyy -= 0.000001;
//                else         yyy += 0.000001;
//            } else {
//                break;
//            }
//        }
//
//        if (i == 30) {
//            ier = 1000;
//            return 0.0;
//        }
//        if (iay > 0) ty = -ty;
//
//    }
//
//    qx = prob_std_norm(x, &d);
//    qy = prob_std_norm(y, &d);
//
//    p = 0.5 * (qx + qy) - tx - ty;
//
//    if ((x * y <= 0) && ((x * y != 0) || (x + y < 0))) p -= 0.5;
//
//    p = (0.0 < p) ? p:0.0;
//    p = (p < 1.0) ? p:1.0;
//
//    if (p < C3) p = 0.0;
//    if ((1 - p) < C3) p = 1.0;
//
//    return p;
//
//}
///*}}}*/




/*
 * UERTST
 */
/*{{{*/
void UERTST(ier, c)
  int ier;
  char *c;
{
  if (ier < 32)
    Rprintf(" *** I M S L (UERTST) *** NON-DEFINED  %s   %d  (ier = %d)\n",
           c, ier-0, ier);
  else if (ier < 64)
    Rprintf(" *** I M S L (UERTST) *** WARNING  %s   %d  (ier = %d)\n", c,
           ier-32, ier);
  else if (ier < 128)
    Rprintf(" *** I M S L (UERTST) *** WARNING(with fix)  %s   %d  (ier = %d)\n", 
           c, ier-64, ier);    
  else
    Rprintf(" *** I M S L (UERTST) *** TERMINAL  %s   %d  (ier = %d)\n", c,
           ier-128, ier);    

}
/*}}}*/




/*
 * prob_std_norm
 */
/*{{{*/
double prob_std_norm(double z, double *dens) {

  double abs_z;
  double t, prob;

  abs_z = fabs(z);
  t     = 1.0 / (1.0 + 0.2316419 * abs_z);

  if (abs_z <= 18) {
    *dens = 0.3989423*exp(-z*z/2.0);
  } else {
    *dens = 0.0;
  }

  prob = 1.0 - (*dens) * t * 
    ((((1.330274 * t - 1.821256) * t + 1.781478) * t - 0.3565638) * 
     t + 0.3193815);

  if (z < 0) prob = 1.0 - prob;

  return prob;

}
/*}}}*/




/*******************************************************************************
 * END OF FILE
 ******************************************************************************/


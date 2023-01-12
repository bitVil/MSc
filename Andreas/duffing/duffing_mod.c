/*
!---------------------------------------------------------------------- 
!---------------------------------------------------------------------- 
!   Duffing Oscillator (Parlitz and Lauterborn 1985)
!---------------------------------------------------------------------- 
!---------------------------------------------------------------------- 
*/

#include "auto_f2c.h"
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_errno.h>


enum {x1i, x2i, l1i, l2i};
enum {omegai, di, fi, qi, Ki};



int func (integer ndim, const doublereal *u, const integer *icp, const doublereal *p, 
          integer ijac, doublereal *dfdt, doublereal *dfdu, doublereal *dfdp)
{
    double
        x1 = u[x1i],
        x2 = u[x2i],
        l1 = u[l1i],
        l2 = u[l2i],
        omega= ((double *)p)[omegai],
        d = ((double *)p)[di],
        f = ((double *)p)[fi],
        q = ((double *)p)[qi],
        K = ((double *)p)[Ki];

    double lnormsquare = (l1 * l1 + l2 * l2 );

    dfdt[x1i] = x2;
    dfdt[x2i] = - d * x2 - x1 - x1 * x1 * x1 + f * l1 - q * tanh(K * x2) ;
    dfdt[l1i] = l1 + omega * l2 - l1 * lnormsquare;
    dfdt[l2i] = l2 - omega * l1 - l2 * lnormsquare;

    return 0;
}

int stpnt(integer ndim, const doublereal t, doublereal *u, doublereal *par)
{
    double
        x1 = 0.0,
        x2 = 0.0,
        l1 = sin(2*M_PI*t),
        l2 = cos(2*M_PI*t),
        omega= 0.5,
        d = 0.1,
        f = 0.0,
        q = 0.0,
        K = 100;


    par[omegai] = omega;
    par[di] = d;
    par[fi] = f;
    par[qi] = q;
    par[Ki] = K;
    par[10] = 2*M_PI/omega;

    u[x1i] = x1;
    u[x2i] = x2;
    u[l1i] = l1;
    u[l2i] = l2;

    return 0;
} 

int pvls(integer ndim, const doublereal *u, doublereal *par)
{
    return 0;
} 

int bcnd (integer ndim, const doublereal *par, const integer *icp,
          integer nbc, const doublereal *u0, const doublereal *u1, integer ijac,
          doublereal *fb, doublereal *dbc)
{
  return 0;
}

int icnd (integer ndim, const doublereal *par, const integer *icp,
          integer nint, const doublereal *u, const doublereal *uold,
          const doublereal *udot, const doublereal *upold, integer ijac,
          doublereal *fi, doublereal *dint)
{
  return 0;
}

int fopt (integer ndim, const doublereal *u, const integer *icp,
          const doublereal *par, integer ijac,
          doublereal *fs, doublereal *dfdu, doublereal *dfdp)
{
  return 0;
}

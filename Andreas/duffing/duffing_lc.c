/*
!---------------------------------------------------------------------- 
!---------------------------------------------------------------------- 
!   Duffing Oscillator (Parlitz and Lauterborn 1985)
!---------------------------------------------------------------------- 
!---------------------------------------------------------------------- 
*/

#include "auto_f2c.h"

enum {x1, x2, l1, l2};
enum {omega, d, f, PERIOD=10};



int func (integer ndim, const doublereal *u, const integer *icp, const doublereal *p, 
          integer ijac, doublereal *dfdt, doublereal *dfdu, doublereal *dfdp)
{
    double lnormsquare = ( u[l1] * u[l1] + u[l2] * u[l2] );

    dfdt[x1] = u[x2];
    dfdt[x2] = - p[d] * u[x2] - u[x1] - u[x1] * u[x1] * u[x1] + p[f] * u[l1];
    dfdt[l1] = u[l1] * (1.0 - lnormsquare) + p[omega] * u[l2];
    dfdt[l2] = u[l2] * (1.0 - lnormsquare) - p[omega] * u[l1];

    return 0;
}

int stpnt(integer ndim, const doublereal t, doublereal *u, doublereal *p)
{

    p[omega]  = 1.0;
    p[d]      = 0.2;
    p[f]      = 0.0;
    p[PERIOD] = 2*M_PI/p[omega];

    u[x1] = 0.0;
    u[x2] = 0.0;
    u[l1] = sin(2*M_PI*t);
    u[l2] = cos(2*M_PI*t);

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

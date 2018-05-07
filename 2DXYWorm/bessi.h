/***********************************************************************
*                                                                      *
*     Program to calculate the first kind modified Bessel function     *
*  of integer order N, for any REAL X, using the function BESSI(N,X).  *
*                                                                      *
* -------------------------------------------------------------------- *
*    SAMPLE RUN:                                                       *
*                                                                      *
*    (Calculate Bessel function for N=2, X=0.75).                      *
*                                                                      *
*    Bessel function of order 2 for X =  0.7500:                       *
*                                                                      *
*         Y = 0.073667                                                 *
*                                                                      *
* -------------------------------------------------------------------- *
*    Reference: From Numath Library By Tuan Dang Trong in Fortran 77.  *
*                                                                      *
*                               C++ Release 1.1 By J-P Moreau, Paris.  *
*                                        (www.jpmoreau.fr)             *
*                                                                      *
*    Version 1.1: corected value of P4 in BESSIO (P4=1.2067492 and not *
*                 1.2067429) Aug. 2011.                                *
***********************************************************************/
#ifndef __BESSI__
#define __BESSI__
extern  double BESSI(int, double);
#endif
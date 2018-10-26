/**
* gauss-quad.h
*
* Quadrature points and weights of Gaussian quadrature with legendre poly-
* nomials.
*
* Version 1.0 @ 08/14/2018
*
* Jaerin Lee
* Applied Superconductivity Laboratory
* Dept. of Electrical and Computer Engineering
* Seoul National University
*/

#ifndef __GAUSS_QUAD_H__
#define __GAUSS_QUAD_H__

#define QUAD_ORDER 8

/* 8-point */
#define GL8X0 -.9602898564975363
#define GL8X1 -.7966664774136267
#define GL8X2 -.5255324099163290
#define GL8X3 -.1834346424956498
#define GL8X4 .1834346424956498
#define GL8X5 .5255324099163290
#define GL8X6 .7966664774136267
#define GL8X7 .9602898564975363

#define GL8W0 .1012285362903706
#define GL8W1 .2223810344533744
#define GL8W2 .3137066458778874
#define GL8W3 .3626837833783621
#define GL8W4 .3626837833783621
#define GL8W5 .3137066458778874
#define GL8W6 .2223810344533744
#define GL8W7 .1012285362903706

/* Degree of order of the Gaussian-Legendre quadrature. */
const double x[QUAD_ORDER] = { GL8X0, GL8X1, GL8X2, GL8X3, GL8X4, GL8X5, GL8X6, GL8X7 };
const double w[QUAD_ORDER] = { GL8W0, GL8W1, GL8W2, GL8W3, GL8W4, GL8W5, GL8W6, GL8W7 };

#endif

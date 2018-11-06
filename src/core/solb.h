/**
 * solb.h
 *
 * Version 1.0 @ 08/14/2018
 *
 * Jaerin Lee
 * Applied Superconductivity Laboratory
 * Dept. of Electrical and Computer Engineering
 * Seoul National University
 */

#ifndef __SOLB_H__
#define __SOLB_H__

#define DEBUG

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* For test codes. */
#ifdef DEBUG
	#include <time.h>
#endif

#include "mkl.h"

#include "physics.h"
#include "topology.h"
#include "gauss-quad.h"

/* External dependencies. */
extern const double x[];
extern const double w[];

/* Public interfaces. */
mag_field_2d_t solb_single(const top_solenoid_t *sol, double r, double z);

#endif

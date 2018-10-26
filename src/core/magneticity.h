/**
 * magneticity.h
 *
 * Defines data structures and rules defining physics.
 *
 * Version 1.0 @ 08/14/2018
 *
 * Jaerin Lee
 * Applied Superconductivity Laboratory
 * Dept. of Electrical and Computer Engineering
 * Seoul National University
 */

#ifndef __MAGNETICITY_H__
#define __MAGNETICITY_H__

#include <stdio.h>

typedef struct _mag_field_2d_t
{
	double Br;
	double Bz;

	_mag_field_2d_t()
	{
		Br = 0.0;
		Bz = 0.0;
	}

	_mag_field_2d_t(double Br, double Bz)
	{
		this->Br = Br;
		this->Bz = Bz;
	}

	void print()
	{
		printf("[Field Strength]\nBr : %lf (T)\nBz : %lf (T)\n", Br, Bz);
	}
} mag_field_2d_t;

#endif

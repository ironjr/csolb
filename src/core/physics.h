/**
 * physics.h
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

#ifndef __PHYSICS_H__
#define __PHYSICS_H__

#include <stdio.h>

typedef struct _vec2d_t
{
    double r;
    double z;

	_vec2d_t()
	{
		r = 0.0;
		z = 0.0;
	}

	_vec2d_t(double r, double z)
	{
		this->r = r;
		this->z = z;
	}

    _vec2d_t(const struct _vec2d_t *vec)
    {
        this->r = vec->r;
        this->z = vec->z;
    }
} vec2d_t;

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

    _mag_field_2d_t(const struct _mag_field_2d_t *mag)
    {
        this->Br = mag->Br;
        this->Bz = mag->Bz;
    }

	void print()
	{
		printf("[Field Strength]\nBr : %lf (T)\nBz : %lf (T)\n", Br, Bz);
	}
} mag_field_2d_t;

#endif

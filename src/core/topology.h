/**
 * toopolgy.h
 *
 * Defines various magnet topology as a current system.
 *
 * Version 1.0 @ 08/14/2018
 *
 * Jaerin Lee
 * Applied Superconductivity Laboratory
 * Dept. of Electrical and Computer Engineering
 * Seoul National University
 */

#ifndef __TOPOLOGY_H__
#define __TOPOLOGY_H__

#include <stdio.h>

/* Solenoid */
typedef struct _top_solenoid_t
{
	double a1;
	double a2;
	double b1;
	double b2;
	double j;

	_top_solenoid_t(double a1, double a2, double b1, double b2, double j)
	{
		this->a1 = a1;
		this->a2 = a2;
		this->b1 = b1;
		this->b2 = b2;
		this->j = j;
	}
    
    _top_solenoid_t(const struct _top_solenoid_t *sol)
    {
        if (sol == NULL)
        {
            _top_solenoid_t();
            return;
        }
        this->a1 = sol->a1;
        this->a2 = sol->a2;
		this->b1 = sol->b1;
		this->b2 = sol->b2;
		this->j = sol->j;
	}

    _top_solenoid_t()
    {
		this->a1 = 0;
		this->a2 = 0;
		this->b1 = 0;
		this->b2 = 0;
		this->j = 0;
    }

	void print()
	{
		printf("[Solenoid Dimension]\nA1 : %lf\nA2 : %lf\nB1 : %lf\nB2 : %lf\nJs : %lf\n",
			a1, a2, b1, b2, j);
	}
} top_solenoid_t;

#endif

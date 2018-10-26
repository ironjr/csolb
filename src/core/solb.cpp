/**
 * solb.cpp
 *
 * Alternative version of SolB written in C with LAPACKE. Calculate magnetic
 * field by an arbitrary solenoid
 *
 * Version 1.0 @ 08/14/2018
 *
 * Jaerin Lee
 * Applied Superconductivity Laboratory
 * Dept. of Electrical and Computer Engineering
 * Seoul National University
 */

#include "solb.h"

#define M_PI 3.14159265358979323846
#define NEAR_CENTER_THRESHOLD 1e-6
#define ERROR_REF 1e-10

/* Function premitives. */
mag_field_2d_t *solBInternal(const top_solenoid_t *sol, double r, double z);

mag_field_2d_t *solBSingle(const top_solenoid_t *sol, double r, double z)
{
	/* Handle bad inputs. */
	if (sol == NULL)
	{
		fprintf(stderr, "No solenoid has been specified.");
		return new mag_field_2d_t();
	}

	/* Caching dimensions of the solenoid. */
	double a1 = sol->a1;
	double a2 = sol->a2;
	double j = sol->j;

	/* Handling bad solenoid dimensions. */
	if (a2 < a1 || sol->b2 < sol->b1)
	{
		fprintf(stderr, "Wrong solenoid dimension.");
		return new mag_field_2d_t();
	}

	if (r > a1 && r < a2)
	{
		top_solenoid_t *inner = new top_solenoid_t(*sol);
		top_solenoid_t *outer = new top_solenoid_t(*sol);
		inner->a2 = r;
		outer->a1 = r;

		mag_field_2d_t *resInner = solBInternal(inner, r, z);
		double dtmp0 = 0.5e-7 * j * M_PI;
		double dtmp1 = dtmp0 * (r - a1);
		/* Note 1 */
		double dtmp2 = (r / a1 < NEAR_CENTER_THRESHOLD) ? 0 : 0.5 / r;
		double Br = resInner->Br * dtmp1 * dtmp2;
		double Bz = resInner->Bz * dtmp1;

		mag_field_2d_t *resOuter = solBInternal(outer, r, z);
		dtmp1 = dtmp0 * (a2 - r);
		Br += resOuter->Br * dtmp1 * dtmp2;
		Bz += resOuter->Bz * dtmp1;

		delete resInner;
		delete resOuter;
		delete inner;
		delete outer;

		return new mag_field_2d_t(Br, Bz);
	}
	else
	{
		mag_field_2d_t *res = solBInternal(sol, r, z);
		double Br = res->Br;
		double Bz = res->Bz;
		double dtmp = 0.5e-7 * j * (a2 - a1) * M_PI;

		/* Note 1: Exception for near-center field. Field calculation in this way
		 * yields 0 / 0 division causing severe computational error term. */
		if (r / a1 < NEAR_CENTER_THRESHOLD)
			Br = 0.0;
		else
			Br *= 0.5 * dtmp / r;
		Bz *= dtmp;

		delete res;

		return new mag_field_2d_t(Br, Bz);
	}
}

mag_field_2d_t *solBInternal(const top_solenoid_t *sol, double r, double z)
{
	/* Caching dimensions of the solenoid. */
	double a1 = sol->a1;
	double a2 = sol->a2;

	/* Calculate quadrature points. */
	double a[QUAD_ORDER];
	double dtmp = (a1 + a2) * 0.5;
	cblas_dcopy(QUAD_ORDER, x, 1, a, 1);
	cblas_daxpby(QUAD_ORDER, 1, &dtmp, 0, (a2 - a1) * 0.5, a, 1); // End of using dtmp

	/* Intermediate values. */
	double datmp0[QUAD_ORDER];
	double datmp1[QUAD_ORDER];
	double datmp2[QUAD_ORDER];
	double aprsq[QUAD_ORDER];
	double amrsq[QUAD_ORDER];
	cblas_dcopy(QUAD_ORDER, a, 1, datmp0, 1);
	cblas_daxpy(QUAD_ORDER, 1, &r, 0, datmp0, 1);
	vdSqr(QUAD_ORDER, datmp0, aprsq);
	cblas_dcopy(QUAD_ORDER, a, 1, datmp0, 1);
	cblas_daxpy(QUAD_ORDER, -1, &r, 0, datmp0, 1);
	vdSqr(QUAD_ORDER, datmp0, amrsq);

	/* For a single solenoid, total 2 * Gaussian quadrature points amount of
	 * calculations are needed. */

	/* First iteration. */
	double h = sol->b1;

	/* Preparation of common parameters.
	 * r1sq = (a + r) ** 2 + (z - h) ** 2	| make this
	 * r2sq = (a - r) ** 2 + (z - h) ** 2	| make this
	 * ksq  = 4 * a * r / r1sq				|
	 * kpsq = r2sq / r1sq					| make this
	 * csq  = 4 * a * r / (a + r) ** 2		|
	 * cpsq = (a - r) ** 2 / (a + r) ** 2	| make this
	 */
	double kpsq[QUAD_ORDER];
	double kp[QUAD_ORDER];
	dtmp = z - h;
	dtmp *= dtmp;
	vdLinearFrac(QUAD_ORDER, amrsq, aprsq, 1, dtmp, 1, dtmp, kpsq);
	vdSqrt(QUAD_ORDER, kpsq, kp);
	double cpsq[QUAD_ORDER];
	vdDiv(QUAD_ORDER, amrsq, aprsq, cpsq);
	double r1sq[QUAD_ORDER];
	cblas_dcopy(QUAD_ORDER, aprsq, 1, r1sq, 1);
	double r1[QUAD_ORDER];
	cblas_daxpy(QUAD_ORDER, 1, &dtmp, 0, r1sq, 1);
	vdSqrt(QUAD_ORDER, r1sq, r1);

	/* Calculation of initial values of iterative methods.
	 * alpha0	= 1							|
	 * beta0	= kp						|
	 * S0		= ksq						| not needed
	 * SG		= 0							|
	 * delta0	= cpsq / kp					|
	 * epsilon0 = csq / cpsq				|
	 * zeta0    = 0							|
	 */
	double *beta0a = kp;
	double delta0a[QUAD_ORDER];
	double epsilon0a[QUAD_ORDER];
	vdLinearFrac(QUAD_ORDER, cpsq, kp, 1, 0, 1, 0, delta0a);
	vdLinearFrac(QUAD_ORDER, cpsq, cpsq, -1, 1, 1, 0, epsilon0a);

	double alphaInf[QUAD_ORDER];
	double zetaInf[QUAD_ORDER];
	double SGInf[QUAD_ORDER];
	for (int i = 0; i < QUAD_ORDER; ++i)
	{
		/* Initialization of iterative method. */
		double alpha = 1;
		double beta = beta0a[i];
		double delta = delta0a[i];
		double epsilon = epsilon0a[i];
		double zeta = 0;
		double SG = 0;

		/* Main iterative loop.
		 * alphaNext = AM(alpha, beta)
		 * betaNext = GM(alpha, beta)
		 * deltaNext = (betaNext / (4 * alphaNext)) * (2 + delta + 1 / delta)
		 * epsilonNext = delta * (epsilon + zeta) / (1 + delta)
		 * zetaNext = AM(epsilon, zeta)
		 * Garrett's Sum = SIGMA_i_0_to_inf(2 ** (i - 1) * (alpha - beta) ** 2)
		 */
		int j = 0;
		double error = 1;
		while (1)
		{
			/* Update Garrett's sum. */
			double temp = (alpha - beta);
			temp *= temp;
			temp *= double(1 << j);
			SG += temp;
			if (error < ERROR_REF && abs(1 - delta) < ERROR_REF) break;

			/* Evaluate the error. */
			error = temp;

			/* Update iterative variables. */
			temp = sqrt(alpha * beta);
			alpha = (alpha + beta) * 0.5;
			beta = temp;
			temp = (epsilon + zeta) / 2;
			epsilon = (delta * epsilon + zeta) / (1 + delta);
			zeta = temp;
			delta = (2 + delta + 1 / delta) * beta / (4 * alpha);

			++j;
		}

		alphaInf[i] = alpha;
		zetaInf[i] = zeta;
		SGInf[i] = SG * 0.5;
	}

	/* Calculation of B. */
	vdMul(QUAD_ORDER, w, r1, datmp0);
	vdDiv(QUAD_ORDER, SGInf, alphaInf, datmp1);
	double BrDiff = -cblas_ddot(QUAD_ORDER, datmp0, 1, datmp1, 1);
	
	vdAdd(QUAD_ORDER, a, a, datmp0);
	vdLinearFrac(QUAD_ORDER, a, datmp0, 1, -r, 0, 1, datmp1);
	vdMul(QUAD_ORDER, datmp1, zetaInf, datmp2);
	vdAdd(QUAD_ORDER, datmp0, datmp2, datmp1);
	vdLinearFrac(QUAD_ORDER, a, datmp1, 1, r, 0, 1, datmp0);
	vdMul(QUAD_ORDER, datmp0, r1, datmp2);
	vdMul(QUAD_ORDER, datmp2, alphaInf, datmp0);
	vdDiv(QUAD_ORDER, datmp1, datmp0, datmp2);
	double BzDiff = cblas_ddot(QUAD_ORDER, w, 1, datmp2, 1) * (z - h);

	/* Second iteration. */
	h = sol->b2; // DIFF

	/* Preparation of common parameters.
	 * r1sq = (a + r) ** 2 + (z - h) ** 2	| make this
	 * r2sq = (a - r) ** 2 + (z - h) ** 2	| make this
	 * ksq  = 4 * a * r / r1sq				|
	 * kpsq = r2sq / r1sq					| make this
	 * csq  = 4 * a * r / (a + r) ** 2		|
	 * cpsq = (a - r) ** 2 / (a + r) ** 2	| make this
	 */
	dtmp = z - h;
	dtmp *= dtmp;
	vdLinearFrac(QUAD_ORDER, amrsq, aprsq, 1, dtmp, 1, dtmp, kpsq);
	vdSqrt(QUAD_ORDER, kpsq, kp);
	vdDiv(QUAD_ORDER, amrsq, aprsq, cpsq);
	cblas_dcopy(QUAD_ORDER, aprsq, 1, r1sq, 1);
	cblas_daxpy(QUAD_ORDER, 1, &dtmp, 0, r1sq, 1);
	vdSqrt(QUAD_ORDER, r1sq, r1);

	/* Calculation of initial values of iterative methods.
	 * alpha0	= 1							|
	 * beta0	= kp						|
	 * S0		= ksq						| not needed
	 * SG		= 0							|
	 * delta0	= cpsq / kp					|
	 * epsilon0 = csq / cpsq				|
	 * zeta0    = 0							|
	 */
	beta0a = kp;
	vdLinearFrac(QUAD_ORDER, cpsq, kp, 1, 0, 1, 0, delta0a);
	vdLinearFrac(QUAD_ORDER, cpsq, cpsq, -1, 1, 1, 0, epsilon0a);

	for (int i = 0; i < QUAD_ORDER; ++i)
	{
		/* Initialization of iterative method. */
		double alpha = 1;
		double beta = beta0a[i];
		double delta = delta0a[i];
		double epsilon = epsilon0a[i];
		double zeta = 0;
		double SG = 0;

		/* Main iterative loop.
		 * alphaNext = AM(alpha, beta)
		 * betaNext = GM(alpha, beta)
		 * deltaNext = (betaNext / (4 * alphaNext)) * (2 + delta + 1 / delta)
		 * epsilonNext = delta * (epsilon + zeta) / (1 + delta)
		 * zetaNext = AM(epsilon, zeta)
		 * Garrett's Sum = SIGMA_i_0_to_inf(2 ** (i - 1) * (alpha - beta) ** 2)
		 */
		int j = 0;
		double error = 1;
		while (1)
		{
			/* Update Garrett's sum. */
			double temp = (alpha - beta);
			temp *= temp;
			temp *= double(1 << j);
			SG += temp;
			if (error < ERROR_REF && abs(1 - delta) < ERROR_REF) break;

			/* Evaluate the error. */
			error = temp;

			/* Update iterative variables. */
			temp = sqrt(alpha * beta);
			alpha = (alpha + beta) * 0.5;
			beta = temp;
			temp = (epsilon + zeta) / 2;
			epsilon = (delta * epsilon + zeta) / (1 + delta);
			zeta = temp;
			delta = (2 + delta + 1 / delta) * beta / (4 * alpha);

			++j;
		}

		alphaInf[i] = alpha;
		zetaInf[i] = zeta;
		SGInf[i] = SG * 0.5;
	}

	/* Calculation of B. */
	vdMul(QUAD_ORDER, w, r1, datmp0);
	vdDiv(QUAD_ORDER, SGInf, alphaInf, datmp1);
	BrDiff += cblas_ddot(QUAD_ORDER, datmp0, 1, datmp1, 1); // DIFF
	
	vdAdd(QUAD_ORDER, a, a, datmp0);
	vdLinearFrac(QUAD_ORDER, a, datmp0, 1, -r, 0, 1, datmp1);
	vdMul(QUAD_ORDER, datmp1, zetaInf, datmp2);
	vdAdd(QUAD_ORDER, datmp0, datmp2, datmp1);
	vdLinearFrac(QUAD_ORDER, a, datmp1, 1, r, 0, 1, datmp0);
	vdMul(QUAD_ORDER, datmp0, r1, datmp2);
	vdMul(QUAD_ORDER, datmp2, alphaInf, datmp0);
	vdDiv(QUAD_ORDER, datmp1, datmp0, datmp2);
	BzDiff -= cblas_ddot(QUAD_ORDER, w, 1, datmp2, 1) * (z - h); // DIFF

	mag_field_2d_t *ans = new mag_field_2d_t(BrDiff, BzDiff);
	return ans;
}

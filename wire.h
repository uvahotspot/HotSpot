#ifndef __WIRE_H
#define __WIRE_H

#include "flp.h"
#include "math.h"

/* 
 * a simple wire-delay model based on 
 * 
 * 1) Equations 1-3 of the following paper: 
 * R. H. J. M. Otten and R. K. Brayton, "Planning for 
 * performance," in DAC '98: Proceedings of the 35th annual 
 * conference on Design automation, pp. 122--127, 1998.
 *
 * and
 *
 * 2) Equations 30-33 of the following paper: 
 * K. Banerjee and A. Mehrotra, "Global (interconnect) 
 * warming," IEEE Circuits and Devices Magazine, vol. 17, 
 * pp. 16--32, September 2001. The lcrit and sopt values
 * for the global metal layers defined below are taken from 
 * Table 4 of this paper.
 *
 * The resistance and capacitance values of wires per unit 
 * length defined below (WIRE_R*, WIRE_C*) are taken from 
 * Table 2 of the following paper:
 *
 * V. Agarwal, S. W. Keckler, and D. Burger, "The effect 
 * of technology scaling on microarchitectural structures,"
 * Tech. Rep. TR-00-02, University of Texas at Austin Computer 
 * Sciences, May 2001.
 */

/* we are modeling a 130 nm wire	*/
#define TECHNODE130

/* metal layer	*/
#define WIRE_GLOBAL			0
#define WIRE_INTER			1

/* 
 * constants dependent on the swing assumptions of the
 * repeaters inserted. usually, a 50% swing is assumed.
 * in such a case, a = 0.4, b = 0.7
 */
#define	WIRE_A				0.4
#define	WIRE_B				0.7

/* 
 * *_G are values for global metal while *_I are
 * for intermediate metal layers respectively. 
 * The units used are listed below:
 * WIRE_R_* - mohm/u
 * WIRE_C_* - fF/u
 * WIRE_LCRIT_* - mm
 * WIRE_SOPT_* - dimensionless
 */
#if	defined(TECHNODE180)
	#define	WIRE_R_G		36.0
	#define	WIRE_C_G		0.350
	#define	WIRE_R_I		107.0
	#define	WIRE_C_I		0.333
	#define	WIRE_LCRIT_G	3.0
	#define	WIRE_SOPT_G		179.0
#elif defined(TECHNODE130)
	#define	WIRE_R_G		61.0
	#define	WIRE_C_G		0.359
	#define	WIRE_R_I		188.0
	#define	WIRE_C_I		0.336
	#define	WIRE_LCRIT_G	2.4
	#define	WIRE_SOPT_G		146.0
#elif defined(TECHNODE100)
	#define	WIRE_R_G		103.0
	#define	WIRE_C_G		0.361
	#define	WIRE_R_I		316.0
	#define	WIRE_C_I		0.332
	#define	WIRE_LCRIT_G	2.12
	#define	WIRE_SOPT_G		96.0
#elif defined(TECHNODE70)
	#define	WIRE_R_G		164.0
	#define	WIRE_C_G		0.360
	#define	WIRE_R_I		500.0
	#define	WIRE_C_I		0.331
	#define	WIRE_LCRIT_G	1.2
	#define	WIRE_SOPT_G		82.0
#elif defined(TECHNODE50)
	#define	WIRE_R_G		321.0
	#define	WIRE_C_G		0.358
	#define	WIRE_R_I		1020.0
	#define	WIRE_C_I		0.341
	#define	WIRE_LCRIT_G	0.99
	#define	WIRE_SOPT_G		48.0
#endif	

/* 
 * lcrit and sopt values for intermediate layer
 * are derived from the values for global layer.
 * lcrit * sqrt(wire_r * wire_c) is a constant
 * for a given layer. So is the expression
 * sopt * sqrt(wire_r / wire_c) (Equation 2 and
 * Theorem 2 from Brayton et. al). Using this
 * for the global layer, one can find the constants
 * for the intermediate layer and combining this
 * with the wire_r and wire_c values for the
 * intermediate layer, its lcrit and sopt values
 * can be calculated.
 */
#define	WIRE_LCRIT_I		(WIRE_LCRIT_G * \
							 sqrt((WIRE_R_G / WIRE_R_I) * (WIRE_C_G / WIRE_C_I)))
#define	WIRE_SOPT_I			(WIRE_SOPT_G * \
							 sqrt((WIRE_R_G / WIRE_R_I) / (WIRE_C_G / WIRE_C_I)))

/* 
 * delay (in secs) per meter of a wire using 
 * Equation 2 and Theorem 2 of Brayton et al. 
 * the following expression is derived by 
 * setting cp = co in the delay equation and
 * substituting for sqrt(ro*co) in the same
 * from Equation 2. The 1.0e-9 is added to take
 * care of the units
 */
#define WIRE_DELAY_G		(WIRE_LCRIT_G * WIRE_R_G * WIRE_C_G * 1.0e-9 * \
							(2.0 * WIRE_A + sqrt (2.0 * WIRE_A * WIRE_B)))
#define WIRE_DELAY_I		(WIRE_LCRIT_I * WIRE_R_I * WIRE_C_I * 1.0e-9 * \
							(2.0 * WIRE_A + sqrt (2.0 * WIRE_A * WIRE_B)))

/* function for wire-length to wire-delay conversion	*/
double wire_length2delay(double length, int layer);

#endif

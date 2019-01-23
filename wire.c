#include "wire.h"
#include "util.h"

/* 
 * function for wire-length to wire-delay conversion.
 * WIRE_DELAY_G and WIRE_DELAY_I are delay values
 * per unit length. see wire.h for the details of the
 * wire model and the computation of these constants
 */
double wire_length2delay(double length, int layer)
{
	if (layer == WIRE_GLOBAL)
		return WIRE_DELAY_G * length;
	else if (layer == WIRE_INTER)
		return WIRE_DELAY_I * length;
	else
		fatal("unknown metal layer\n");
	return 0.0;	
}

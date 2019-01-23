#ifndef __PACKAGE_H_
#define __PACKAGE_H_

#include "flp.h"
#include "util.h"
#include "temperature.h"

/* santity check convection resistance boundary values.
	 if outside this range, most likely some settings of 
	 the heat sink or the fan are unreasonable.
*/
#define R_CONVEC_HIGH	50.0 /* upper-bound value for natural convection without heat sink */
#define R_CONVEC_LOW	0.01 /* lower-bound value for best cooling solution, e.g. liquid cooling */

/* fluid (air) properties	*/
#define AIR_DSTY		1.059	/* density */
#define AIR_SPECHT		1007 /* specific heat */
#define AIR_COND	0.028	/* thermal conductivity */
#define AIR_DIFF	2.6e-5	/* thermal diffusivity */
#define AIR_DYNVISC		2.0e-5	/* dynamic viscosity */
#define AIR_KINVISC		1.9e-5	/* kinetic viscosity */

#define PRANTDL_NUM	0.73	/* Prantdl Number of air */
#define REY_THRESHOLD	3500 /* laminar/turbulent reynolds threshold */ 

/* parameters specific to natural convection and radiation */
#define VOL_EXP_COEFF	3.0e-3 /* volume expansion coefficient */
#define GRAVITY	9.8	/* gravity acceleration rate */
#define STEFAN	5.67e-8 /* Stefan-Boltzmann constant */
#define EMISSIVITY	0.95 /* emissivity, typical value for most heatsink materials */
#define SMALL_FOR_CONVEC 0.01 /* initial small temperature diff b/w heatsink and ambient, to avoid sigular matrices */
#define NATURAL_CONVEC_TOL 0.01 /* r_convec convergence criterion for natural convection calculations */
#define MAX_SINK_TEMP 1000	/* max avg sink temperature during iteration for natural convection. if exceeded, report thermal runaway*/

/* flow resistance constants */
#define KC	0.42
#define KE	1.0
#define PI	3.1416

/* fan constants */
#define FAN_ALPHA	1 /* 1 for laminar flow, >=2 for turbulent flow */
#define FAN_BETA 4 /* 4 for laminar flow, 4.6 for turbulent flow */
#define FAN_POWER_COEFF 2.8 /* fan_power=b*omega^2.8 */
#define RPM_TO_RAD 0.105	/* convert rpm to radians per sec */
#define RAD_TO_RPM	9.55 /* convert radians per sec to rpm */

/* package model configuration	*/
typedef struct package_config_t_st
{
	/* natural convection or not - 0: forced convection, 1: natural convection*/
	int natural_convec;
	
	/* airflow type - 0: lateral airflow from sink side, 1: impinging airflow	from sink top*/
	int flow_type;
	
	/* heatsink type - 0: fin-channel sink, 1: pin-fin sink */
	int sink_type;
	
	/* sink specs */
		/* sink base size is defined in thermal_config */
		/* 1) fin-channel sink */
	double fin_height;
	double fin_width;
	double channel_width;
		/* 2) pin-fin sink */
	double pin_height;
	double pin_diam;
	double pin_dist; /* distance between pins */ 
	
	/* fan specs */
	double fan_radius;
	double motor_radius;
	int rpm;
}package_config_t;

/* defaults	*/
package_config_t default_package_config(void);
/* 
 * parse a table of name-value string pairs and add the configuration
 * parameters to 'package_config'
 */
void package_config_add_from_strs(package_config_t *package_config, str_pair *package_table, int size);
/* 
 * convert package_config into a table of name-value pairs. returns the no.
 * of parameters converted
 */
int package_config_to_strs(package_config_t *package_config, str_pair *package_table, int max_entries);

/* airflow parameters	*/
typedef struct convection_t_st
{
	double n_fin; /* number of fins */
	double sur_area_fin; /* total surface area */ 

	double n_pin; /* number of pins */
	double sur_area_pin; /* total surface area */ 

	/* flow charateristics */
	double reynolds;
	double nusselt;
	double h_coeff; /* heat transfer coefficient */
	double v; /* air velocity */
	double r_th; /* lumped convection thermal resistance */
}convection_t;

/*  airflow parameter routines	*/
void calculate_flow(convection_t *p, package_config_t *package_config, thermal_config_t *thermal_config);

/* debug print	*/
void debug_print_convection(convection_t *p);

/* parse thermal config table and update r_conv if needed */
int package_model(thermal_config_t *thermal_config, str_pair *table, int size, double sink_temp);

#endif

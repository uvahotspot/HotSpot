#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef _MSC_VER
#define strcasecmp    _stricmp
#define strncasecmp   _strnicmp
#else
#include <strings.h>
#endif
#include <math.h>

#include "package.h"
#include "temperature.h"
#include "flp.h"
#include "util.h"

/* default package configuration parameters	*/
package_config_t default_package_config(void)
{
	package_config_t config;
	/* 0: forced convection, 1: natural convection */
	config.natural_convec = 0;
	
	/* airflow type - 0: lateral airflow from sink side, 1: impinging airflow	from sink top*/
	config.flow_type = 0;
	
	/* heatsink type - 0: fin-channel sink, 1: pin-fin sink */
	config.sink_type = 0;
	
	/* default sink specs */
		/* sink base size is defined in thermal_config
		 * 1) fin-channel sink 
		 */
	config.fin_height = 0.03;
	config.fin_width = 0.001;
	config.channel_width = 0.002; 
		/* 2) pin-fin sink */
	config.pin_height = 0.02;
	config.pin_diam = 0.002;
	config.pin_dist = 0.005;
	
	/* fan specs */
	config.fan_radius = 0.03;
	config.motor_radius = 0.01;
	config.rpm = 1000;
	
	return config;
}

/* 
 * parse a table of name-value string pairs and add the configuration
 * parameters to package 'config'
 */
void package_config_add_from_strs(package_config_t *config, str_pair *table, int size)
{
	int idx;
	if ((idx = get_str_index(table, size, "natural_convec")) >= 0)
		if(sscanf(table[idx].value, "%d", &config->natural_convec) != 1)
			fatal("invalid format for heatsink configuration parameter natural_convec\n");
	if ((idx = get_str_index(table, size, "flow_type")) >= 0)
		if(sscanf(table[idx].value, "%d", &config->flow_type) != 1)
			fatal("invalid format for heatsink configuration parameter flow_type\n");
	if ((idx = get_str_index(table, size, "sink_type")) >= 0)
		if(sscanf(table[idx].value, "%d", &config->sink_type) != 1)
			fatal("invalid format for heatsink configuration parameter sink_type\n");
	if ((idx = get_str_index(table, size, "fin_height")) >= 0)
		if(sscanf(table[idx].value, "%lf", &config->fin_height) != 1)
			fatal("invalid format for heatsink configuration parameter fin_height\n");
	if ((idx = get_str_index(table, size, "fin_width")) >= 0)
		if(sscanf(table[idx].value, "%lf", &config->fin_width) != 1)
			fatal("invalid format for heatsink configuration parameter fin_width\n");
	if ((idx = get_str_index(table, size, "channel_width")) >= 0)
		if(sscanf(table[idx].value, "%lf", &config->channel_width) != 1)
			fatal("invalid format for heatsink configuration parameter channel_width\n");
	if ((idx = get_str_index(table, size, "pin_height")) >= 0)
		if(sscanf(table[idx].value, "%lf", &config->pin_height) != 1)
			fatal("invalid format for heatsink configuration parameter pin_height\n");
	if ((idx = get_str_index(table, size, "pin_diam")) >= 0)
		if(sscanf(table[idx].value, "%lf", &config->pin_diam) != 1)
			fatal("invalid format for heatsink configuration parameter pin_diam\n");
	if ((idx = get_str_index(table, size, "pin_dist")) >= 0)
		if(sscanf(table[idx].value, "%lf", &config->pin_dist) != 1)
			fatal("invalid format for heatsink configuration parameter pin_dist\n");
	if ((idx = get_str_index(table, size, "fan_radius")) >= 0)
		if(sscanf(table[idx].value, "%lf", &config->fan_radius) != 1)
			fatal("invalid format for fan configuration parameter fan_radius\n");
	if ((idx = get_str_index(table, size, "motor_radius")) >= 0)
		if(sscanf(table[idx].value, "%lf", &config->motor_radius) != 1)
			fatal("invalid format for fan configuration parameter motor_radius\n");
	if ((idx = get_str_index(table, size, "rpm")) >= 0)
		if(sscanf(table[idx].value, "%d", &config->rpm) != 1)
			fatal("invalid format for configuration  parameter rpm\n");
	
	if ((config->fin_height <= 0) || (config->fin_width <= 0) || (config->channel_width <= 0) || 
		(config->pin_height <= 0) || (config->pin_diam <= 0) || (config->pin_dist <= 0) || 
		(config->fan_radius <= 0) || (config->motor_radius <= 0) || (config->rpm <= 0))
		fatal("heatsink/fan dimensions and fan speed should be greater than zero\n");
	if ((config->natural_convec != 0) && (config->natural_convec != 1))
		fatal("invalid convection mode\n");
	if ((config->flow_type != 0) && (config->flow_type != 1))
		fatal("invalid air flow type\n");
	if ((config->sink_type != 0) && (config->sink_type != 1))
		fatal("invalid heatsink type\n");
}

/* 
 * convert config into a table of name-value pairs. returns the no.
 * of parameters converted
 */
int package_config_to_strs(package_config_t *package_config, str_pair *package_table, int max_entries)
{
	if (max_entries < 12)
		fatal("not enough entries in table\n");
	
	sprintf(package_table[0].name, "natural_convec");
	sprintf(package_table[1].name, "flow_type");
	sprintf(package_table[2].name, "sink_type");
	sprintf(package_table[3].name, "fin_height");
	sprintf(package_table[4].name, "fin_width");
	sprintf(package_table[5].name, "channel_width");
	sprintf(package_table[6].name, "pin_height");
	sprintf(package_table[7].name, "pin_diam");
	sprintf(package_table[8].name, "pin_dist");
	sprintf(package_table[9].name, "fan_radius");
	sprintf(package_table[10].name, "motor_radius");
	sprintf(package_table[11].name, "rpm");

	sprintf(package_table[0].value, "%d", package_config->natural_convec);
	sprintf(package_table[1].value, "%d", package_config->flow_type);
	sprintf(package_table[2].value, "%d", package_config->sink_type);
	sprintf(package_table[3].value, "%lg", package_config->fin_height);
	sprintf(package_table[4].value, "%lg", package_config->fin_width);
	sprintf(package_table[5].value, "%lg", package_config->channel_width);
	sprintf(package_table[6].value, "%lg", package_config->pin_height);
	sprintf(package_table[7].value, "%lg", package_config->pin_diam);
	sprintf(package_table[8].value, "%lg", package_config->pin_dist);
	sprintf(package_table[9].value, "%lg", package_config->fan_radius);
	sprintf(package_table[10].value, "%lg", package_config->motor_radius);
	sprintf(package_table[11].value, "%d", package_config->rpm);

	return 12;
}

/* calculate forced air flow and package thermal parameters	
 * references:															
 * 1-1	lateral flow, fin-channel heat sinks:		
 *			-laminar flow: P. terrtstra et al. "Analytical Forced Convection Modeling of Plate-Fin Heat Sinks". IEEE SEMI-THERM, pp 34-41, 1999
 *			-turbulent flow: Y. A. Cengel. "Heat and Mass Transfer: A Practical Approach", Mcgraw-Hill Inc. New York 2007 
 * 1-2	lateral flow, pin-fin heat sinks:					
 * 			R. Ribando et al. "Estimating the Convection Coefficient for Flow Through Banks of Fins on a Heat Sink". U.Va. MAE314 Course Notes, Spring 2007
 * 2-1	impinging flow for both fin-channel and pin-fin heat sinks:
 *			H. A. El-Sheikh et al. "Heat Transfer from Pin-Fin Heat Sinks under Multiple Impinging Jets", IEEE Trans. on Adv. Packaging, 23(1):113-120, Feb. 2000
 * 3-1	fan model:
 *			F. P. Bleier. "Fan Handbook: Selection, Application and Design". McGraw-Hill Inc. New York, 1998
 * 3-2	motor model of the fan:
 * 			Y. Zhang et al. "SODA: Sensitivity-Based Optimization of Disk Architecture". IEEE/ACM DAC, pp. 865-870, 2007 
 */
 
void calculate_flow(convection_t *p, package_config_t *config, thermal_config_t *thermal_config)
{
	/* local variables */
	double n_fin, n_pin;
	double sur_area_fin, sur_area_pin;
	double reynolds, nusselt, h_coeff, v, r_th;
	double dh, a_fan, dr, r_approx, vol_v, rey_star, m, eta, f, c1, c2, a_hs;
	double t1, t2; /* temporary variables for long formulas*/
	
	double s_sink = thermal_config->s_sink;
	double k_sink = thermal_config->k_sink;
	
	int flow_type = config->flow_type;
	int sink_type = config->sink_type;
	double fin_height = config->fin_height;
	double fin_width = config->fin_width;
	double channel_width = config->channel_width;
	double pin_height = config->pin_height;
	double pin_diam = config->pin_diam;
	double pin_dist = config->pin_dist;
	double fan_radius = config->fan_radius;
	double motor_radius = config->motor_radius;
	double rpm = config->rpm;
	
	double temp_val = (s_sink-pin_diam)/(pin_diam+pin_dist);
	
	n_fin = ceil((s_sink-fin_width)/(fin_width+channel_width)-0.5);
	sur_area_fin = s_sink*(s_sink+2.0*n_fin*fin_height);
	n_pin = ceil(temp_val*temp_val-0.5);
	sur_area_pin = s_sink*s_sink+PI*pin_diam*pin_height*n_pin;
		
	/* calculate volumetric air speed out of the fan */
	dr = sqrt(fan_radius * fan_radius - motor_radius * motor_radius);
	r_approx = motor_radius + dr; /* approximated average fan radius */
	a_fan = PI * dr * dr; /* total fan blade area */
	
	/* refer to DAC'07 SODA paper by Zhang, Gurumurthi and Stan 
	 * crudely approximating an IC fan motor with a hard drive spindle motor model
	 * the principle is that the dragging momentum on the blades from the air equals to the torque of the motor at steady state
	 * so, torque=b*(omega^alpha)=drag_force*radius (1)
	 * where b=0.5*pi*air_density*C_d*(radius^beta)
	 * drag coeff C_d=drag_force/(0.5*air_density*air_velocity^2*total_blade_area)
	 * manipulate both side of (1), derive volmetric velocity from the fan as... 
	 */
	vol_v = a_fan * sqrt(0.25 * PI * pow(r_approx,FAN_BETA-1) * pow(rpm * RPM_TO_RAD,FAN_ALPHA) / a_fan);
	
	/* calculate the actual air velocity through heatsink: vol_velocity/area_duct_sink */
	if (flow_type==0) { /* lateral airflow */
		dh = 2.0*channel_width*s_sink/(channel_width+s_sink); /* hydraulic diameter */
		if (sink_type==0) { /* fin-channel sink */
			v = vol_v / ((n_fin-1)*channel_width*fin_height);
		}
		else { /* pin-fin sink */
			v = vol_v / ((sqrt(n_pin)-1)*pin_dist*pin_height);
		}
	}
	else { /* impinging flow */
		dh = 2.0 * s_sink / sqrt(PI); /* equivalent air nozzle diameter */ 
		v = vol_v / (s_sink*s_sink-n_pin*PI*(pin_diam*0.5)*(pin_diam*0.5));
	}
	
	/* Reynolds number */
	reynolds = AIR_DSTY * v * dh / AIR_DYNVISC;
	
	/* calculate nusselt number, heat transfer coeff
	 * and equivalent overall lumped convection resistance 
	 */
	if (flow_type==0) { /* lateral airflow */
		if (sink_type==0) { /* fin-channel sink */
			if (reynolds <= REY_THRESHOLD) { /* laminar flow */
				rey_star = AIR_DSTY*v*channel_width*channel_width/(AIR_DYNVISC*s_sink);
				t1 = pow(rey_star*PRANTDL_NUM*0.5,-3);
				t2 = pow(0.664*sqrt(rey_star)*pow(PRANTDL_NUM,0.33)*sqrt(1+3.65/sqrt(rey_star)),-3);
				nusselt = pow(t1+t2,-0.33); /* nusselt number */
				h_coeff = nusselt*AIR_COND/channel_width; /* heat transfer coefficient */
				m = sqrt(2*h_coeff/(k_sink*fin_width)); /* fin parameter */
				eta = tanh(m*s_sink)/(m*s_sink); /* fin efficiency */
				r_th = 1.0/(h_coeff*(channel_width+2*eta*s_sink)*n_fin*s_sink);
			}
			else { /* turbulent flow */
				f = 1.0/pow((0.79*log(reynolds)-1.64),2.0);
				nusselt = 0.125*f*reynolds*pow(PRANTDL_NUM,0.33);
				h_coeff = nusselt*AIR_COND/dh;
				m = sqrt(2*h_coeff/(k_sink*fin_width)); /* fin parameter */
				eta = tanh(m*s_sink)/(m*s_sink); /* fin efficiency */
				r_th = 1.0/(h_coeff*(channel_width+2*eta*s_sink)*n_fin*s_sink);
			}
		}				
		else { /* pin-fin sink */
			c1 = 0.6;
			c2 = 0.8+0.8333*sqrt(sqrt(n_pin)-1)/pow((1+0.03*pow((sqrt(n_pin)-1),2.0)),0.25);
		  nusselt = 0.4+0.9*(1.25*c2*c1*sqrt(reynolds)+0.001*reynolds)*pow(PRANTDL_NUM,0.33);
			h_coeff = nusselt*AIR_COND/pin_diam;
			r_th = 1.0/(h_coeff*sur_area_pin);
		}
	}
	else { /* impinging flow */
		if (sink_type==0) { /* fin-channel sink */
			a_hs = sur_area_fin;
		}
		else { /* pin-fin sink */
			a_hs = sur_area_pin;
		}
		dh = 2.0 * s_sink / sqrt(PI);
		nusselt = 1.92*pow(reynolds,0.716)*pow(PRANTDL_NUM,0.4)*pow(a_hs/(s_sink*s_sink),-0.698);
		h_coeff = nusselt*AIR_COND/dh;
		r_th = 1.0/(h_coeff*a_hs);
	}
	p->n_fin = n_fin;
	p->sur_area_fin = sur_area_fin;
	p->n_pin = n_pin;
	p->sur_area_pin = sur_area_pin;
	p->reynolds = reynolds;
	p->nusselt = nusselt;
	p->h_coeff = h_coeff;
	p->v = v;
	p->r_th = r_th;
}

/* calculate convection parameters for natural convection.
 * reference: Y. A. Cengel. "Heat and Mass Transfer: A Practical Approach", Mcgraw-Hill Inc. New York 2007 
 */
void calc_natural_convec(convection_t *p, package_config_t *config, thermal_config_t *thermal_config, double sink_temp) 
{
	/* local variables */
	double rayleigh;
	double n_fin, n_pin;
	double w;
	double sur_area, sur_area_fin, sur_area_pin;
	double nusselt, h_coeff, r_th, r_th_rad;
	
	double s_sink = thermal_config->s_sink;
	double ambient = thermal_config->ambient;
	
	int sink_type = config->sink_type;
	double fin_height = config->fin_height;
	double fin_width = config->fin_width;
	double channel_width = config->channel_width;
	double pin_height = config->pin_height;
	double pin_diam = config->pin_diam;
	double pin_dist = config->pin_dist;
	
	double temp_val = (s_sink-pin_diam)/(pin_diam+pin_dist);
	
	n_fin = ceil((s_sink-fin_width)/(fin_width+channel_width)-0.5);
	sur_area_fin = s_sink*(s_sink+2.0*n_fin*fin_height);
	n_pin = ceil(temp_val*temp_val-0.5);
	sur_area_pin = s_sink*s_sink+PI*pin_diam*pin_height*n_pin;
	
	if (sink_type==0) { /* fin-channel sink */
		sur_area = sur_area_fin;
		w = channel_width;
	} else { /* pin-fin heatsink */
		sur_area = sur_area_pin;
		w = pin_dist;
	}
	
	/* CAUTION: equations are derived for fin-channel heatsink, not validated for pin-fin heatsink */
	rayleigh = GRAVITY*VOL_EXP_COEFF*(sink_temp-ambient)*pow(w,3.0)*PRANTDL_NUM/(AIR_KINVISC*AIR_KINVISC);
	nusselt = pow((576/((rayleigh*w/s_sink)*(rayleigh*w/s_sink))+2.873/sqrt(rayleigh*w/s_sink)),-0.5);
	h_coeff = nusselt*AIR_COND/w;
	r_th = 1.0/(h_coeff*sur_area);
	
	/* thermal radiation*/
	r_th_rad = (sink_temp-ambient)/(EMISSIVITY*STEFAN*(pow(sink_temp,4.0)-pow(ambient,4.0))*sur_area);
	
	/* overall thermal resistance = natural convection in parallel with thermal radiation */
	r_th = r_th*r_th_rad/(r_th+r_th_rad);
	
	p->n_fin = n_fin;
	p->sur_area_fin = sur_area_fin;
	p->n_pin = n_pin;
	p->sur_area_pin = sur_area_pin;
	p->nusselt = nusselt;
	p->h_coeff = h_coeff;
	p->r_th = r_th;
}

/* debug print	*/
void debug_print_convection(convection_t *p)
{
	fprintf(stdout, "printing airflow information...\n");
	fprintf(stdout, "n_fin: %f\n", p->n_fin);
	fprintf(stdout, "sur_area_fin: %f\n", p->sur_area_fin);
	fprintf(stdout, "n_pin: %f\n", p->n_pin);
	fprintf(stdout, "sur_area_pin: %f\n", p->sur_area_pin);
	fprintf(stdout, "reynolds: %f\n", p->reynolds);
	fprintf(stdout, "nusselt: %f\n", p->nusselt);
	fprintf(stdout, "h_coeff: %f\n", p->h_coeff);
	fprintf(stdout, "v: %f\n", p->v);
	fprintf(stdout, "r_th: %f\n", p->r_th);
}

/* initialize and calculate package parameters and update r_convec */
int package_model(thermal_config_t *thermal_config, str_pair *table, int size, double sink_temp)
{
	int idx;
	int natural_convec;
	
	str_pair package_table[MAX_ENTRIES];
	int package_size;
	
	/* package config parameters */
	package_config_t package_config;
	convection_t p;
	
	package_size = 0;
	
	/* get defaults */
	package_config = default_package_config();
	/* parse the package config file name */
	if ((idx = get_str_index(table, size, "package_config_file")) >= 0) {
		if(sscanf(table[idx].value, "%s", thermal_config->package_config_file) != 1)
			fatal("invalid format for configuration  parameter package_config_file\n");
	}
	/* read package config file	*/
	if (strcmp(thermal_config->package_config_file, NULLFILE)) 
		package_size += read_str_pairs(&package_table[package_size], MAX_ENTRIES, thermal_config->package_config_file);
			
	/* modify according to package config file	*/
	package_config_add_from_strs(&package_config, package_table, package_size);
	
	natural_convec = package_config.natural_convec; 
		
	/* calculate flow parameters into p */
	if (!natural_convec) /* forced convection */
		calculate_flow(&p, &package_config, thermal_config);
	else /* natural convection */
		calc_natural_convec(&p, &package_config, thermal_config, sink_temp);
		
	/* print flow parameters for debug */
	#if VERBOSE > 1
	debug_print_convection(&p);
	#endif
	/* assign new r_convec calculated from the package model */
	thermal_config->r_convec = p.r_th;
	
	return natural_convec;
}

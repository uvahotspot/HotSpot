#ifndef __TEMPERATURE_H_
#define __TEMPERATURE_H_

#include "flp.h"
#include "util.h"
#include "microchannel.h"
#include "materials.h"

/* temperature sanity check threshold for natural convection and thermal runaway*/
#define TEMP_HIGH	500.0

/* model type	*/
#define BLOCK_MODEL		0
#define GRID_MODEL		1
#define BLOCK_MODEL_STR	"block"
#define GRID_MODEL_STR	"grid"

/* grid to block mapping mode	*/
#define	GRID_AVG	0
#define	GRID_MIN	1
#define	GRID_MAX	2
#define	GRID_CENTER	3
#define	GRID_AVG_STR	"avg"
#define	GRID_MIN_STR	"min"
#define	GRID_MAX_STR	"max"
#define	GRID_CENTER_STR	"center"

/* temperature-leakage loop constants */
#define LEAKAGE_MAX_ITER 100 /* max thermal-leakage iteration number, if exceeded, report thermal runaway*/
#define LEAK_TOL	0.01 /* thermal-leakage temperature convergence criterion */

/* number of extra nodes due to the model:
 * 4 spreader nodes, 4 heat sink nodes under
 * the spreader (center), 4 peripheral heat
 * sink nodes (north, south, east and west)
 * and a separate node for the ambient
 */
#define EXTRA		12
/* spreader nodes	*/
#define	SP_W		0
#define	SP_E		1
#define	SP_N		2
#define	SP_S		3
/* central sink nodes (directly under the spreader) */
#define SINK_C_W	4
#define SINK_C_E	5
#define SINK_C_N	6
#define SINK_C_S	7
/* peripheral sink nodes	*/
#define	SINK_W		8
#define	SINK_E		9
#define	SINK_N		10
#define	SINK_S		11

/* secondary extra nodes */
#define EXTRA_SEC		16
/* package substrate nodes	*/
#define	SUB_W		12
#define	SUB_E		13
#define	SUB_N		14
#define	SUB_S		15
/* solder ball nodes	*/
#define	SOLDER_W		16
#define	SOLDER_E		17
#define	SOLDER_N		18
#define	SOLDER_S		19
/* central PCB nodes (directly under the solder balls) */
#define PCB_C_W	20
#define PCB_C_E	21
#define PCB_C_N	22
#define PCB_C_S	23
/* peripheral PCB nodes	*/
#define	PCB_W		24
#define	PCB_E		25
#define	PCB_N		26
#define	PCB_S		27

/* physical constants, now become configurable in the config file	*/
//#define RHO_SI	0.01	/* thermal resistivity of silicon between 300K-400K in (mK)/W	*/
//#define	RHO_CU	0.0025	/* thermal resistivity of copper between 300K-400K in (mK)/W	*/
//#define RHO_INT	0.25	/* thermal resistivity of the interface material in (mK)/W	*/
//#define K_SI	(1.0/RHO_SI)	/* thermal conductivity of silicon	*/
//#define K_CU	(1.0/RHO_CU)	/* thermal conductivity of copper	*/
//#define K_INT	(1.0/RHO_INT)	/* thermal conductivity of the interface material	*/
//#define SPEC_HEAT_SI	1.75e6	/* specfic heat of silicon in J/(m^3K)	*/
//#define SPEC_HEAT_CU	3.55e6	/* specific heat of copper in J/(m^3K)	*/
//#define SPEC_HEAT_INT	4e6		/* specific heat of the interface material in J/(m^3K)	*/

/* typical thermal properties for secondary path layers */
/* */
#define RHO_METAL	0.0025
#define RHO_DIELECTRIC	1.0
#define RHO_C4	0.8
#define RHO_UNDERFILL 0.03
#define	RHO_SUB	0.5
#define RHO_SOLDER	0.06
#define RHO_PCB	0.333
#define K_METAL	(1.0/RHO_METAL)
#define K_DIELECTRIC (1.0/RHO_DIELECTRIC)
#define K_C4	(1.0/RHO_C4)
#define K_UNDERFILL (1.0/RHO_UNDERFILL)
#define K_SUB	(1.0/RHO_SUB)
#define K_SOLDER	(1.0/RHO_SOLDER)
#define K_PCB	(1.0/RHO_PCB)
//FIXME! the following values need to be found
#define SPEC_HEAT_METAL	3.55e6	/* specfic heat of silicon in J/(m^3K)	*/
#define SPEC_HEAT_DIELECTRIC	2.2e6
#define SPEC_HEAT_C4	1.65e6
#define SPEC_HEAT_UNDERFILL	2.65e6
#define SPEC_HEAT_SUB	1.6e6	/* specific heat of copper in J/(m^3K)	*/
#define SPEC_HEAT_SOLDER	2.1e6		/* specific heat of the interface material in J/(m^3K)	*/
#define SPEC_HEAT_PCB	1.32e6

/* model specific constants	*/
/* changed from 1/2 to 1/3 due to the difference from traditional Elmore Delay scenario */
#define C_FACTOR	0.333		/* fitting factor to match floworks (due to lumping)	*/

/* constants related to transient temperature calculation	*/
#define MIN_STEP	1e-7	/* 0.1 us	*/

/* BLAS/LAPACK definitions	*/
#define MA_NONE		0
#define MA_INTEL	1
#define MA_AMD		2
#define MA_APPLE	3
#define MA_SUN		4

#if (MATHACCEL == MA_INTEL)
	#include <mkl_cblas.h>
	#include <mkl_lapack.h>
#elif (MATHACCEL == MA_AMD)
	#include <acml.h>
#elif (MATHACCEL == MA_APPLE)
	#include <vecLib/cblas.h>
	#include <vecLib/clapack.h>
#elif (MATHACCEL == MA_SUN)
	#include <sunperf.h>
#endif

/* thermal model configuration	*/
typedef struct thermal_config_t_st
{
	/* chip specs	*/
	double t_chip;	/* chip thickness in meters	*/
	double k_chip;	/* chip thermal conductivity */
	double p_chip;	/* chip specific heat */
	double thermal_threshold;	/* temperature threshold for DTM (Kelvin)*/

	/* heat sink specs	*/
	double c_convec;	/* convection capacitance in J/K */
	double r_convec;	/* convection resistance in K/W	*/
	double s_sink;	/* heatsink side in meters	*/
	double t_sink;	/* heatsink thickness in meters	*/
	double k_sink;	/* heatsink thermal conductivity */
	double p_sink;	/* heatsink specific heat */

	/* heat spreader specs	*/
	double s_spreader;	/* spreader side in meters	*/
	double t_spreader;	/* spreader thickness in meters	*/
	double k_spreader;	/* spreader thermal conductivity */
	double p_spreader;	/* spreader specific heat */

	/* interface material specs	*/
	double t_interface;	/* interface material thickness in meters	*/
	double k_interface;	/* interface material thermal conductivity */
	double p_interface; /* interface material specific heat */

	/* secondary path specs */
	int model_secondary;
	double r_convec_sec;
	double c_convec_sec;
	int n_metal;
	double t_metal;
	double t_c4;
	double s_c4;
	int n_c4;
	double s_sub;
	double t_sub;
	double s_solder;
	double t_solder;
	double s_pcb;
	double t_pcb;

	/* others	*/
	double ambient;			/* ambient temperature in kelvin	*/
	/* initial temperatures	from file	*/
	char init_file[STR_SIZE];
	double init_temp;		/* if init_file is NULL	*/
	/* steady state temperatures to file	*/
	char steady_file[STR_SIZE];
	double sampling_intvl;	/* interval per call to compute_temp	*/
	double base_proc_freq;	/* in Hz	*/
	int dtm_used;			/* flag to guide the scaling of init Ts	*/
	/* model type - block or grid */
	char model_type[STR_SIZE];

	/* temperature-leakage loop */
	int leakage_used;
	int leakage_mode;

	/* package model */
	int package_model_used; /* flag to indicate whether package model is used */
	char package_config_file[STR_SIZE]; /* package/fan configurations */

	/* parameters specific to block model	*/
	int block_omit_lateral;	/* omit lateral resistance?	*/

	/* parameters specific to grid model	*/
	int grid_rows;			/* grid resolution - no. of rows	*/
	int grid_cols;			/* grid resolution - no. of cols	*/
	/* layer configuration from	file */
	char grid_layer_file[STR_SIZE];
	/* output grid temperatures instead of block temperatures */
	char grid_steady_file[STR_SIZE];
	/* mapping mode between grid and block models	*/
	char grid_map_mode[STR_SIZE];
  /* transient grid temperatures to file */
  char grid_transient_file[STR_SIZE];

	int detailed_3D_used; //BU_3D: Added parameter to check for heterogenous R-C model
}thermal_config_t;

/* defaults	*/
thermal_config_t default_thermal_config(void);
/*
 * parse a table of name-value string pairs and add the configuration
 * parameters to 'config'
 */
void thermal_config_add_from_strs(thermal_config_t *config, materials_list_t *materials_list, str_pair *table, int size);
/*
 * convert config into a table of name-value pairs. returns the no.
 * of parameters converted
 */
int thermal_config_to_strs(thermal_config_t *config, str_pair *table, int max_entries);

/* package parameters	*/
typedef struct package_RC_t_st
{
	/* lateral resistances	*/
	/* peripheral spreader nodes */
	double r_sp1_x;
	double r_sp1_y;
	/* sink's inner periphery */
	double r_hs1_x;
	double r_hs1_y;
	double r_hs2_x;
	double r_hs2_y;
	/* sink's outer periphery */
	double r_hs;

	/* vertical resistances */
	/* peripheral spreader nodes */
	double r_sp_per_x;
	double r_sp_per_y;
	/* sink's inner periphery	*/
	double r_hs_c_per_x;
	double r_hs_c_per_y;
	/* sink's outer periphery	*/
	double r_hs_per;

	/* vertical capacitances	*/
	/* peripheral spreader nodes */
	double c_sp_per_x;
	double c_sp_per_y;
	/* sink's inner periphery	*/
	double c_hs_c_per_x;
	double c_hs_c_per_y;
	/* sink's outer periphery	*/
	double c_hs_per;

	/* vertical R's and C's to ambient	*/
	/* sink's inner periphery	*/
	double r_amb_c_per_x;
	double c_amb_c_per_x;
	double r_amb_c_per_y;
	double c_amb_c_per_y;
	/* sink's outer periphery	*/
	double r_amb_per;
	double c_amb_per;

	/* secondary path R's and C's */

	/* lateral resistances	*/
	/* peripheral package substrate nodes */
	double r_sub1_x;
	double r_sub1_y;
	/* peripheral solder ball nodes */
	double r_solder1_x;
	double r_solder1_y;
	/* PCB's inner periphery */
	double r_pcb1_x;
	double r_pcb1_y;
	double r_pcb2_x;
	double r_pcb2_y;
	/* PCB's outer periphery */
	double r_pcb;

	/* vertical resistances */
	/* peripheral package substrate nodes */
	double r_sub_per_x;
	double r_sub_per_y;
	/* peripheral solder ball nodes */
	double r_solder_per_x;
	double r_solder_per_y;
	/* PCB's inner periphery	*/
	double r_pcb_c_per_x;
	double r_pcb_c_per_y;
	/* PCB's outer periphery	*/
	double r_pcb_per;

	/* vertical capacitances	*/
	/* peripheral package substrate nodes */
	double c_sub_per_x;
	double c_sub_per_y;
	/* peripheral solder ballnodes */
	double c_solder_per_x;
	double c_solder_per_y;
	/* PCB's inner periphery	*/
	double c_pcb_c_per_x;
	double c_pcb_c_per_y;
	/* PCB's outer periphery	*/
	double c_pcb_per;

	/* vertical R's and C's to ambient at PCB	*/
	/* PCB's inner periphery	*/
	double r_amb_sec_c_per_x;
	double c_amb_sec_c_per_x;
	double r_amb_sec_c_per_y;
	double c_amb_sec_c_per_y;
	/* PCB's outer periphery	*/
	double r_amb_sec_per;
	double c_amb_sec_per;

}package_RC_t;
/* package parameter routines	*/
void populate_package_R(package_RC_t *p, thermal_config_t *config, double width, double height);
void populate_package_C(package_RC_t *p, thermal_config_t *config, double width, double height);
/* debug print	*/
void debug_print_package_RC(package_RC_t *p);

/* slope function pointer - used as a call back by the transient solver	*/
typedef void (*slope_fn_ptr)(void *model, void *y, void *p, void *dy);

/* hotspot thermal model - can be a block or grid model	*/
struct block_model_t_st;
struct grid_model_t_st;
typedef struct RC_model_t_st
{
	union
	{
		struct block_model_t_st *block;
		struct grid_model_t_st *grid;
	};
	/* block model or grid model	*/
	int type;
	thermal_config_t *config;
}RC_model_t;

/* constructor/destructor	*/
/* placeholder is an empty floorplan frame with only the names of the functional units	*/
RC_model_t *alloc_RC_model(thermal_config_t *config, flp_t *placeholder, microchannel_config_t *microchannel_config, materials_list_t *materials_list,
	                         int do_detailed_3D, int use_microchannels);
void delete_RC_model(RC_model_t *model);

/* initialization	*/
void populate_R_model(RC_model_t *model, flp_t *flp);
void populate_C_model(RC_model_t *model, flp_t *flp);

/* hotspot main interfaces - temperature.c	*/
void steady_state_temp(RC_model_t *model, double *power, double *temp);
void compute_temp(RC_model_t *model, double *power, double *temp, double time_elapsed);
/* differs from 'dvector()' in that memory for internal nodes is also allocated	*/
double *hotspot_vector(RC_model_t *model);
/* copy 'src' to 'dst' except for a window of 'size'
 * elements starting at 'at'. useful in floorplan
 * compaction
 */
void trim_hotspot_vector(RC_model_t *model, double *dst, double *src,
						 int at, int size);
/* update the model's node count	*/
void resize_thermal_model(RC_model_t *model, int n_units);
void set_temp (RC_model_t *model, double *temp, double val);
void dump_temp (RC_model_t *model, double *temp, char *file);
void copy_temp (RC_model_t *model, double *dst, double *src);
void read_temp (RC_model_t *model, double *temp, char *file, int clip);
void dump_power(RC_model_t *model, double *power, char *file);
void read_power (RC_model_t *model, double *power, char *file);
double find_max_temp(RC_model_t *model, double *temp);
double find_avg_temp(RC_model_t *model, double *temp);
/* debug print	*/
void debug_print_model(RC_model_t *model);

/* other functions used by the above interfaces	*/

/* LUP decomposition	*/
void lupdcmp(double**a, int n, int *p, int spd);

/* get the thermal resistance and capacitance values	*/
double getr(double conductivity, double thickness, double area);
double getcap(double sp_heat, double thickness, double area);

/* LU forward and backward substitution	*/
void lusolve(double **a, int n, int *p, double *b, double *x, int spd);

/* 4th order Runge Kutta solver with adaptive step sizing */
double rk4(void *model, double *y, void *p, int n, double *h, double *yout, slope_fn_ptr f);

#if SUPERLU > 0
/* Backward Euler solver with adaptive step sizing */
double backward_euler(SuperMatrix *G, diagonal_matrix_t *C, double *T, double *P, double *h, double *Tout);
#endif

/* matrix and vector routines	*/
void matmult(double **c, double **a, double **b, int n);
/* same as above but 'a' is a diagonal matrix stored as a 1-d array	*/
void diagmatmult(double **c, double *a, double **b, int n);
void matvectmult(double *vout, double **m, double *vin, int n);
/* same as above but 'm' is a diagonal matrix stored as a 1-d array	*/
void diagmatvectmult(double *vout, double *m, double *vin, int n);
/*
 * inv = m^-1, inv, m are n by n matrices.
 * the spd flag indicates that m is symmetric
 * and positive definite
 */
void matinv(double **inv, double **m, int n, int spd);

/* dst = src1 + scale * src2	*/
void scaleadd_dvector (double *dst, double *src1, double *src2, int n, double scale);

/* temperature-aware leakage calculation */
double calc_leakage(int mode, double h, double w, double temp);

/* calculate average heatsink temperature for natural convection package model */
double calc_sink_temp(RC_model_t *model, double *temp);

#endif

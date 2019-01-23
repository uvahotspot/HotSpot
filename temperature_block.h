#ifndef __TEMPERATURE_BLOCK_H_
#define __TEMPERATURE_BLOCK_H_

#include "temperature.h"

/* functional block layers	*/
/* total	*/
#define NL	4
/* silicon is always layer 0	*/
/* interface layer	*/
#define IFACE	1
/* heat spreader	*/
#define HSP	2
/* heat sink */
#define HSINK 3

/* block thermal model	*/
typedef struct block_model_t_st
{
	/* floorplan	*/
	flp_t *flp;

	/* configuration	*/
	thermal_config_t config;

	/* main matrices	*/
	/* conductance matrix */
	double **b;
	/* LUP decomposition of b	*/
	double **lu; 
	int *p;
	/* diagonal capacitance matrix stored as a 1-d vector	*/
	double *a; 
	/* inverse of the above	*/
	double *inva;
	/* c = inva * b	*/
	double **c;

	/* package parameters	*/
	package_RC_t pack;

	/* intermediate vectors and matrices	*/
	double *gx, *gy;
	double *gx_int, *gy_int;
	double *gx_sp, *gy_sp;
	double *gx_hs, *gy_hs;
	double *g_amb;
	double *t_vector;
	double **len, **g;
	int **border;

	/* total no. of nodes	*/
	int n_nodes;
	/* total no. of blocks	*/
	int n_units;
	/* to allow for resizing	*/
	int base_n_units;

	/* flags	*/
	int r_ready;	/* are the R's initialized?	*/
	int c_ready;	/* are the C's initialized?	*/
}block_model_t;

/* constructor/destructor	*/
/* placeholder is an empty floorplan frame with only the names of the functional units	*/
block_model_t *alloc_block_model(thermal_config_t *config, flp_t *placeholder);
void delete_block_model(block_model_t *model);

/* initialization	*/
void populate_R_model_block(block_model_t *model, flp_t *flp);
void populate_C_model_block(block_model_t *model, flp_t *flp);

/* hotspot main interfaces - temperature.c	*/
void steady_state_temp_block(block_model_t *model, double *power, double *temp);
void compute_temp_block(block_model_t *model, double *power, double *temp, double time_elapsed);
/* differs from 'dvector()' in that memory for internal nodes is also allocated	*/
double *hotspot_vector_block(block_model_t *model);
/* copy 'src' to 'dst' except for a window of 'size'
 * elements starting at 'at'. useful in floorplan
 * compaction
 */
void trim_hotspot_vector_block(block_model_t *model, double *dst, double *src, 
						 	   int at, int size);
/* update the model's node count	*/						 
void resize_thermal_model_block(block_model_t *model, int n_units);						 
void set_temp_block (block_model_t *model, double *temp, double val);
void dump_temp_block (block_model_t *model, double *temp, char *file);
void copy_temp_block (block_model_t *model, double *dst, double *src);
void read_temp_block (block_model_t *model, double *temp, char *file, int clip);
void dump_power_block(block_model_t *model, double *power, char *file);
void read_power_block (block_model_t *model, double *power, char *file);
double find_max_temp_block(block_model_t *model, double *temp);
double find_avg_temp_block(block_model_t *model, double *temp);
double calc_sink_temp_block(block_model_t *model, double *temp, thermal_config_t *config); //for natural convection package model
/* debug print	*/
void debug_print_block(block_model_t *model);

#endif

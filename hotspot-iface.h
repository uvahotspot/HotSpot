#ifndef __HOTSPOT_IFACE_H_
#define __HOTSPOT_IFACE_H_

#define STR_SIZE	512
#define MAX_UNITS	8192
#define NULLFILE	"(null)"
#define BLOCK_MODEL		0
#define GRID_MODEL		1

/* floorplan structure	*/
typedef struct unit_t_st
{
	char name[STR_SIZE];
	double width;
	double height;
	double leftx;
	double bottomy;
}unit_t;
typedef struct flp_t_st
{
	unit_t *units;
	int n_units;
  	double **wire_density;
} flp_t;

/* floorplan routines	*/

/* reads the floorplan from a file and allocates memory. 
 * 'read_connects' is a boolean flag indicating if 
 * connectivity information should also be read (usually 
 * set to FALSE).
 */ 
flp_t *read_flp(char *file, int read_connects);

/* deletes floorplan from memory. 'compacted' is a 
 * boolean flag indicating if this is a floorplan
 * compacted by HotFloorplan (usually set to FALSE).
 */ 
void free_flp(flp_t *flp, int compacted);

/* thermal model configuration structure	*/
typedef struct thermal_config_t_st
{
	double t_chip;	/* chip thickness in meters	*/
	double k_chip;	/* chip thermal conductivity */
	double p_chip;	/* chip specific heat */
	double thermal_threshold;	/* temperature threshold for DTM (Kelvin)*/
	double c_convec;	/* convection capacitance in J/K */
	double r_convec;	/* convection resistance in K/W	*/
	double s_sink;	/* heatsink side in meters	*/
	double t_sink;	/* heatsink thickness in meters	*/
	double k_sink;	/* heatsink thermal conductivity */
	double p_sink;	/* heatsink specific heat */
	double s_spreader;	/* spreader side in meters	*/
	double t_spreader;	/* spreader thickness in meters	*/
	double k_spreader;	/* spreader thermal conductivity */
	double p_spreader;	/* spreader specific heat */
	double t_interface;	/* interface material thickness in meters	*/
	double k_interface;	/* interface material thermal conductivity */
	double p_interface; /* interface material specific heat */
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
	double ambient;			/* ambient temperature in kelvin	*/
	char init_file[STR_SIZE];
	double init_temp;		/* if init_file is NULL	*/
	char steady_file[STR_SIZE];
	double sampling_intvl;	/* interval per call to compute_temp	*/
	double base_proc_freq;	/* in Hz	*/
	int dtm_used;			/* flag to guide the scaling of init Ts	*/
	char model_type[STR_SIZE];
	int leakage_used;
	int leakage_mode;
	int package_model_used; /* flag to indicate whether package model is used */
	char package_config_file[STR_SIZE]; /* package/fan configurations */ 
	int block_omit_lateral;	/* omit lateral resistance?	*/
	int grid_rows;			/* grid resolution - no. of rows	*/
	int grid_cols;			/* grid resolution - no. of cols	*/
	char grid_layer_file[STR_SIZE];
	char grid_steady_file[STR_SIZE];
	char grid_map_mode[STR_SIZE];
	int detailed_3D_used; //BU_3D: Added parameter to check for heterogenous R-C model 
}thermal_config_t;

/* thermal configuration routines */

/* returns a thermal configuration structure with
 * the default set of parameters
 */ 
thermal_config_t default_thermal_config(void);

/* thermal model structure	*/
struct block_model_t_st;
struct grid_model_t_st;
typedef struct RC_model_t_st
{
	union
	{
		struct block_model_t_st *block;
		struct grid_model_t_st *grid;
	};
	int type;
	thermal_config_t *config;
}RC_model_t;

/* thermal model routines */

/* creates a new thermal model. 'config' can be obtained 
 * from the 'default_thermal_config' function and 
 * 'placeholder' can be obtained from the 'read_flp' 
 * function
 */ 
RC_model_t *alloc_RC_model(thermal_config_t *config, flp_t *placeholder, int do_detailed_3D);

/* deletes the thermal model and frees up memory	*/
void delete_RC_model(RC_model_t *model);

/* populates the thermal resistances of the model. This is
 * a prerequisite for computing transient or steady state
 * temperatures.
 */ 
void populate_R_model(RC_model_t *model, flp_t *flp);

/* populates the thermal capacitances of the model. This is
 * a prerequisite for computing transient temperatures.
 */ 
void populate_C_model(RC_model_t *model, flp_t *flp);

/* memory allocator for the power and temperature vectors	*/
double *hotspot_vector(RC_model_t *model);

/* destructor for a vector allocated using 'hotspot_vector'	*/
void free_dvector(double *v);

/* outputs the 'temp' vector onto 'file'. 'temp' must
 * be allocated using ' hotspot_vector'.
 */
void dump_temp (RC_model_t *model, double *temp, char *file);

/* sets all the temperatures of the 'temp' vector to the
 * value 'val'. 'temp' must be allocated using '
 * hotspot_vector'.
 */
void set_temp (RC_model_t *model, double *temp, double val);

/* read the 'temp' vector from 'file'. The format of the 
 * file should be the same as the one output by the 
 * 'dump_temp' function. 'temp' must be allocated using 
 * 'hotspot_vector'. 'clip' is a boolean flag indicating
 * whether to clip the peak temperature of the vector to
 * the thermal threshold 'model->config->thermal_threshold'
 * (usually set to FALSE).
 */
void read_temp (RC_model_t *model, double *temp, char *file, int clip);

/* computation of the steady state temperatures. 'power'
 * and 'temp' must be allocated using 'hotspot_vector'.
 * 'populate_R_model' must be called before this.
 * 'power' should contain the input power numbers. 'temp'
 * will contain the output temperature numbers after the 
 * call.
 */ 
void steady_state_temp(RC_model_t *model, double *power, double *temp);

/* computation of the transient temperatures. 'power'
 * and 'temp' must be allocated using 'hotspot_vector'.
 * 'populate_R_model' and 'populate_C_model' must be
 * called before this. 'power' should  contain the 
 * input power numbers and 'temp' should contain the 
 * current temperatures. 'time_elapsed' is the duration
 * of the transient simulation. 'temp' will contain the
 * output temperature numbers after the call.
 */ 
void compute_temp(RC_model_t *model, double *power, double *temp, double time_elapsed);

#endif

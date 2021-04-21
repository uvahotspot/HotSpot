#ifndef __FLP_H_
#define __FLP_H_

#include "util.h"

#define MAX_UNITS		8192
#define MAX_MOVES		16

/* types of cuts	*/
#define CUT_NONE		-1
#define CUT_VERTICAL	-2
#define	CUT_HORIZONTAL	-3

/* wrap around L2 with extra arms	*/
#define L2_LEFT			0
#define L2_RIGHT		1
#define L2_ARMS			2
#define L2_LEFT_STR		"_left"
#define L2_RIGHT_STR	"_right"
/* L2 sizing ratio of arm width to base height	*/
#define WRAP_L2_RATIO	5

/*
 * chip edge has true dead space, which is
 * modeled by the following blocks
 */
#define RIM_LEFT		1
#define RIM_RIGHT		2
#define RIM_TOP			4
#define RIM_BOTTOM		8
#define RIM_PREFIX		"RIM"
#define RIM_LEFT_STR	RIM_PREFIX"_left"
#define RIM_RIGHT_STR	RIM_PREFIX"_right"
#define RIM_TOP_STR		RIM_PREFIX"_top"
#define RIM_BOTTOM_STR	RIM_PREFIX"_bottom"

/* prefix denoting dead block	*/
#define DEAD_PREFIX		"_"

/* flags denoting orientation	*/
/* rotated orientations	*/
#define	ROT_0		0x01	/* normal	*/
#define	ROT_90		0x02	/* 90 degrees anticlockwise */
#define	ROT_180		0x04	/* 180 degrees anticlockwise */
#define	ROT_270		0x08	/* 270 degrees anticlockwise */
/* flipped + rotated orientations	*/
#define	FLIP_0		0x10	/* flip about y axis of ROT_0	*/
#define	FLIP_90		0x20	/* flip about y axis of ROT_90	*/
#define	FLIP_180	0x40	/* flip about y axis of ROT_180	*/
#define	FLIP_270	0x80	/* flip about y axis of ROT_270	*/
#define ORIENTS_N	8		/* total no. of orientations	*/
#define ORIENT_NDEF	0xDF	/* undefined orientation	*/

/* type for holding the above flags	*/
typedef unsigned char orient_t;

/* forward declarations	*/
struct RC_model_t_st;	/* see temperature.h	*/
struct shape_t_st;		/* see shape.h	*/

/* configuration parameters for the floorplan	*/
typedef struct flp_config_t_st
{
	/* wrap around L2?	*/
	int wrap_l2;
	/* name of L2 to look for	*/
	char l2_label[STR_SIZE];

	/* model dead space around the rim of the chip? */
	int model_rim;
	double rim_thickness;

	/* area ratio below which to ignore dead space	*/
	double compact_ratio;

	/*
	 * no. of discrete orientations for a shape curve.
	 * should be an even number greater than 1
	 */
	int n_orients;

	/* annealing parameters	*/
	double P0;		/* initial acceptance probability	*/
	double Davg;	/* average change (delta) in cost	*/
	double Kmoves;		/* no. of moves to try in each step	*/
	double Rcool;	/* ratio for the cooling schedule */
	double Rreject;	/* ratio of rejects at which to stop annealing */
	int Nmax;		/* absolute max no. of annealing steps	*/

	/* weights for the metric: lambdaA * A + lambdaT * T + lambdaW * W	*/
	double lambdaA;
	double lambdaT;
	double lambdaW;
} flp_config_t;

/* unplaced unit	*/
typedef struct unplaced_t_st
{
  char name[STR_SIZE];
  /* can be rotated?	*/
  int rotable;
  double area;
  /* minimum and maximum aspect ratios	*/
  double min_aspect;
  double max_aspect;
  /* shape curve for this unit	*/
  struct shape_t_st *shape;
}unplaced_t;

/* input description for floorplanning	*/
typedef struct flp_desc_t_st
{
  unplaced_t *units;
  /* density of wires between units	*/
  double **wire_density;
  /* configuration parameters	*/
  flp_config_t config;
  int n_units;
}flp_desc_t;

/* placed functional unit */
typedef struct unit_t_st
{
	char name[STR_SIZE];
	double width;
	double height;
	double leftx;
	double bottomy;
	//BU_3D: new parameters for each functional unit. Allows for specific R-C values at a block-level.
	double specificheat;
	double resistivity;
	int hasRes;
	int hasSh;
	//end->BU_3D
}unit_t;

/* floorplan data structure	*/
typedef struct flp_t_st
{
	unit_t *units;
	int n_units;
  	/* density of wires between units	*/
  	double **wire_density;
} flp_t;

/* flp_config routines	*/

/* default flp_config	*/
flp_config_t default_flp_config(void);
/*
 * parse a table of name-value string pairs and add the configuration
 * parameters to 'config'
 */
void flp_config_add_from_strs(flp_config_t *config, str_pair *table, int size);
/*
 * convert config into a table of name-value pairs. returns the no.
 * of parameters converted
 */
int flp_config_to_strs(flp_config_t *config, str_pair *table, int max_entries);

/* flp_desc routines	*/

/* read floorplan description and allocate memory	*/
flp_desc_t *read_flp_desc(char *file, flp_config_t *config);
void free_flp_desc(flp_desc_t *flp_desc);
/* debug print	*/
void print_unplaced(unplaced_t *unit);
void print_flp_desc(flp_desc_t *flp_desc);

/* flp routines	*/

/* create a floorplan placeholder from description	*/
flp_t *flp_placeholder(flp_desc_t *flp_desc);
/* skip floorplanning and read floorplan directly from file */
flp_t *read_flp(char *file, int read_connects, int initialize_connects);
/*
 * main flooplanning routine - allocates
 * memory internally. returns the number
 * of compacted blocks
 */
int floorplan(flp_t *flp, flp_desc_t *flp_desc,
			  struct RC_model_t_st *model, double *power);
/*
 * print the floorplan in a FIG like format
 * that can be read by tofig.pl to produce
 * an xfig output
 */
void print_flp_fig (flp_t *flp);
/* debug print	*/
void print_flp (flp_t *flp, int print_connects);
/* print the statistics about this floorplan	*/
void print_flp_stats(flp_t *flp, struct RC_model_t_st *model,
					 char *l2_label, char *power_file,
					 char *connects_file);
/* wrap the L2 around this floorplan	*/
void flp_wrap_l2(flp_t *flp, flp_desc_t *flp_desc);
/* wrap the rim blocks around - returns the no. of blocks	*/
int flp_wrap_rim(flp_t *flp, double rim_thickness);
/* translate the floorplan to new origin (x,y)	*/
void flp_translate(flp_t *flp, double x, double y);
/* scale the floorplan by a factor 'factor'	*/
void flp_scale(flp_t *flp, double factor);
/*
 * change the orientation of the floorplan by
 * rotating and/or flipping. the target orientation
 * is specified in 'target'. 'width', 'height', 'xorig'
 * and 'yorig' are those of 'flp' respectively.
 */
void flp_change_orient(flp_t *flp, double xorig, double yorig,
					   double width, double height, orient_t target);

/*
 * create a non-uniform grid-like floorplan equivalent to this.
 * this function is mainly useful when using the HotSpot block
 * model to model floorplans of drastically differing aspect
 * ratios and granularity. an example for such a floorplan
 * would be the standard ev6 floorplan that comes with HotSpot,
 * where the register file is subdivided into say 128 entries.
 * the HotSpot block model could result in inaccuracies while
 * trying to model such floorplans of differing granularity.
 * if such inaccuracies occur, use this function to create an
 * equivalent floorplan that can be modeled accurately in
 * HotSpot. 'map' is an output parameter to store the 2-d array
 * allocated by the function.
 */
flp_t *flp_create_grid(flp_t *flp, int ***map);
/* free the map allocated by flp_create_grid	*/
void free_blkgrid_map(flp_t *flp, int **map);
/* translate power numbers to the grid created by flp_create_grid	*/
void xlate_power_blkgrid(flp_t *flp, flp_t *grid, \
						 double *bpower, double *gpower, int **map);
/* the metric used to evaluate the floorplan	*/
double flp_evaluate_metric(flp_t *flp, struct RC_model_t_st *model, double *power,
						   double lambdaA, double lambdaT, double lambdaW);
/* dump the floorplan onto a file	*/
void dump_flp(flp_t *flp, char *file, int dump_connects);
/* memory uninitialization	*/
void free_flp(flp_t *flp, int compacted, int free_connects);


/* placed floorplan access routines	*/

/* get unit index from its name	*/
int get_blk_index(flp_t *flp, char *name);
/* are the units horizontally adjacent?	*/
int is_horiz_adj(flp_t *flp, int i, int j);
/* are the units vertically adjacent?	*/
int is_vert_adj (flp_t *flp, int i, int j);
/* shared length between units	*/
double get_shared_len(flp_t *flp, int i, int j);
/* total chip width	*/
double get_total_width(flp_t *flp);
/* total chip height */
double get_total_height(flp_t *flp);
/* x and y origins	*/
double get_minx(flp_t *flp);
double get_miny(flp_t *flp);
/* precondition: L2 should have been wrapped around	*/
double get_core_width(flp_t *flp, char *l2_label);
/* precondition: L2 should have been wrapped around	*/
double get_core_height(flp_t *flp, char *l2_label);
/* other queries	*/
double get_manhattan_dist(flp_t *flp, int i, int j);
double get_total_area(flp_t *flp);
double get_core_area(flp_t *flp, char *l2_label);
double get_core_occupied_area(flp_t *flp, char *l2_label);
double get_wire_metric(flp_t *flp);

#endif

#include <stdio.h>
#include <string.h>
#include <string.h>
#ifdef _MSC_VER
#define strcasecmp    _stricmp
#define strncasecmp   _strnicmp
#else
#include <strings.h>
#endif
#include <stdlib.h>
#include <math.h>

#include "flp.h"
#include "npe.h"
#include "shape.h"
#include "util.h"
#include "temperature.h"
#include "temperature_block.h"

/*
 * this is the metric function used for the floorplanning.
 * in order to enable a different metric, just change the
 * return statement of this function to return an appropriate
 * metric. The current metric used is a linear function of
 * area (A), temperature (T) and wire length (W):
 * lambdaA * A + lambdaT * T  + lambdaW * W
 * thermal model and power density are passed as parameters
 * since temperature is used in the metric.
 */
double flp_evaluate_metric(flp_t *flp, RC_model_t *model, double *power,
						   double lambdaA, double lambdaT, double lambdaW)
{
	double tmax, area, wire_length;
	double *temp;

	temp = hotspot_vector(model);
	populate_R_model(model, flp);
	steady_state_temp(model, power, temp);
	tmax = find_max_temp(model, temp);
	area = get_total_area(flp);
	wire_length = get_wire_metric(flp);
	free_dvector(temp);

	/* can return any arbitrary function of area, tmax and wire_length	*/
	return (lambdaA * area + lambdaT * tmax + lambdaW * wire_length);
}


/* default flp_config	*/
flp_config_t default_flp_config(void)
{
	flp_config_t config;

	/* wrap around L2?	*/
	config.wrap_l2 = TRUE;
	strcpy(config.l2_label, "L2");

	/* model dead space around the rim of the chip? */
	config.model_rim = FALSE;
	config.rim_thickness = 50e-6;

	/* area ratio below which to ignore dead space	*/
	config.compact_ratio = 0.005;

	/*
	 * no. of discrete orientations for a shape curve.
	 * should be an even number greater than 1
	 */
	config.n_orients = 300;

	/* annealing parameters	*/
	config.P0 = 0.99;		/* initial acceptance probability	*/
	/*
	 * average change (delta) in cost. varies according to
	 * the metric. need not be very accurate. just the right
	 * order of magnitude is enough. for instance, if the
	 * metric is flp area, this Davg is the average difference
	 * in the area of successive slicing floorplan attempts.
	 * since the areas are in the order of mm^2, the delta
	 * is also in the same ball park.
	 */
	config.Davg = 1.0;		/* for our a*A + b*T + c*W metric	*/
	config.Kmoves = 7.0;	/* no. of moves to try in each step	*/
	config.Rcool = 0.99;	/* ratio for the cooling schedule */
	config.Rreject = 0.99;	/* ratio of rejects at which to stop annealing */
	config.Nmax = 1000;		/* absolute max no. of annealing steps	*/

	/* weights for the metric: lambdaA * A + lambdaT * T + lambdaW * W
	 * the weights incorporate two things:
	 * 1) the conversion of A to mm^2, T to K and W to mm.
	 * 2) weighing the relative importance of A, T and K
	 */
	config.lambdaA = 5.0e+6;
	config.lambdaT = 1.0;
	config.lambdaW = 350;

	return config;
}

/*
 * parse a table of name-value string pairs and add the configuration
 * parameters to 'config'
 */
void flp_config_add_from_strs(flp_config_t *config, str_pair *table, int size)
{
	int idx;
	if ((idx = get_str_index(table, size, "wrap_l2")) >= 0)
		if(sscanf(table[idx].value, "%d", &config->wrap_l2) != 1)
			fatal("invalid format for configuration  parameter wrap_l2\n");
	if ((idx = get_str_index(table, size, "l2_label")) >= 0)
		if(sscanf(table[idx].value, "%s", config->l2_label) != 1)
			fatal("invalid format for configuration  parameter l2_label\n");
	if ((idx = get_str_index(table, size, "model_rim")) >= 0)
		if(sscanf(table[idx].value, "%d", &config->model_rim) != 1)
			fatal("invalid format for configuration  parameter model_rim\n");
	if ((idx = get_str_index(table, size, "rim_thickness")) >= 0)
		if(sscanf(table[idx].value, "%lf", &config->rim_thickness) != 1)
			fatal("invalid format for configuration  parameter rim_thickness\n");
	if ((idx = get_str_index(table, size, "compact_ratio")) >= 0)
		if(sscanf(table[idx].value, "%lf", &config->compact_ratio) != 1)
			fatal("invalid format for configuration  parameter compact_ratio\n");
	if ((idx = get_str_index(table, size, "n_orients")) >= 0)
		if(sscanf(table[idx].value, "%d", &config->n_orients) != 1)
			fatal("invalid format for configuration  parameter n_orients\n");
	if ((idx = get_str_index(table, size, "P0")) >= 0)
		if(sscanf(table[idx].value, "%lf", &config->P0) != 1)
			fatal("invalid format for configuration  parameter P0\n");
	if ((idx = get_str_index(table, size, "Davg")) >= 0)
		if(sscanf(table[idx].value, "%lf", &config->Davg) != 1)
			fatal("invalid format for configuration  parameter Davg\n");
	if ((idx = get_str_index(table, size, "Kmoves")) >= 0)
		if(sscanf(table[idx].value, "%lf", &config->Kmoves) != 1)
			fatal("invalid format for configuration  parameter Kmoves\n");
	if ((idx = get_str_index(table, size, "Rcool")) >= 0)
		if(sscanf(table[idx].value, "%lf", &config->Rcool) != 1)
			fatal("invalid format for configuration  parameter Rcool\n");
	if ((idx = get_str_index(table, size, "Rreject")) >= 0)
		if(sscanf(table[idx].value, "%lf", &config->Rreject) != 1)
			fatal("invalid format for configuration  parameter Rreject\n");
	if ((idx = get_str_index(table, size, "Nmax")) >= 0)
		if(sscanf(table[idx].value, "%d", &config->Nmax) != 1)
			fatal("invalid format for configuration  parameter Nmax\n");
	if ((idx = get_str_index(table, size, "lambdaA")) >= 0)
		if(sscanf(table[idx].value, "%lf", &config->lambdaA) != 1)
			fatal("invalid format for configuration  parameter lambdaA\n");
	if ((idx = get_str_index(table, size, "lambdaT")) >= 0)
		if(sscanf(table[idx].value, "%lf", &config->lambdaT) != 1)
			fatal("invalid format for configuration  parameter lambdaT\n");
	if ((idx = get_str_index(table, size, "lambdaW")) >= 0)
		if(sscanf(table[idx].value, "%lf", &config->lambdaW) != 1)
			fatal("invalid format for configuration  parameter lambdaW\n");

	if (config->rim_thickness <= 0)
		fatal("rim thickness should be greater than zero\n");
	if ((config->compact_ratio < 0) || (config->compact_ratio > 1))
		fatal("compact_ratio should be between 0 and 1\n");
	if ((config->n_orients <= 1) || (config->n_orients & 1))
		fatal("n_orients should be an even number greater than 1\n");
	if (config->Kmoves < 0)
		fatal("Kmoves should be non-negative\n");
	if ((config->P0 < 0) || (config->P0 > 1))
		fatal("P0 should be between 0 and 1\n");
	if ((config->Rcool < 0) || (config->Rcool > 1))
		fatal("Rcool should be between 0 and 1\n");
	if ((config->Rreject < 0) || (config->Rreject > 1))
		fatal("Rreject should be between 0 and 1\n");
	if (config->Nmax < 0)
		fatal("Nmax should be non-negative\n");
}

/*
 * convert config into a table of name-value pairs. returns the no.
 * of parameters converted
 */
int flp_config_to_strs(flp_config_t *config, str_pair *table, int max_entries)
{
	if (max_entries < 15)
		fatal("not enough entries in table\n");

	sprintf(table[0].name, "wrap_l2");
	sprintf(table[1].name, "l2_label");
	sprintf(table[2].name, "model_rim");
	sprintf(table[3].name, "rim_thickness");
	sprintf(table[4].name, "compact_ratio");
	sprintf(table[5].name, "n_orients");
	sprintf(table[6].name, "P0");
	sprintf(table[7].name, "Davg");
	sprintf(table[8].name, "Kmoves");
	sprintf(table[9].name, "Rcool");
	sprintf(table[10].name, "Rreject");
	sprintf(table[11].name, "Nmax");
	sprintf(table[12].name, "lambdaA");
	sprintf(table[13].name, "lambdaT");
	sprintf(table[14].name, "lambdaW");

	sprintf(table[0].value, "%d", config->wrap_l2);
	sprintf(table[1].value, "%s", config->l2_label);
	sprintf(table[2].value, "%d", config->model_rim);
	sprintf(table[3].value, "%lg", config->rim_thickness);
	sprintf(table[4].value, "%lg", config->compact_ratio);
	sprintf(table[5].value, "%d", config->n_orients);
	sprintf(table[6].value, "%lg", config->P0);
	sprintf(table[7].value, "%lg", config->Davg);
	sprintf(table[8].value, "%lg", config->Kmoves);
	sprintf(table[9].value, "%lg", config->Rcool);
	sprintf(table[10].value, "%lg", config->Rreject);
	sprintf(table[11].value, "%d", config->Nmax);
	sprintf(table[12].value, "%lg", config->lambdaA);
	sprintf(table[13].value, "%lg", config->lambdaT);
	sprintf(table[14].value, "%lg", config->lambdaW);

	return 15;
}

/*
 * copy L2 connectivity from 'from' of flp_desc to 'to'
 * of flp. 'size' elements are copied. the arms are not
 * connected amidst themselves or with L2 base block
 */
void copy_l2_info (flp_t *flp, int to, flp_desc_t *flp_desc, int from, int size)
{
	int j, count;

	for(count=0; count < L2_ARMS + 1; count++, to++) {
		/* copy names */
		strcpy(flp->units[to].name, flp_desc->units[from].name);
		for(j=0; j < size; j++) {
			/* rows	*/
			flp->wire_density[to][j] = flp_desc->wire_density[from][j];
			/* columns	*/
			flp->wire_density[j][to] = flp_desc->wire_density[j][from];
		}
	}
	/* fix the names of the arms	*/
	strcat(flp->units[to-L2_ARMS+L2_LEFT].name, L2_LEFT_STR);
	strcat(flp->units[to-L2_ARMS+L2_RIGHT].name, L2_RIGHT_STR);
}


/* create a floorplan placeholder from description	*/
flp_t *flp_placeholder(flp_desc_t *flp_desc)
{
	int i, j, count, n_dead;
	flp_t *flp;

	/* wrap L2 around?	*/
	int wrap_l2 = FALSE;
	if (flp_desc->config.wrap_l2 &&
		!strcasecmp(flp_desc->units[flp_desc->n_units-1].name, flp_desc->config.l2_label))
		wrap_l2 = TRUE;

	flp = (flp_t *) calloc (1, sizeof(flp_t));
	if(!flp)
		fatal("memory allocation error\n");
	/*
	 * number of dead blocks = no. of core blocks - 1.
	 * (one per vertical or horizontal cut). if L2 is
	 * wrapped around, core blocks = flp_desc->n_units-1
	 */
	n_dead = flp_desc->n_units - !!(wrap_l2) - 1;
	flp->n_units = flp_desc->n_units + n_dead;

	/* wrap L2 around - extra arms are added */
	if (wrap_l2)
		flp->n_units += L2_ARMS;

	/*
	 * model the dead space in the edge. let us make
	 * one dead block per corner edge of a block. so,
	 * no. of rim blocks could be at most 2*n+2 where
	 * n is the total no. of blocks (the worst case
	 * is just all blocks lined up side-by-side)
	 */
	if (flp_desc->config.model_rim)
		flp->n_units += (2*flp->n_units + 2);

	flp->units = (unit_t *) calloc (flp->n_units, sizeof(unit_t));
	flp->wire_density = (double **) calloc(flp->n_units, sizeof(double *));
	if (!flp->units || !flp->wire_density)
		fatal("memory allocation error\n");
	for (i=0; i < flp->n_units; i++) {
	  flp->wire_density[i] = (double *) calloc(flp->n_units, sizeof(double));
	  if (!flp->wire_density[i])
	  	fatal("memory allocation error\n");
	}

	/* copy connectivity (only for non-dead core blocks) */
	for(i=0; i < flp_desc->n_units-!!(wrap_l2); i++) {
	  strcpy(flp->units[i].name, flp_desc->units[i].name);
	  for (j=0; j < flp_desc->n_units-!!(wrap_l2); j++) {
	  	flp->wire_density[i][j] = flp_desc->wire_density[i][j];
	  }
	}

	/* name the dead blocks	*/
	for(count=0; count < n_dead; count++, i++)
		sprintf(flp->units[i].name, DEAD_PREFIX"%d", count);

	/* L2 connectivity info	*/
	if (wrap_l2)
		copy_l2_info(flp, i, flp_desc, flp_desc->n_units-1, flp_desc->n_units-1);

	return flp;
}

/*
 * note that if wrap_l2 is true, L2 is beyond the boundary in flp_desc
 * but flp contains it within its boundaries.
 */
void restore_dead_blocks(flp_t *flp, flp_desc_t *flp_desc,
						 int compacted, int wrap_l2,
						 int model_rim, int rim_blocks)
{
	int i, j, idx=0;
	/* remove L2 and rim blocks and restore the compacted blocks */
	if(model_rim)
		flp->n_units -= rim_blocks;
	if (wrap_l2)
		flp->n_units -= (L2_ARMS+1);
	flp->n_units += compacted;

	/* reinitialize the dead blocks	*/
	for(i=0; i < flp_desc->n_units-1; i++) {
		idx = flp_desc->n_units + i;
		sprintf(flp->units[idx].name, DEAD_PREFIX"%d", i);
		flp->units[idx].leftx = flp->units[idx].bottomy = 0;
		flp->units[idx].width = flp->units[idx].height = 0;
		for(j=0; j < flp->n_units; j++)
			flp->wire_density[idx][j] = flp->wire_density[j][idx] = 0;
	}
}

/* translate the floorplan to new origin (x,y)	*/
void flp_translate(flp_t *flp, double x, double y)
{
	int i;
	double minx = flp->units[0].leftx;
	double miny = flp->units[0].bottomy;

	for (i=1; i < flp->n_units; i++) {
		if (minx > flp->units[i].leftx)
			minx = flp->units[i].leftx;
		if (miny > flp->units[i].bottomy)
			miny = flp->units[i].bottomy;
	}
	for (i=0; i < flp->n_units; i++) {
		flp->units[i].leftx += (x - minx);
		flp->units[i].bottomy += (y - miny);
	}
}

/* scale the floorplan by a factor 'factor'	*/
void flp_scale(flp_t *flp, double factor)
{
	int i;
	double minx = flp->units[0].leftx;
	double miny = flp->units[0].bottomy;

	for (i=1; i < flp->n_units; i++) {
		if (minx > flp->units[i].leftx)
			minx = flp->units[i].leftx;
		if (miny > flp->units[i].bottomy)
			miny = flp->units[i].bottomy;
	}
	for(i=0; i < flp->n_units; i++) {
		flp->units[i].leftx = (flp->units[i].leftx - minx) * factor + minx;
		flp->units[i].bottomy = (flp->units[i].bottomy - miny) * factor + miny;
		flp->units[i].width *= factor;
		flp->units[i].height *= factor;
	}
}

/*
 * change the orientation of the floorplan by
 * rotating and/or flipping. the target orientation
 * is specified in 'target'. 'width', 'height', 'xorig'
 * and 'yorig' are those of 'flp' respectively.
 */
void flp_change_orient(flp_t *flp, double xorig, double yorig,
					   double width, double height, orient_t target)
{
	int i;

	for(i=0; i < flp->n_units; i++) {
		double leftx, bottomy, rightx, topy;
		/* all co-ordinate calculations are
		 * done assuming (0,0) as the center.
		 * so, shift accordingly
		 */
		leftx = flp->units[i].leftx  - (xorig + width / 2.0);
		bottomy = flp->units[i].bottomy - (yorig + height / 2.0);
		rightx = leftx + flp->units[i].width;
		topy = bottomy + flp->units[i].height;
		/* when changing orientation, leftx and
		 * bottomy of a rectangle could change
		 * to one of the other three corners.
		 * also, signs of the co-ordinates
		 * change according to the rotation
		 * or reflection. Further x & y are
		 * swapped for rotations that are
		 * odd multiples of 90 degrees
		 */
		switch(target) {
			case ROT_0:
					flp->units[i].leftx = leftx;
					flp->units[i].bottomy = bottomy;
					break;
			case ROT_90:
					flp->units[i].leftx = -topy;
					flp->units[i].bottomy = leftx;
					swap_dval(&(flp->units[i].width), &(flp->units[i].height));
					break;
			case ROT_180:
					flp->units[i].leftx = -rightx;
					flp->units[i].bottomy = -topy;
					break;
			case ROT_270:
					flp->units[i].leftx = bottomy;
					flp->units[i].bottomy = -rightx;
					swap_dval(&(flp->units[i].width), &(flp->units[i].height));
					break;
			case FLIP_0:
					flp->units[i].leftx = -rightx;
					flp->units[i].bottomy = bottomy;
					break;
			case FLIP_90:
					flp->units[i].leftx = bottomy;
					flp->units[i].bottomy = leftx;
					swap_dval(&(flp->units[i].width), &(flp->units[i].height));
					break;
			case FLIP_180:
					flp->units[i].leftx = leftx;
					flp->units[i].bottomy = -topy;
					break;
			case FLIP_270:
					flp->units[i].leftx = -topy;
					flp->units[i].bottomy = -rightx;
					swap_dval(&(flp->units[i].width), &(flp->units[i].height));
					break;
			default:
					fatal("unknown orientation\n");
					break;
		}
		/* translate back to original origin	*/
		flp->units[i].leftx += (xorig + width / 2.0);
		flp->units[i].bottomy += (yorig + height / 2.0);
	}
}

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
 * HotSpot. 'map', if non-NULL, is an output parameter to store
 * the 2-d array allocated by the function.
 */
flp_t *flp_create_grid(flp_t *flp, int ***map)
{
	double x[MAX_UNITS], y[MAX_UNITS];
	int i, j, n, xsize=0, ysize=0, count=0, found, **ptr;
	flp_t *grid;

	/* sort the units' boundary co-ordinates	*/
	for(i=0; i < flp->n_units; i++) {
		double r, t;
		r = flp->units[i].leftx + flp->units[i].width;
		t = flp->units[i].bottomy + flp->units[i].height;
		if(bsearch_insert_double(x, xsize, flp->units[i].leftx))
			xsize++;
		if(bsearch_insert_double(y, ysize, flp->units[i].bottomy))
			ysize++;
		if(bsearch_insert_double(x, xsize, r))
			xsize++;
		if(bsearch_insert_double(y, ysize, t))
			ysize++;
	}

	/*
	 * the grid formed by the lines from x and y arrays
	 * is our desired floorplan. allocate memory for it
	 */
	grid = (flp_t *) calloc (1, sizeof(flp_t));
	if(!grid)
		fatal("memory allocation error\n");
	grid->n_units = (xsize-1) * (ysize-1);
	grid->units = (unit_t *) calloc (grid->n_units, sizeof(unit_t));
	grid->wire_density = (double **) calloc(grid->n_units, sizeof(double *));
	if (!grid->units || !grid->wire_density)
		fatal("memory allocation error\n");
	for (i=0; i < grid->n_units; i++) {
	  grid->wire_density[i] = (double *) calloc(grid->n_units, sizeof(double));
	  if (!grid->wire_density[i])
	  	fatal("memory allocation error\n");
	}
	/* mapping between blocks of 'flp' to those of 'grid'	*/
	ptr = (int **) calloc(flp->n_units, sizeof(int *));
	if (!ptr)
		fatal("memory allocation error\n");
	/*
	 * ptr is a 2-d array with each row of possibly different
	 * length. the size of each row is stored in its first element.
	 * here, it is basically the mapping between 'flp' to 'grid'
	 * i.e., for each flp->unit, it stores the set of grid->units
	 * it maps to.
	 */
	for(i=0; i < flp->n_units; i++) {
		ptr[i] = (int *) calloc(grid->n_units+1, sizeof(int));
		if(!ptr[i])
	  		fatal("memory allocation error\n");
	}

	/*
	 * now populate the 'grid' blocks and map the blocks
	 * from 'flp' to 'grid'. for each block, identify the
	 * intervening lines that chop it into grid cells and
	 * assign the names of those cells from that of the
	 * block
	 */
	for(i=0; i < flp->n_units; i++) {
		double *xstart, *xend, *ystart, *yend;
		double *ptr1, *ptr2;
		int grid_num=0;
		if (!bsearch_double(x, xsize, flp->units[i].leftx, &xstart))
			fatal("invalid sorted arrays\n");
		if (!bsearch_double(x, xsize, flp->units[i].leftx+flp->units[i].width, &xend))
			fatal("invalid sorted arrays\n");
		if (!bsearch_double(y, ysize, flp->units[i].bottomy, &ystart))
			fatal("invalid sorted arrays\n");
		if (!bsearch_double(y, ysize, flp->units[i].bottomy+flp->units[i].height, &yend))
			fatal("invalid sorted arrays\n");
		for(ptr1 = xstart; ptr1 < xend; ptr1++)
			for(ptr2 = ystart; ptr2 < yend; ptr2++) {
				/* add this grid block if it has not been added already	*/
				for(n=0, found=FALSE; n < count; n++) {
					if (grid->units[n].leftx == ptr1[0] && grid->units[n].bottomy == ptr2[0]) {
						found = TRUE;
						break;
					}
				}
				if(!found) {
					sprintf(grid->units[count].name, "%s_%d", flp->units[i].name, grid_num);
					grid->units[count].leftx = ptr1[0];
					grid->units[count].bottomy = ptr2[0];
					grid->units[count].width = ptr1[1]-ptr1[0];
					grid->units[count].height = ptr2[1]-ptr2[0];
					/* map between position in 'flp' to that in 'grid'	*/
					ptr[i][++ptr[i][0]] = count;
					grid_num++;
          count++;
        }
      }
  }


	/* sanity check	*/
	if(count != (xsize-1) * (ysize-1))
		fatal("mismatch in the no. of units\n");

	/* fill-in the wire densities	*/
	for(i=0; i < flp->n_units; i++)
		for(j=0; j < flp->n_units; j++) {
			int p, q;
			for(p=1; p <= ptr[i][0]; p++)
				for(q=1; q <= ptr[j][0]; q++)
					grid->wire_density[ptr[i][p]][ptr[j][q]] = flp->wire_density[i][j];
		}

	/* output the map	*/
	if (map)
		(*map) = ptr;
	else
		free_blkgrid_map(flp, ptr);

	return grid;
}

/* free the map allocated by flp_create_grid	*/
void free_blkgrid_map(flp_t *flp, int **map)
{
	int i;

	for(i=0; i < flp->n_units; i++)
		free(map[i]);
	free(map);
}

/* translate power numbers to the grid created by flp_create_grid	*/
void xlate_power_blkgrid(flp_t *flp, flp_t *grid, \
						 double *bpower, double *gpower, int **map)
{
	int i, p;

	for(i=0; i < flp->n_units; i++)
		for(p=1; p <= map[i][0]; p++)
			/* retain the power density	*/
			gpower[map[i][p]] = bpower[i] / (flp->units[i].width * flp->units[i].height) *\
								grid->units[map[i][p]].width * grid->units[map[i][p]].height;
}

/*
 * wrap the L2 around this floorplan. L2's area information
 * is obtained from flp_desc. memory for L2 and its arms has
 * already been allocated in the flp. note that flp & flp_desc
 * have L2 hidden beyond the boundary at this point
 */
void flp_wrap_l2(flp_t *flp, flp_desc_t *flp_desc)
{
	/*
	 * x is the width of the L2 arms
	 * y is the height of the bottom portion
	 */
	double x, y, core_width, core_height, total_side, core_area, l2_area;
	unit_t *l2, *l2_left, *l2_right;

	/* find L2 dimensions so that the total chip becomes a square	*/
	core_area = get_total_area(flp);
	core_width = get_total_width(flp);
	core_height = get_total_height(flp);
	/* flp_desc has L2 hidden beyond the boundary	*/
	l2_area = flp_desc->units[flp_desc->n_units].area;
	total_side = sqrt(core_area + l2_area);
	/*
	 * width of the total chip after L2 wrapping is equal to
	 * the width of the core plus the width of the two arms
	 */
	x = (total_side - core_width) / 2.0;
	y = total_side - core_height;
	/*
	 * we are trying to solve the equation
	 * (2*x+core_width) * (y+core_height)
	 * = l2_area + core_area
	 * for x and y. it is possible that the values
	 * turnout to be negative if we restrict the
	 * total chip to be a square. in that case,
	 * theoretically, any value of x in the range
	 * (0, l2_area/(2*core_height)) and the
	 * corresponding value of y or any value of y
	 * in the range (0, l2_area/core_width) and the
	 * corresponding value of x would be a solution
	 * we look for a solution with a reasonable
	 * aspect ratio. i.e., we constrain kx = y (or
	 * ky = x  depending on the aspect ratio of the
	 * core) where k = WRAP_L2_RATIO. solving the equation
	 * with this constraint, we get the following
	 */
	if ( x <= 0 || y <= 0.0) {
		double sum;
		if (core_width >= core_height) {
			sum = WRAP_L2_RATIO * core_width + 2 * core_height;
			x = (sqrt(sum*sum + 8*WRAP_L2_RATIO*l2_area) - sum) / (4*WRAP_L2_RATIO);
			y = WRAP_L2_RATIO * x;
		} else {
			sum = core_width + 2 * WRAP_L2_RATIO * core_height;
			y = (sqrt(sum*sum + 8*WRAP_L2_RATIO*l2_area) - sum) / (4*WRAP_L2_RATIO);
			x = WRAP_L2_RATIO * y;
		}
		total_side = 2 * x + core_width;
	}

	/* fix the positions of core blocks	*/
	flp_translate(flp, x, y);

	/* restore the L2 blocks	*/
	flp->n_units += (L2_ARMS+1);
	/* copy L2 info again from flp_desc but from beyond the boundary	*/
	copy_l2_info(flp, flp->n_units-L2_ARMS-1, flp_desc,
				 flp_desc->n_units, flp_desc->n_units);

	/* fix the positions of the L2  blocks. connectivity
	 * information has already been fixed (in flp_placeholder).
	 * bottom L2 block - (leftx, bottomy) is already (0,0)
	 */
	l2 = &flp->units[flp->n_units-1-L2_ARMS];
	l2->width = total_side;
	l2->height = y;
	l2->leftx = l2->bottomy = 0;

	/* left L2 arm */
	l2_left = &flp->units[flp->n_units-L2_ARMS+L2_LEFT];
	l2_left->width = x;
	l2_left->height = core_height;
	l2_left->leftx = 0;
	l2_left->bottomy = y;

	/* right L2 arm */
	l2_right = &flp->units[flp->n_units-L2_ARMS+L2_RIGHT];
	l2_right->width = x;
	l2_right->height = core_height;
	l2_right->leftx = x + core_width;
	l2_right->bottomy = y;
}

/*
 * wrap the rim strips around. each edge has rim blocks
 * equal to the number of blocks abutting that edge. at
 * the four corners, the rim blocks are extended by the
 * rim thickness in a clockwise fashion
 */
int flp_wrap_rim(flp_t *flp, double rim_thickness)
{
	double width, height;
	int i, j = 0, k, n = flp->n_units;
	unit_t *unit;

	width = get_total_width(flp) + 2 * rim_thickness;
	height = get_total_height(flp) + 2 * rim_thickness;
	flp_translate(flp, rim_thickness, rim_thickness);

	for (i = 0; i < n; i++) {
		/* shortcut	*/
		unit = &flp->units[i];

		/* block is on the western border	*/
		if (eq(unit->leftx, rim_thickness)) {
			sprintf(flp->units[n+j].name, "%s_%s",
					RIM_LEFT_STR, unit->name);
			flp->units[n+j].width = rim_thickness;
			flp->units[n+j].height = unit->height;
			flp->units[n+j].leftx = 0;
			flp->units[n+j].bottomy = unit->bottomy;
			/* northwest corner	*/
			if (eq(unit->bottomy + unit->height, height-rim_thickness))
				flp->units[n+j].height += rim_thickness;
			j++;
		}

		/* block is on the eastern border	*/
		if (eq(unit->leftx + unit->width, width-rim_thickness)) {
			sprintf(flp->units[n+j].name, "%s_%s",
					RIM_RIGHT_STR, unit->name);
			flp->units[n+j].width = rim_thickness;
			flp->units[n+j].height = unit->height;
			flp->units[n+j].leftx = unit->leftx + unit->width;
			flp->units[n+j].bottomy = unit->bottomy;
			/* southeast corner	*/
			if (eq(unit->bottomy, rim_thickness)) {
				flp->units[n+j].height += rim_thickness;
				flp->units[n+j].bottomy = 0;
			}
			j++;
		}

		/* block is on the northern border 	*/
		if (eq(unit->bottomy + unit->height, height-rim_thickness)) {
			sprintf(flp->units[n+j].name, "%s_%s",
					RIM_TOP_STR, unit->name);
			flp->units[n+j].width = unit->width;
			flp->units[n+j].height = rim_thickness;
			flp->units[n+j].leftx = unit->leftx;
			flp->units[n+j].bottomy = unit->bottomy + unit->height;
			/* northeast corner	*/
			if (eq(unit->leftx + unit->width, width-rim_thickness))
				flp->units[n+j].width += rim_thickness;
			j++;
		}

		/* block is on the southern border	*/
		if (eq(unit->bottomy, rim_thickness)) {
			sprintf(flp->units[n+j].name, "%s_%s",
					RIM_BOTTOM_STR, unit->name);
			flp->units[n+j].width = unit->width;
			flp->units[n+j].height = rim_thickness;
			flp->units[n+j].leftx = unit->leftx;
			flp->units[n+j].bottomy = 0;
			/* southwest corner	*/
			if (eq(unit->leftx, rim_thickness)) {
				flp->units[n+j].width += rim_thickness;
				flp->units[n+j].leftx = 0;
			}
			j++;
		}
	}

	flp->n_units += j;

	/* update all the rim wire densities */
	for(i=n; i < n+j; i++)
		for(k=0; k <= i; k++)
			flp->wire_density[i][k] = flp->wire_density[k][i] = 0;

	return j;
}

/*
 * floorplanning using simulated annealing.
 * precondition: flp is a pre-allocated placeholder.
 * returns the number of compacted blocks in the selected
 * floorplan
 */
int floorplan(flp_t *flp, flp_desc_t *flp_desc,
			  RC_model_t *model, double *power)
{
	NPE_t *expr, *next, *best;	/* Normalized Polish Expressions */
	tree_node_stack_t *stack;	/* for NPE evaluation	*/
	tree_node_t *root;			/* shape curve tree	*/
	double cost, new_cost, best_cost, sum_cost, T, Tcold;
	int i, steps, downs, n, rejects, compacted, rim_blocks = 0;
	int original_n = flp->n_units;
	int wrap_l2;

	/* to maintain the order of power values during
	 * the compaction/shifting around of blocks
	 */
	double *tpower = hotspot_vector(model);

	/* shortcut	*/
	flp_config_t cfg = flp_desc->config;

	/*
	 * make the rim strips disappear for slicing tree
	 * purposes. can be restored at the end
	 */
	if (cfg.model_rim)
		flp->n_units = (flp->n_units - 2) / 3;

	/* wrap L2 around?	*/
	wrap_l2 = FALSE;
	if (cfg.wrap_l2 &&
		!strcasecmp(flp_desc->units[flp_desc->n_units-1].name, cfg.l2_label)) {
		wrap_l2 = TRUE;
		/* make L2 disappear too */
		flp_desc->n_units--;
		flp->n_units -= (L2_ARMS+1);
	}

	/* initialization	*/
	expr = NPE_get_initial(flp_desc);
	stack = new_tree_node_stack();
	init_rand();

	/* convert NPE to flp	*/
	root = tree_from_NPE(flp_desc, stack, expr);
	/* compacts too small dead blocks	*/
	compacted = tree_to_flp(root, flp, TRUE, cfg.compact_ratio);
	/* update the tpower vector according to the compaction	*/
	trim_hotspot_vector(model, tpower, power, flp->n_units, compacted);
	free_tree(root);
	if(wrap_l2)
		flp_wrap_l2(flp, flp_desc);
	if(cfg.model_rim)
		rim_blocks = flp_wrap_rim(flp, cfg.rim_thickness);

	resize_thermal_model(model, flp->n_units);
	#if VERBOSE > 2
	print_flp(flp, TRUE);
	#endif
	cost = flp_evaluate_metric(flp, model, tpower, cfg.lambdaA, cfg.lambdaT, cfg.lambdaW);
	/* restore the compacted blocks	*/
	restore_dead_blocks(flp, flp_desc, compacted, wrap_l2, cfg.model_rim, rim_blocks);

	best = NPE_duplicate(expr);	/* best till now	*/
	best_cost = cost;

	/* simulated annealing	*/
	steps = 0;
	/* initial annealing temperature	*/
	T = -cfg.Davg / log(cfg.P0);
	/*
	 * final annealing temperature - we stop when there
	 * are fewer than (1-cfg.Rreject) accepts.
	 * of those accepts, assuming half are uphill moves,
	 * we want the temperature so that the probability
	 * of accepting uphill moves is as low as
	 * (1-cfg.Rreject)/2.
	 */
	Tcold = -cfg.Davg / log ((1.0 - cfg.Rreject) / 2.0);
	#if VERBOSE > 0
	fprintf(stdout, "initial cost: %g\tinitial T: %g\tfinal T: %g\n", cost, T, Tcold);
	#endif
	/*
	 * stop annealing if temperature has cooled down enough or
	 * max no. of iterations have been tried
	 */
	while (T >= Tcold && steps < cfg.Nmax) {
		/* shortcut	*/
		n = cfg.Kmoves * flp->n_units;
		i = downs = rejects = 0;
		sum_cost = 0;
		/* try enough total or downhill moves per T */
		while ((i < 2 * n) && (downs < n)) {
			next = make_random_move(expr);

			/* convert NPE to flp	*/
			root = tree_from_NPE(flp_desc, stack, next);
			compacted = tree_to_flp(root, flp, TRUE, cfg.compact_ratio);
			/* update the tpower vector according to the compaction	*/
			trim_hotspot_vector(model, tpower, power, flp->n_units, compacted);
			free_tree(root);
			if(wrap_l2)
				flp_wrap_l2(flp, flp_desc);
			if(cfg.model_rim)
				rim_blocks = flp_wrap_rim(flp, cfg.rim_thickness);

			resize_thermal_model(model, flp->n_units);
			#if VERBOSE > 2
			print_flp(flp, TRUE);
			#endif
			new_cost = flp_evaluate_metric(flp, model, tpower, cfg.lambdaA, cfg.lambdaT, cfg.lambdaW);
			restore_dead_blocks(flp, flp_desc, compacted, wrap_l2, cfg.model_rim, rim_blocks);

			#if VERBOSE > 1
			fprintf(stdout, "count: %d\tdowns: %d\tcost: %g\t",
					i, downs, new_cost);
			#endif

			/* move accepted?	*/
			if (new_cost < cost || 	/* downhill always accepted	*/
				/* boltzmann probability function	*/
			    rand_fraction() < exp(-(new_cost-cost)/T)) {

				free_NPE(expr);
				expr = next;

				/* downhill move	*/
				if (new_cost < cost) {
					downs++;
					/* found new best	*/
					if (new_cost < best_cost) {
						free_NPE(best);
						best = NPE_duplicate(expr);
						best_cost = new_cost;
					}
				}

				#if VERBOSE > 1
				fprintf(stdout, "accepted\n");
				#endif
				cost = new_cost;
				sum_cost += cost;
			} else {	/* rejected move	*/
				rejects++;
				free_NPE(next);
				#if VERBOSE > 1
				fprintf(stdout, "rejected\n");
				#endif
			}
			i++;
		}
		#if VERBOSE > 0
		fprintf(stdout, "step: %d\tT: %g\ttries: %d\taccepts: %d\trejects: %d\t",
				steps, T, i, (i-rejects), rejects);
		fprintf(stdout, "avg. cost: %g\tbest cost: %g\n",
		 		(i-rejects)?(sum_cost / (i-rejects)):sum_cost, best_cost);
		#endif

		/* stop annealing if there are too little accepts */
		if(((double)rejects/i) > cfg.Rreject)
			break;

		/* annealing schedule	*/
		T *= cfg.Rcool;
		steps++;
	}

	/* best floorplan found	*/
	root = tree_from_NPE(flp_desc, stack, best);
	#if VERBOSE > 0
	{
		int pos = min_area_pos(root->curve);
		print_tree_relevant(root, pos, flp_desc);
	}
	#endif
	compacted = tree_to_flp(root, flp, TRUE, cfg.compact_ratio);
	/* update the power vector according to the compaction	*/
	trim_hotspot_vector(model, power, power, flp->n_units, compacted);
	free_tree(root);
	/*  restore L2 and rim */
	if(wrap_l2) {
		flp_wrap_l2(flp, flp_desc);
		flp_desc->n_units++;
	}
	if(cfg.model_rim)
		rim_blocks = flp_wrap_rim(flp, cfg.rim_thickness);
	resize_thermal_model(model, flp->n_units);
	#if VERBOSE > 2
	print_flp(flp, TRUE);
	#endif

	free_NPE(expr);
	free_NPE(best);
	free_tree_node_stack(stack);
	free_dvector(tpower);

	/*
	 * return the number of blocks compacted finally
	 * so that any deallocator can take care of memory
	 * accordingly.
	 */
	return (original_n - flp->n_units);
}

/* functions duplicated from flp_desc.c */
/*
 * find the number of units from the
 * floorplan file
 */
int flp_count_units(FILE *fp)
{
    char str1[LINE_SIZE], str2[LINE_SIZE];
	char name[STR_SIZE];
	double leftx, bottomy, width, height;
	double cp, res;
	char *ptr;
    int count = 0;

	fseek(fp, 0, SEEK_SET);
	while(!feof(fp)) {
		fgets(str1, LINE_SIZE, fp);
		if (feof(fp))
			break;
		strcpy(str2, str1);

		/* ignore comments and empty lines	*/
		ptr = strtok(str1, " \r\t\n");
		if (!ptr || ptr[0] == '#')
			continue;

		/* functional block placement information	*/
		if (sscanf(str2, "%s%lf%lf%lf%lf%lf%lf", name, &leftx, &bottomy,
		  		   &width, &height, &cp, &res) > 4)
			count++;
	}
	return count;
}

flp_t *flp_alloc_init_mem(int count, int use_wire_density)
{
	int i;
	flp_t *flp;
	flp = (flp_t *) calloc (1, sizeof(flp_t));
	if(!flp)
		fatal("memory allocation error\n");
	flp->units = (unit_t *) calloc(count, sizeof(unit_t));
  if(!flp->units)
    fatal("memory allocation error\n");
  if(use_wire_density) {
    flp->wire_density = (double **) calloc(count, sizeof(double *));
	   if (!flp->wire_density)
		   fatal("memory allocation error\n");

   for (i=0; i < count; i++) {
     flp->wire_density[i] = (double *) calloc(count, sizeof(double));
 	  if (!flp->wire_density[i])
 	  	fatal("memory allocation error\n");
 	 }
  }
	flp->n_units = count;
	return flp;
}

/* populate block information	*/
void flp_populate_blks(flp_t *flp, FILE *fp)
{
	int i=0;
	char str[LINE_SIZE], copy[LINE_SIZE];
	char name1[STR_SIZE], name2[STR_SIZE];
	double width, height, leftx, bottomy, cp, res, temp;
	double wire_density;
	char *ptr;

	fseek(fp, 0, SEEK_SET);
	while(!feof(fp)) {		/* second pass	*/
		fgets(str, LINE_SIZE, fp);
		if (feof(fp))
			break;
		strcpy(copy, str);

		/* ignore comments and empty lines	*/
		ptr = strtok(str, " \r\t\n");
		if (!ptr || ptr[0] == '#')
			continue;
		cp = res = 0.0;
		if (sscanf(copy, "%s%lf%lf%lf%lf%lf%lf%lf", name1, &width, &height,
					   &leftx, &bottomy, &cp, &res, &temp) > 7)
			fatal("invalid floorplan file format\n");
    else if (sscanf(copy, "%s%lf%lf%lf%lf%lf%lf", name1, &width, &height,
					   &leftx, &bottomy, &cp, &res) == 7)
		{
				strcpy(flp->units[i].name, name1);
				flp->units[i].width = width;
				flp->units[i].height = height;
				flp->units[i].leftx = leftx;
				flp->units[i].bottomy = bottomy;
				flp->units[i].specificheat = cp; //BU_3D set specific heat
				flp->units[i].resistivity = res; //BU_3D set resistivity
				flp->units[i].hasRes = TRUE;
				flp->units[i].hasSh = TRUE;
				i++;
		}
		/* Default Option */
		else if(sscanf(copy, "%s%lf%lf%lf%lf%lf%lf", name1, &width, &height,
					   &leftx, &bottomy, &cp, &res) == 5)
		{
				strcpy(flp->units[i].name, name1);
				flp->units[i].width = width;
				flp->units[i].height = height;
				flp->units[i].leftx = leftx;
				flp->units[i].bottomy = bottomy;
				flp->units[i].specificheat = 0.0;	//BU_3D else - set specific heat to default value 0.0
				flp->units[i].resistivity = 0.0; 	//BU_3D set resistivity to default value 0.0
				flp->units[i].hasRes = FALSE; 		//BU_3D set boolean for has resistivity flag to false
				flp->units[i].hasSh = FALSE; 		//BU_3D set boolean for has specific heat flag to false
				i++;
		}
		/* skip connectivity info	*/
		else if(sscanf(copy, "%s%s%lf", name1, name2, &wire_density) != 3)
			fatal("invalid floorplan file format\n");
	}
	if (i != flp->n_units)
	  fatal("mismatch of number of units\n");
}

/* populate connectivity info	*/
void flp_populate_connects(flp_t *flp, FILE *fp)
{
	char str1[LINE_SIZE], str2[LINE_SIZE];
	char name1[STR_SIZE], name2[STR_SIZE];
	/* dummy fields	*/
	double f1, f2, f3, f4, f5, f6;
	double wire_density;
	char *ptr;
	int x, y, temp;

	/* initialize wire_density	*/
	for(x=0; x < flp->n_units; x++)
		for(y=0; y < flp->n_units; y++)
			flp->wire_density[x][y] = 0.0;

	fseek(fp, 0, SEEK_SET);
	while(!feof(fp)) {
		fgets(str1, LINE_SIZE, fp);
		if (feof(fp))
			break;
		strcpy(str2, str1);

		/* ignore comments and empty lines	*/
		ptr = strtok(str1, " \r\t\n");
		if (!ptr || ptr[0] == '#')
			continue;

		/* lines with unit positions	*/
		if (sscanf(str2, "%s%lf%lf%lf%lf%lf%lf", name1, &f1, &f2, &f3, &f4, &f5, &f6) == 7 ||
		  	/* flp_desc like lines. ignore them	*/
		  	sscanf(str2, "%s%lf%lf%lf%d", name1, &f1, &f2, &f3, &temp) == 5)
			continue;

		 /* lines with connectivity info	*/
		else if (sscanf(str2, "%s%s%lf", name1, name2, &wire_density) == 3) {
			x = get_blk_index(flp, name1);
			y = get_blk_index(flp, name2);

			if (x == y)
				fatal("block connected to itself?\n");

			if (!flp->wire_density[x][y] && !flp->wire_density[y][x])
				flp->wire_density[x][y] = flp->wire_density[y][x] = wire_density;
			else if((flp->wire_density[x][y] != flp->wire_density[y][x]) ||
			        (flp->wire_density[x][y] != wire_density)) {
				sprintf(str2, "wrong connectivity information for blocks %s and %s\n",
				        name1, name2);
				fatal(str2);
			}
		} else
		  	fatal("invalid floorplan file format\n");
	} /* end while	*/
}

flp_t *read_flp(char *file, int read_connects, int initialize_connects)
{
	char str[STR_SIZE];
	FILE *fp;
	flp_t *flp;
	int count, i, j;

	if (!strcasecmp(file, "stdin"))
		fp = stdin;
	else
		fp = fopen (file, "r");

	if (!fp) {
		sprintf(str, "error opening file %s\n", file);
		fatal(str);
	}

	/* 1st pass - find n_units	*/
	count = flp_count_units(fp);
	if(!count)
		fatal("no units specified in the floorplan file\n");

	/* allocate initial memory */
	flp = flp_alloc_init_mem(count, read_connects || initialize_connects);

	/* 2nd pass - populate block info	*/
	flp_populate_blks(flp, fp);

	/* 3rd pass - populate connectivity info    */
	if (read_connects)
		flp_populate_connects(flp, fp);
	/* older version - no connectivity	*/
	else if(initialize_connects)
    for (i=0; i < flp->n_units; i++)
			for (j=0; j < flp->n_units; j++)
				flp->wire_density[i][j] = 1.0;

	if(fp != stdin)
		fclose(fp);

	/* make sure the origin is (0,0)	*/
	flp_translate(flp, 0, 0);
	return flp;
}

void dump_flp(flp_t *flp, char *file, int dump_connects)
{
	char str[STR_SIZE];
	int i, j;
	FILE *fp;

	if (!strcasecmp(file, "stdout"))
		fp = stdout;
	else if (!strcasecmp(file, "stderr"))
		fp = stderr;
	else
		fp = fopen (file, "w");

	if (!fp) {
		sprintf(str, "error opening file %s\n", file);
		fatal(str);
	}
	/* functional unit placement info	*/
	for(i=0; i < flp->n_units; i++)  {
		fprintf(fp, "%s\t%.11f\t%.11f\t%.11f\t%.11f\n",
				flp->units[i].name, flp->units[i].width, flp->units[i].height,
				flp->units[i].leftx, flp->units[i].bottomy);
	}

	if (dump_connects) {
		fprintf(fp, "\n");
		/* connectivity information	*/
		for(i=1; i < flp->n_units; i++)
			for(j=0; j < i; j++)
				if (flp->wire_density[i][j])
					fprintf(fp, "%s\t%s\t%.3f\n", flp->units[i].name,
							flp->units[j].name, flp->wire_density[i][j]);
	}

	if(fp != stdout && fp != stderr)
		fclose(fp);
}

void free_flp(flp_t *flp, int compacted, int free_connects)
{
	int i;
  if(free_connects) {
  	for (i=0; i < flp->n_units + compacted; i++) {
  		free(flp->wire_density[i]);
  	}
    free(flp->wire_density);
  }
	free(flp->units);
	free(flp);
}

void print_flp_fig (flp_t *flp)
{
	int i;
	double leftx, bottomy, rightx, topy;

	fprintf(stdout, "FIG starts\n");
	for (i=0; i< flp->n_units; i++) {
		leftx = flp->units[i].leftx;
		bottomy = flp->units[i].bottomy;
		rightx = flp->units[i].leftx + flp->units[i].width;
		topy = flp->units[i].bottomy + flp->units[i].height;
		fprintf(stdout, "%.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f\n",
			    leftx, bottomy, leftx, topy, rightx, topy, rightx, bottomy,
				leftx, bottomy);
		fprintf(stdout, "%s\n", flp->units[i].name);
	}
	fprintf(stdout, "FIG ends\n");
}

/* debug print	*/
void print_flp (flp_t *flp, int print_connects)
{
	int i, j;

	fprintf(stdout, "printing floorplan information for %d blocks\n", flp->n_units);
	fprintf(stdout, "name\tarea\twidth\theight\tleftx\tbottomy\trightx\ttopy\n");
	for (i=0; i< flp->n_units; i++) {
		double area, width, height, leftx, bottomy, rightx, topy;
		char *name;
		name = flp->units[i].name;
		width = flp->units[i].width;
		height = flp->units[i].height;
		area = width * height;
		leftx = flp->units[i].leftx;
		bottomy = flp->units[i].bottomy;
		rightx = flp->units[i].leftx + flp->units[i].width;
		topy = flp->units[i].bottomy + flp->units[i].height;
		fprintf(stdout, "%s\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\n",
			    name, area, width, height, leftx, bottomy, rightx, topy);
	}
  if(print_connects) {
  	fprintf(stdout, "printing connections:\n");
  	for (i=0; i< flp->n_units; i++)
  		for (j=i+1; j < flp->n_units; j++)
  			if (flp->wire_density[i][j])
  				fprintf(stdout, "%s\t%s\t%lg\n", flp->units[i].name,
  						flp->units[j].name, flp->wire_density[i][j]);
  }
}

/* print the statistics about this floorplan.
 * note that connects_file is NULL if wire
 * information is already populated
 */
void print_flp_stats(flp_t *flp, RC_model_t *model,
					 char *l2_label, char *power_file,
					 char *connects_file)
{
	double core, total, occupied;	/* area	*/
	double width, height, aspect, total_w, total_h;
	double wire_metric;
	double peak, avg;		/* temperature	*/
	double *power, *temp;
	FILE *fp = NULL;
	char str[STR_SIZE];

	if (connects_file) {
		fp = fopen(connects_file, "r");
		if (!fp) {
			sprintf(str, "error opening file %s\n", connects_file);
			fatal(str);
		}
		flp_populate_connects(flp, fp);
	}

	power = hotspot_vector(model);
	temp = hotspot_vector(model);
	read_power(model, power, power_file);

	core = get_core_area(flp, l2_label);
	total = get_total_area(flp);
	total_w = get_total_width(flp);
	total_h = get_total_height(flp);
	occupied = get_core_occupied_area(flp, l2_label);
	width = get_core_width(flp, l2_label);
	height = get_core_height(flp, l2_label);
	aspect = (height > width) ? (height/width) : (width/height);
	wire_metric = get_wire_metric(flp);

	populate_R_model(model, flp);
	steady_state_temp(model, power, temp);
	peak = find_max_temp(model, temp);
	avg = find_avg_temp(model, temp);

	fprintf(stdout, "printing summary statistics about the floorplan\n");
	fprintf(stdout, "total area:\t%g\n", total);
	fprintf(stdout, "total width:\t%g\n", total_w);
	fprintf(stdout, "total height:\t%g\n", total_h);
	fprintf(stdout, "core area:\t%g\n", core);
	fprintf(stdout, "occupied area:\t%g\n", occupied);
	fprintf(stdout, "area utilization:\t%.3f\n", occupied / core * 100.0);
	fprintf(stdout, "core width:\t%g\n", width);
	fprintf(stdout, "core height:\t%g\n", height);
	fprintf(stdout, "core aspect ratio:\t%.3f\n", aspect);
	fprintf(stdout, "wire length metric:\t%.3f\n", wire_metric);
	fprintf(stdout, "peak temperature:\t%.3f\n", peak);
	fprintf(stdout, "avg temperature:\t%.3f\n", avg);

	free_dvector(power);
	free_dvector(temp);
	if (fp)
		fclose(fp);
}

int get_blk_index(flp_t *flp, char *name)
{
	int i;
	char msg[STR_SIZE];

	if (!flp)
		fatal("null pointer in get_blk_index\n");

	for (i = 0; i < flp->n_units; i++) {
    if (!strcasecmp(name, flp->units[i].name)) {
			return i;
		}
	}

	sprintf(msg, "block %s not found\n", name);
	fatal(msg);
	return -1;
}

int is_horiz_adj(flp_t *flp, int i, int j)
{
	double x1, x2, x3, x4;
	double y1, y2, y3, y4;

	if (i == j)
		return FALSE;

	x1 = flp->units[i].leftx;
	x2 = x1 + flp->units[i].width;
	x3 = flp->units[j].leftx;
	x4 = x3 + flp->units[j].width;

	y1 = flp->units[i].bottomy;
	y2 = y1 + flp->units[i].height;
	y3 = flp->units[j].bottomy;
	y4 = y3 + flp->units[j].height;

	/* diagonally adjacent => not adjacent */
	if (eq(x2,x3) && eq(y2,y3))
		return FALSE;
	if (eq(x1,x4) && eq(y1,y4))
		return FALSE;
	if (eq(x2,x3) && eq(y1,y4))
		return FALSE;
	if (eq(x1,x4) && eq(y2,y3))
		return FALSE;

	if (eq(x1,x4) || eq(x2,x3))
		if ((y3 >= y1 && y3 <= y2) || (y4 >= y1 && y4 <= y2) ||
		    (y1 >= y3 && y1 <= y4) || (y2 >= y3 && y2 <= y4))
			return TRUE;

	return FALSE;
}

int is_vert_adj (flp_t *flp, int i, int j)
{
	double x1, x2, x3, x4;
	double y1, y2, y3, y4;

	if (i == j)
		return FALSE;

	x1 = flp->units[i].leftx;
	x2 = x1 + flp->units[i].width;
	x3 = flp->units[j].leftx;
	x4 = x3 + flp->units[j].width;

	y1 = flp->units[i].bottomy;
	y2 = y1 + flp->units[i].height;
	y3 = flp->units[j].bottomy;
	y4 = y3 + flp->units[j].height;

	/* diagonally adjacent => not adjacent */
	if (eq(x2,x3) && eq(y2,y3))
		return FALSE;
	if (eq(x1,x4) && eq(y1,y4))
		return FALSE;
	if (eq(x2,x3) && eq(y1,y4))
		return FALSE;
	if (eq(x1,x4) && eq(y2,y3))
		return FALSE;

	if (eq(y1,y4) || eq(y2,y3))
		if ((x3 >= x1 && x3 <= x2) || (x4 >= x1 && x4 <= x2) ||
		    (x1 >= x3 && x1 <= x4) || (x2 >= x3 && x2 <= x4))
			return TRUE;

	return FALSE;
}

double get_shared_len(flp_t *flp, int i, int j)
{
	double p11, p12, p21, p22;
	p11 = p12 = p21 = p22 = 0.0;

	if (i==j)
		return FALSE;

	if (is_horiz_adj(flp, i, j)) {
		p11 = flp->units[i].bottomy;
		p12 = p11 + flp->units[i].height;
		p21 = flp->units[j].bottomy;
		p22 = p21 + flp->units[j].height;
	}

	if (is_vert_adj(flp, i, j)) {
		p11 = flp->units[i].leftx;
		p12 = p11 + flp->units[i].width;
		p21 = flp->units[j].leftx;
		p22 = p21 + flp->units[j].width;
	}

	return (MIN(p12, p22) - MAX(p11, p21));
}

double get_total_width(flp_t *flp)
{
	int i;
	double min_x = flp->units[0].leftx;
	double max_x = flp->units[0].leftx + flp->units[0].width;

	for (i=1; i < flp->n_units; i++) {
		if (flp->units[i].leftx < min_x)
			min_x = flp->units[i].leftx;
		if (flp->units[i].leftx + flp->units[i].width > max_x)
			max_x = flp->units[i].leftx + flp->units[i].width;
	}

	return (max_x - min_x);
}

double get_total_height(flp_t *flp)
{
	int i;
	double min_y = flp->units[0].bottomy;
	double max_y = flp->units[0].bottomy + flp->units[0].height;

	for (i=1; i < flp->n_units; i++) {
		if (flp->units[i].bottomy < min_y)
			min_y = flp->units[i].bottomy;
		if (flp->units[i].bottomy + flp->units[i].height > max_y)
			max_y = flp->units[i].bottomy + flp->units[i].height;
	}

	return (max_y - min_y);
}

double get_minx(flp_t *flp)
{
	int i;
	double min_x = flp->units[0].leftx;

	for (i=1; i < flp->n_units; i++)
		if (flp->units[i].leftx < min_x)
			min_x = flp->units[i].leftx;

	return min_x;
}

double get_miny(flp_t *flp)
{
	int i;
	double min_y = flp->units[0].bottomy;

	for (i=1; i < flp->n_units; i++)
		if (flp->units[i].bottomy < min_y)
			min_y = flp->units[i].bottomy;

	return min_y;
}

/* precondition: L2 should have been wrapped around	*/
double get_core_width(flp_t *flp, char *l2_label)
{
	int i;
	double min_x = LARGENUM;
	double max_x = -LARGENUM;

	for (i=0; i < flp->n_units; i++) {
		/* core is that part of the chip excluding the l2 and rim	*/
		if (strstr(flp->units[i].name, l2_label) != flp->units[i].name &&
			strstr(flp->units[i].name, RIM_PREFIX) != flp->units[i].name) {
			if (flp->units[i].leftx < min_x)
				min_x = flp->units[i].leftx;
			if (flp->units[i].leftx + flp->units[i].width > max_x)
				max_x = flp->units[i].leftx + flp->units[i].width;
		}
	}

	return (max_x - min_x);
}

/* precondition: L2 should have been wrapped around	*/
double get_core_height(flp_t *flp, char *l2_label)
{
	int i;
	double min_y = LARGENUM;
	double max_y = -LARGENUM;

	for (i=0; i < flp->n_units; i++) {
		/* core is that part of the chip excluding the l2 and rim	*/
		if (strstr(flp->units[i].name, l2_label) != flp->units[i].name &&
			strstr(flp->units[i].name, RIM_PREFIX) != flp->units[i].name) {
			if (flp->units[i].bottomy < min_y)
				min_y = flp->units[i].bottomy;
			if (flp->units[i].bottomy + flp->units[i].height > max_y)
				max_y = flp->units[i].bottomy + flp->units[i].height;
		}
	}

	return (max_y - min_y);
}

double get_total_area(flp_t *flp)
{
	int i;
	double area = 0.0;
	for(i=0; i < flp->n_units; i++)
		area += flp->units[i].width * flp->units[i].height;
	return area;
}

double get_core_area(flp_t *flp, char *l2_label)
{
	int i;
	double area = 0.0;
	for(i=0; i < flp->n_units; i++)
		if (strstr(flp->units[i].name, l2_label) != flp->units[i].name &&
			strstr(flp->units[i].name, RIM_PREFIX) != flp->units[i].name)
			area += flp->units[i].width * flp->units[i].height;
	return area;
}

/* excluding the dead blocks	*/
double get_core_occupied_area(flp_t *flp, char *l2_label)
{
	int i, num;
	double dead_area = 0.0;
	for(i=0; i < flp->n_units; i++) {
		/*
		 * there can be a max of n-1 dead blocks where n is the
		 * number of non-dead blocks (since each cut, vertical
		 * or horizontal, can correspond to a maximum of one
		 * dead block
		 */
		if ((sscanf(flp->units[i].name, DEAD_PREFIX"%d", &num) == 1) &&
			(num < (flp->n_units-1) / 2))
			dead_area += flp->units[i].width * flp->units[i].height;
	}
	return get_core_area(flp, l2_label) - dead_area;
}

double get_wire_metric(flp_t *flp)
{
	int i, j;
	double w = 0.0, dist;

	for (i=0; i < flp->n_units; i++)
		for (j=0; j < flp->n_units; j++)
			if (flp->wire_density[i][j]) {
				dist = get_manhattan_dist(flp, i, j);
				w += flp->wire_density[i][j] * dist;
			}
	return w;
}

double get_manhattan_dist(flp_t *flp, int i, int j)
{
	double x1 = flp->units[i].leftx + flp->units[i].width / 2.0;
	double y1 = flp->units[i].bottomy + flp->units[i].height / 2.0;
	double x2 = flp->units[j].leftx + flp->units[j].width / 2.0;
	double y2 = flp->units[j].bottomy + flp->units[j].height / 2.0;
	return (fabs(x2-x1) + fabs(y2-y1));
}

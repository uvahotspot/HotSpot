#include <stdio.h>
#ifdef _MSC_VER
#define strcasecmp    _stricmp
#define strncasecmp   _strnicmp
#else
#include <strings.h>
#endif
#include <string.h>
#include <stdlib.h>

#include "flp.h"
#include "shape.h"
#include "util.h"

/* get unit index from its name */
int desc_get_blk_index(flp_desc_t *flp_desc, char *name)
{
	int i;
	char msg[STR_SIZE];

	if (!flp_desc)
		fatal("null pointer in get_desc_blk_index\n");

	for (i = 0; i < flp_desc->n_units; i++) {
		if (!strcasecmp(name, flp_desc->units[i].name)) {
			return i;
		}
	}

	sprintf(msg, "block %s not found\n", name);
	fatal(msg);
	return -1;
}


/* 
 * find the number of units from the 
 * floorplan description file
 */
int desc_count_units(FILE *fp)
{
    char str1[LINE_SIZE], str2[LINE_SIZE];
	char name[STR_SIZE];
	double area, min, max;
	int rotable;
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

		/* lines describing functional block size and aspect ratio	*/
		if (sscanf(str2, "%s%lf%lf%lf%d", name, &area, &min, &max, &rotable) == 5)
			count++;
	}
	return count;
}

flp_desc_t *desc_alloc_init_mem(int count, flp_config_t *config)
{
	int i;
	flp_desc_t *flp_desc;
	flp_desc = (flp_desc_t *) calloc (1, sizeof(flp_desc_t));
	if(!flp_desc)
		fatal("memory allocation error\n");
	flp_desc->units = (unplaced_t *) calloc(count, sizeof(unplaced_t));
	flp_desc->wire_density = (double **) calloc(count, sizeof(double *));
	if (!flp_desc->units || !flp_desc->wire_density)
		fatal("memory allocation error\n");
	flp_desc->n_units = count;
	flp_desc->config = *config;

	for (i=0; i < count; i++) {
	  flp_desc->wire_density[i] = (double *) calloc(count, sizeof(double));
	  if (!flp_desc->wire_density[i])
	  	fatal("memory allocation error\n");
	}
	return flp_desc;
}

/* populate block information	*/
void desc_populate_blks(flp_desc_t *flp_desc, FILE *fp)
{
	int i=0;
	char str1[LINE_SIZE], str2[LINE_SIZE]; 
	char name1[STR_SIZE], name2[STR_SIZE];
	double area, min, max, wire_density;
	int rotable;
	char *ptr;
	int wrap_l2 = FALSE;

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

		/* lines describing functional block size and aspect ratio	*/
		if (sscanf(str2, "%s%lf%lf%lf%d", name1, &area, &min, &max, &rotable) == 5) {
			if (min > max)
				fatal("minimum aspect ratio greater than maximum\n");
			if (min <= 0 || max <= 0 || area <= 0)
				fatal("invalid number in floorplan description\n");
		  	/* L2 wrap around	*/
		  	if (flp_desc->config.wrap_l2 && !strcasecmp(name1, flp_desc->config.l2_label)) {
				wrap_l2 = TRUE;
                /* min, max, rotable and shape curve have no meaning here   */
                strcpy(flp_desc->units[flp_desc->n_units - 1].name, name1);
                flp_desc->units[flp_desc->n_units - 1].area = area;
			} else {
		    	strcpy(flp_desc->units[i].name, name1);
				flp_desc->units[i].area = area;
				flp_desc->units[i].min_aspect = min;
				flp_desc->units[i].max_aspect = max;
				flp_desc->units[i].rotable = rotable;
				flp_desc->units[i].shape = shape_from_aspect(area, min, max, rotable, 
															 flp_desc->config.n_orients);
		    	i++;
			}
		} else if (sscanf(str2, "%s%s%lf", name1, name2, &wire_density) != 3)
				fatal("invalid floorplan description file format\n");
	}

	if((wrap_l2 && i != flp_desc->n_units - 1) ||
	  (!wrap_l2 && i != flp_desc->n_units))
		fatal("mismatch of number of units\n");
}

/* populate connectivity info	*/
void desc_populate_connects(flp_desc_t *flp_desc, FILE *fp)
{
	char str1[LINE_SIZE], str2[LINE_SIZE]; 
	char name1[STR_SIZE], name2[STR_SIZE];
	double area, min, max, wire_density;
	int rotable;
	char *ptr;

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
			
		/* ignore lines describing functional block size and aspect ratio	*/
		if (sscanf(str2, "%s%lf%lf%lf%d", name1, &area, &min, &max, &rotable) == 5)
			continue;

		/* lines with connectivity info	*/
		else if (sscanf(str2, "%s%s%lf", name1, name2, &wire_density) == 3) {
			int x, y;
			x = desc_get_blk_index(flp_desc, name1);
			y = desc_get_blk_index(flp_desc, name2);

			if (x == y)
				fatal("block connected to itself?\n");

			if (!flp_desc->wire_density[x][y] && !flp_desc->wire_density[y][x])
				flp_desc->wire_density[x][y] = flp_desc->wire_density[y][x] = wire_density;
			else if((flp_desc->wire_density[x][y] != flp_desc->wire_density[y][x]) ||
			        (flp_desc->wire_density[x][y] != wire_density)) {
				sprintf(str2, "wrong connectivity information for blocks %s and %s\n", 
						name1, name2);
				fatal(str2);
			}
		} else 
		  	fatal("invalid floorplan description file format\n");
	} /* end while	*/
}

/* read floorplan description and allocate memory	*/
flp_desc_t *read_flp_desc(char *file, flp_config_t *config)
{
	char str[STR_SIZE];
	FILE *fp = fopen(file, "r");
	flp_desc_t *flp_desc;
	int count;

	if (!fp) {
		sprintf(str, "error opening file %s\n", file);
		fatal(str);
	}

	/* 1st pass - find n_units	*/
	count = desc_count_units(fp);
	if(!count)
		fatal("no units specified in the floorplan description file\n");

	/* allocate initial memory */
	flp_desc = desc_alloc_init_mem(count, config);

	/* 2nd pass - populate block info	*/
	desc_populate_blks(flp_desc, fp);

	/* 3rd pass - populate connectivity info	*/
	desc_populate_connects(flp_desc, fp);

	fclose(fp);	
	return flp_desc;
}

void free_flp_desc(flp_desc_t *flp_desc)
{
	int i;
	for (i=0; i < flp_desc->n_units; i++) {
		/* wrapped L2 doesn't have a shape curve	*/
		if (flp_desc->units[i].shape)
			free_shape(flp_desc->units[i].shape);
		free(flp_desc->wire_density[i]);
	}
	free(flp_desc->units);
	free(flp_desc->wire_density);
	free(flp_desc);
}

/* debug print	*/
void print_unplaced(unplaced_t *unit)
{
	fprintf(stdout, "printing unit info for %s\n", unit->name);
	fprintf(stdout, "area=%.5f\tmin=%.5f\tmax=%.5f\trotable=%d\n",
	        unit->area, unit->min_aspect, unit->max_aspect, unit->rotable);
	print_shape(unit->shape);
	fprintf(stdout, "\n");
}

void desc_print_wire_density(flp_desc_t *flp_desc)
{
	int i, j;

	fprintf(stdout, "printing wire density matrix\n\t");

	/* header row	*/
	for(i=0; i < flp_desc->n_units; i++)
		fprintf(stdout, "%s\t", flp_desc->units[i].name);
	fprintf(stdout, "\n");

	for(i=0; i < flp_desc->n_units; i++) {
		fprintf(stdout, "%s\t", flp_desc->units[i].name);
		for(j=0; j < flp_desc->n_units; j++)
			fprintf(stdout, "%.3f\t", flp_desc->wire_density[i][j]);
		fprintf(stdout, "\n");
	}
	fprintf(stdout, "\n");
}


/* debug print	*/
void print_flp_desc(flp_desc_t *flp_desc)
{
	int i;

	for(i=0;  i < flp_desc->n_units; i++)
		print_unplaced(&flp_desc->units[i]);

	desc_print_wire_density(flp_desc);
	fprintf(stdout, "\n");
}


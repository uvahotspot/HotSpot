#ifndef __HOTFLOORPLAN_H_
#define __HOTFLOORPLAN_H_

#include "util.h"

/* global configuration parameters for HotFloorplan	*/
typedef struct global_config_t_st
{
	/* floorplan description input file	*/
	char flp_desc[STR_SIZE];
	/* power input file	*/
	char power_in[STR_SIZE];
  /* input material properties file */
  char materials_file[STR_SIZE];
	/* floorplan output file */
	char flp_out[STR_SIZE];
	/* input configuration parameters from file	*/
	char config[STR_SIZE];
	/* output configuration parameters to file	*/
	char dump_config[STR_SIZE];
}global_config_t;

/*
 * parse a table of name-value string pairs and add the configuration
 * parameters to 'config'
 */
void global_config_from_strs(global_config_t *config, str_pair *table, int size);
/*
 * convert config into a table of name-value pairs. returns the no.
 * of parameters converted
 */
int global_config_to_strs(global_config_t *config, str_pair *table, int max_entries);

#endif

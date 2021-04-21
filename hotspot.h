#ifndef __HOTSPOT_H_
#define __HOTSPOT_H_

#include "util.h"

/* global configuration parameters for HotSpot	*/
typedef struct global_config_t_st
{
	/* floorplan input file */
	char flp_file[STR_SIZE];
	/* input power trace file */
	char p_infile[STR_SIZE];
  /* input material properties file */
  char materials_file[STR_SIZE];
  /* output file for the temperature trace */
	char t_outfile[STR_SIZE];
	/* input configuration parameters from file	*/
	char config[STR_SIZE];
	/* output configuration parameters to file	*/
	char dump_config[STR_SIZE];
	/* input microchannel configuration file */
	int use_microchannels;


	/*BU_3D: Option to turn on heterogenous R-C assignment*/
	char detailed_3D[STR_SIZE];

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

/*
 * HotFloorplan is a temperature-aware floorplanning tool
 * that can be easily modified to optimize for any arbitrary
 * metric. This program reads in a floorplan description
 * file that has the area, aspect ratio and connectivity
 * information of a set of functional blocks. It also reads
 * in a file with the average power dissipation values for
 * those blocks and outputs the best floorplan it can find
 * that optimizes the specified metric.
 */
#include <stdio.h>
#include <string.h>
#ifdef _MSC_VER
#define strcasecmp    _stricmp
#define strncasecmp   _strnicmp
#else
#include <strings.h>
#endif

#include "flp.h"
#include "temperature.h"
#include "materials.h"
#include "wire.h"
#include "util.h"
#include "hotfloorplan.h"

void usage(int argc, char **argv)
{
	fprintf(stdout, "Usage: %s -f <file> -p <file> -o <file> [-c <file>] [-d <file>] [options]\n", argv[0]);
	fprintf(stdout, "Finds a thermally-aware floorplan for the given functional blocks.\n");
	fprintf(stdout, "Options:(may be specified in any order, within \"[]\" means optional)\n");
	fprintf(stdout, "   -f <file>\tfloorplan description input file (e.g. ev6.desc)\n");
	fprintf(stdout, "   -p <file>\tpower input file (e.g. avg.p)\n");
	fprintf(stdout, "   -o <file>\tfloorplan output file\n");
	fprintf(stdout, "  [-c <file>]\tinput configuration parameters from file (e.g. hotspot.config)\n");
	fprintf(stdout, "  [-d <file>]\toutput configuration parameters to file\n");
	fprintf(stdout, "  [options]\tzero or more options of the form \"-<name> <value>\",\n");
	fprintf(stdout, "           \toverride the options from config file\n");
}

/*
 * parse a table of name-value string pairs and add the configuration
 * parameters to 'config'
 */
void global_config_from_strs(global_config_t *config, str_pair *table, int size)
{
	int idx;
	if ((idx = get_str_index(table, size, "f")) >= 0) {
		if(sscanf(table[idx].value, "%s", config->flp_desc) != 1)
			fatal("invalid format for configuration  parameter flp_desc\n");
	} else {
		fatal("required parameter flp_desc missing. check usage\n");
	}
	if ((idx = get_str_index(table, size, "p")) >= 0) {
		if(sscanf(table[idx].value, "%s", config->power_in) != 1)
			fatal("invalid format for configuration  parameter power_in\n");
	} else {
		fatal("required parameter power_in missing. check usage\n");
	}
	if ((idx = get_str_index(table, size, "o")) >= 0) {
		if(sscanf(table[idx].value, "%s", config->flp_out) != 1)
			fatal("invalid format for configuration  parameter flp_out\n");
	} else {
		fatal("required parameter flp_out missing. check usage\n");
	}
	if ((idx = get_str_index(table, size, "c")) >= 0) {
		if(sscanf(table[idx].value, "%s", config->config) != 1)
			fatal("invalid format for configuration  parameter config\n");
	} else {
		strcpy(config->config, NULLFILE);
	}
	if ((idx = get_str_index(table, size, "d")) >= 0) {
		if(sscanf(table[idx].value, "%s", config->dump_config) != 1)
			fatal("invalid format for configuration  parameter dump_config\n");
	} else {
		strcpy(config->dump_config, NULLFILE);
	}
  if ((idx = get_str_index(table, size, "materials_file")) >= 0) {
      if(sscanf(table[idx].value, "%s", config->materials_file) != 1)
        fatal("invalid format for configuration parameter materials_file\n");
  } else {
      strcpy(config->materials_file, NULLFILE);
  }
}

/*
 * convert config into a table of name-value pairs. returns the no.
 * of parameters converted
 */
int global_config_to_strs(global_config_t *config, str_pair *table, int max_entries)
{
	if (max_entries < 5)
		fatal("not enough entries in table\n");

	sprintf(table[0].name, "f");
	sprintf(table[1].name, "p");
	sprintf(table[2].name, "o");
	sprintf(table[3].name, "c");
	sprintf(table[4].name, "d");

	sprintf(table[0].value, "%s", config->flp_desc);
	sprintf(table[1].value, "%s", config->power_in);
	sprintf(table[2].value, "%s", config->flp_out);
	sprintf(table[3].value, "%s", config->config);
	sprintf(table[4].value, "%s", config->dump_config);

	return 5;
}

void print_wire_delays(flp_t *flp, double frequency)
{
	int i, j;
	double delay_g, delay_i;
	fprintf(stdout, "printing wire delay between blocks for global and intermediate metal layers:\n");
	fprintf(stdout, "(in %.1f GHz cycles)\n", frequency / 1.0e9);
	fprintf(stdout, "name1\tname2\tglobal\tintermediate\n");
	for (i=0; i < flp->n_units-1; i++)
		for (j=i+1; j < flp->n_units; j++)
			if (flp->wire_density[i][j]) {
				delay_g = wire_length2delay(get_manhattan_dist(flp, i, j), WIRE_GLOBAL);
				delay_i = wire_length2delay(get_manhattan_dist(flp, i, j), WIRE_INTER);
				fprintf(stdout, "%s\t%s\t%.3f\t%.3f\n", flp->units[i].name, flp->units[j].name,
						frequency * delay_g, frequency * delay_i);
			}
}

/* main function for the floorplanner	*/
int main(int argc, char **argv)
{
	flp_desc_t *flp_desc;
	flp_t *flp;
	RC_model_t *model;
	double *power;
	thermal_config_t thermal_config;
	flp_config_t flp_config;
	global_config_t global_config;
  materials_list_t materials_list;
	str_pair table[MAX_ENTRIES];
	int size, compacted;

	if (!(argc >= 7 && argc % 2)) {
		usage(argc, argv);
		return 1;
	}

	size = parse_cmdline(table, MAX_ENTRIES, argc, argv);
	global_config_from_strs(&global_config, table, size);

	/* read configuration file	*/
	if (strcmp(global_config.config, NULLFILE))
		size += read_str_pairs(&table[size], MAX_ENTRIES, global_config.config);

	/*
	 * in the str_pair 'table', earlier entries override later ones.
	 * so, command line options have priority over config file
	 */
	size = str_pairs_remove_duplicates(table, size);

  // fill in material properties
  default_materials(&materials_list);
  if(strncmp(global_config.materials_file, NULLFILE, STR_SIZE)) {
    materials_add_from_file(&materials_list, global_config.materials_file);
  }

	/* get defaults */
	thermal_config = default_thermal_config();
	flp_config = default_flp_config();
	/* modify according to command line / config file	*/
	thermal_config_add_from_strs(&thermal_config, &materials_list, table, size);
	flp_config_add_from_strs(&flp_config, table, size);

	/* dump configuration if specified	*/
	if (strcmp(global_config.dump_config, NULLFILE)) {
		size = global_config_to_strs(&global_config, table, MAX_ENTRIES);
		size += thermal_config_to_strs(&thermal_config, &table[size], MAX_ENTRIES-size);
		size += flp_config_to_strs(&flp_config, &table[size], MAX_ENTRIES-size);
		/* prefix the name of the variable with a '-'	*/
		dump_str_pairs(table, size, global_config.dump_config, "-");
	}

	/* If the grid model is used, things could be really slow!
	 * Also make sure it is not in the 3-d chip mode (specified
	 * with the layer configuration file)
	 */
	if (!strcmp(thermal_config.model_type, GRID_MODEL_STR)) {
		if (strcmp(thermal_config.grid_layer_file, NULLFILE))
			fatal("lcf file specified with the grid model. 3-d chips not supported\n");
		warning("grid model is used. HotFloorplan could be REALLY slow\n");
	}

	/* description of the functional blocks to be floorplanned	*/
	flp_desc = read_flp_desc(global_config.flp_desc, &flp_config);
	/*
	 * just an empty frame with blocks' names.
	 * block positions not known yet.
	 */
	flp = flp_placeholder(flp_desc);
	/* temperature model	*/
	model = alloc_RC_model(&thermal_config, flp, NULL, &materials_list, 0, 0);
	/* input power vector	*/
	power = hotspot_vector(model);
	read_power(model, power, global_config.power_in);

	/* main floorplanning routine	*/
	compacted = floorplan(flp, flp_desc, model, power);
	/*
	 * print the finally selected floorplan in a format that can
	 * be understood by tofig.pl (which converts it to a FIG file)
	 */
	print_flp_fig(flp);
	/* print the floorplan statistics	*/
	if (flp_config.wrap_l2 &&
		!strcasecmp(flp_desc->units[flp_desc->n_units-1].name, flp_config.l2_label))
		print_flp_stats(flp, model, flp_config.l2_label,
						global_config.power_in, global_config.flp_desc);

	/* print the wire delays between blocks	*/
	print_wire_delays(flp, thermal_config.base_proc_freq);

	/* also output the floorplan to a file readable by hotspot	*/
	dump_flp(flp, global_config.flp_out, FALSE);

	free_flp_desc(flp_desc);
	delete_RC_model(model);
	free_dvector(power);

	/* while deallocating, free the compacted blocks too	*/
	free_flp(flp, compacted, TRUE);

	return 0;
}

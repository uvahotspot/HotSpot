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

#include "temperature.h"
#include "temperature_block.h"
#include "temperature_grid.h"
#include "flp.h"
#include "util.h"

/* default thermal configuration parameters	*/
thermal_config_t default_thermal_config(void)
{
	thermal_config_t config;

	/* chip specs	*/
	config.t_chip = 0.15e-3;			/* chip thickness in meters	*/
	config.k_chip = 100.0; /* chip thermal conductivity in W/(m-K) */
	config.p_chip = 1.75e6; /* chip specific heat in J/(m^3-K) */
	/* temperature threshold for DTM (Kelvin)*/
	config.thermal_threshold = 81.8 + 273.15;

	/* heat sink specs	*/
	config.c_convec = 140.4;			/* convection capacitance in J/K */
	config.r_convec = 0.1;				/* convection resistance in K/W	*/
	config.s_sink = 60e-3;				/* heatsink side in m	*/
	config.t_sink = 6.9e-3; 			/* heatsink thickness  in m	*/
	config.k_sink = 400.0; /* heatsink thermal conductivity in W/(m-K) */
	config.p_sink = 3.55e6; /* heatsink specific heat in J/(m^3-K) */


	/* heat spreader specs	*/
	config.s_spreader = 30e-3;			/* spreader side in m	*/
	config.t_spreader = 1e-3;			/* spreader thickness in m	*/
	config.k_spreader = 400.0; /* heat spreader thermal conductivity in W/(m-K) */
	config.p_spreader = 3.55e6; /* heat spreader specific heat in J/(m^3-K) */

	/* interface material specs	*/
	config.t_interface = 20e-6;			/* interface material thickness in m */
	config.k_interface = 4.0; /* interface material thermal conductivity in W/(m-K) */
	config.p_interface = 4.0e6; /* interface material specific heat in J/(m^3-K) */

	/* secondary heat transfer path */
	config.model_secondary = FALSE;
	config.r_convec_sec = 1.0;
	config.c_convec_sec = 140.4; //FIXME! need updated value.
	config.n_metal = 8;
	config.t_metal = 10.0e-6;
	config.t_c4 = 0.0001;
	config.s_c4 = 20.0e-6;
	config.n_c4 = 400;
	config.s_sub = 0.021;
	config.t_sub = 0.001;
	config.s_solder = 0.021;
	config.t_solder = 0.00094;
	config.s_pcb = 0.1;
	config.t_pcb = 0.002;

	/* others	*/
	config.ambient = 45 + 273.15;		/* in kelvin	*/
	/* initial temperatures	from file	*/
	strcpy(config.init_file, NULLFILE);
	config.init_temp = 60 + 273.15;		/* in Kelvin	*/
	/* steady state temperatures to file	*/
	strcpy(config.steady_file, NULLFILE);
  /* transient grid temperatures to file */
  strcpy(config.grid_transient_file, NULLFILE);
 	/* 3.33 us sampling interval = 10K cycles at 3GHz	*/
	config.sampling_intvl = 3.333e-6;
	config.base_proc_freq = 3e9;		/* base processor frequency in Hz	*/
	config.dtm_used = FALSE;			/* set accordingly	*/

	config.leakage_used = 0;
	config.leakage_mode = 0;

	config.package_model_used = 0;
	strcpy(config.package_config_file, NULLFILE);

	/* set block model as default	*/
	strcpy(config.model_type, BLOCK_MODEL_STR);

	/* block model specific parameters	*/
	config.block_omit_lateral = FALSE;	/* omit lateral chip resistances?	*/

	/* grid model specific parameters	*/
	config.grid_rows = 64;				/* grid resolution - no. of rows	*/
	config.grid_cols = 64;				/* grid resolution - no. of cols	*/
	/* layer configuration from	file */
	strcpy(config.grid_layer_file, NULLFILE);
	/* output steady state grid temperatures apart from block temperatures */
	strcpy(config.grid_steady_file, NULLFILE);
	/*
	 * mapping mode between block and grid models.
	 * default: use the temperature of the center
	 * grid cell as that of the entire block
	 */
	strcpy(config.grid_map_mode, GRID_CENTER_STR);

	config.detailed_3D_used = 0;	//BU_3D: by default detailed 3D modeling is disabled.

	return config;
}

/*
 * parse a table of name-value string pairs and add the configuration
 * parameters to 'config'
 */
void thermal_config_add_from_strs(thermal_config_t *config, materials_list_t *materials_list, str_pair *table, int size)
{
	int idx;
	if ((idx = get_str_index(table, size, "t_chip")) >= 0)
		if(sscanf(table[idx].value, "%lf", &config->t_chip) != 1)
			fatal("invalid format for configuration  parameter t_chip\n");
	if ((idx = get_str_index(table, size, "k_chip")) >= 0)
		if(sscanf(table[idx].value, "%lf", &config->k_chip) != 1)
			fatal("invalid format for configuration  parameter k_chip\n");
	if ((idx = get_str_index(table, size, "p_chip")) >= 0)
		if(sscanf(table[idx].value, "%lf", &config->p_chip) != 1)
			fatal("invalid format for configuration  parameter p_chip\n");
	if ((idx = get_str_index(table, size, "thermal_threshold")) >= 0)
		if(sscanf(table[idx].value, "%lf", &config->thermal_threshold) != 1)
			fatal("invalid format for configuration  parameter thermal_threshold\n");
	if ((idx = get_str_index(table, size, "c_convec")) >= 0)
		if(sscanf(table[idx].value, "%lf", &config->c_convec) != 1)
			fatal("invalid format for configuration  parameter c_convec\n");
	if ((idx = get_str_index(table, size, "r_convec")) >= 0)
		if(sscanf(table[idx].value, "%lf", &config->r_convec) != 1)
			fatal("invalid format for configuration  parameter r_convec\n");
	if ((idx = get_str_index(table, size, "s_sink")) >= 0)
		if(sscanf(table[idx].value, "%lf", &config->s_sink) != 1)
			fatal("invalid format for configuration  parameter s_sink\n");
	if ((idx = get_str_index(table, size, "t_sink")) >= 0)
		if(sscanf(table[idx].value, "%lf", &config->t_sink) != 1)
			fatal("invalid format for configuration  parameter t_sink\n");
	if ((idx = get_str_index(table, size, "k_sink")) >= 0)
		if(sscanf(table[idx].value, "%lf", &config->k_sink) != 1)
			fatal("invalid format for configuration  parameter k_sink\n");
	if ((idx = get_str_index(table, size, "p_sink")) >= 0)
		if(sscanf(table[idx].value, "%lf", &config->p_sink) != 1)
			fatal("invalid format for configuration  parameter p_sink\n");
	if ((idx = get_str_index(table, size, "s_spreader")) >= 0)
		if(sscanf(table[idx].value, "%lf", &config->s_spreader) != 1)
			fatal("invalid format for configuration  parameter s_spreader\n");
	if ((idx = get_str_index(table, size, "t_spreader")) >= 0)
		if(sscanf(table[idx].value, "%lf", &config->t_spreader) != 1)
			fatal("invalid format for configuration  parameter t_spreader\n");
	if ((idx = get_str_index(table, size, "k_spreader")) >= 0)
		if(sscanf(table[idx].value, "%lf", &config->k_spreader) != 1)
			fatal("invalid format for configuration  parameter k_spreader\n");
	if ((idx = get_str_index(table, size, "p_spreader")) >= 0)
		if(sscanf(table[idx].value, "%lf", &config->p_spreader) != 1)
			fatal("invalid format for configuration  parameter p_spreader\n");
	if ((idx = get_str_index(table, size, "t_interface")) >= 0)
		if(sscanf(table[idx].value, "%lf", &config->t_interface) != 1)
			fatal("invalid format for configuration  parameter t_interface\n");
	if ((idx = get_str_index(table, size, "k_interface")) >= 0)
		if(sscanf(table[idx].value, "%lf", &config->k_interface) != 1)
			fatal("invalid format for configuration  parameter k_interface\n");
	if ((idx = get_str_index(table, size, "p_interface")) >= 0)
		if(sscanf(table[idx].value, "%lf", &config->p_interface) != 1)
			fatal("invalid format for configuration  parameter p_interface\n");
	if ((idx = get_str_index(table, size, "model_secondary")) >= 0)
		if(sscanf(table[idx].value, "%d", &config->model_secondary) != 1)
			fatal("invalid format for configuration  parameter model_secondary\n");
	if ((idx = get_str_index(table, size, "r_convec_sec")) >= 0)
		if(sscanf(table[idx].value, "%lf", &config->r_convec_sec) != 1)
			fatal("invalid format for configuration  parameter r_convec_sec\n");
	if ((idx = get_str_index(table, size, "c_convec_sec")) >= 0)
		if(sscanf(table[idx].value, "%lf", &config->c_convec_sec) != 1)
			fatal("invalid format for configuration  parameter c_convec_sec\n");
	if ((idx = get_str_index(table, size, "n_metal")) >= 0)
		if(sscanf(table[idx].value, "%d", &config->n_metal) != 1)
			fatal("invalid format for configuration  parameter n_metal\n");
	if ((idx = get_str_index(table, size, "t_metal")) >= 0)
		if(sscanf(table[idx].value, "%lf", &config->t_metal) != 1)
			fatal("invalid format for configuration  parameter t_metal\n");
	if ((idx = get_str_index(table, size, "t_c4")) >= 0)
		if(sscanf(table[idx].value, "%lf", &config->t_c4) != 1)
			fatal("invalid format for configuration  parameter t_c4\n");
	if ((idx = get_str_index(table, size, "s_c4")) >= 0)
		if(sscanf(table[idx].value, "%lf", &config->s_c4) != 1)
			fatal("invalid format for configuration  parameter s_c4\n");
	if ((idx = get_str_index(table, size, "n_c4")) >= 0)
		if(sscanf(table[idx].value, "%d", &config->n_c4) != 1)
			fatal("invalid format for configuration  parameter n_c4\n");
	if ((idx = get_str_index(table, size, "s_sub")) >= 0)
		if(sscanf(table[idx].value, "%lf", &config->s_sub) != 1)
			fatal("invalid format for configuration  parameter s_sub\n");
	if ((idx = get_str_index(table, size, "t_sub")) >= 0)
		if(sscanf(table[idx].value, "%lf", &config->t_sub) != 1)
			fatal("invalid format for configuration  parameter t_sub\n");
	if ((idx = get_str_index(table, size, "s_solder")) >= 0)
		if(sscanf(table[idx].value, "%lf", &config->s_solder) != 1)
			fatal("invalid format for configuration  parameter s_solder\n");
	if ((idx = get_str_index(table, size, "t_solder")) >= 0)
		if(sscanf(table[idx].value, "%lf", &config->t_solder) != 1)
			fatal("invalid format for configuration  parameter t_solder\n");
	if ((idx = get_str_index(table, size, "s_pcb")) >= 0)
		if(sscanf(table[idx].value, "%lf", &config->s_pcb) != 1)
			fatal("invalid format for configuration  parameter s_pcb\n");
	if ((idx = get_str_index(table, size, "t_pcb")) >= 0)
		if(sscanf(table[idx].value, "%lf", &config->t_pcb) != 1)
			fatal("invalid format for configuration  parameter t_pcb\n");
	if ((idx = get_str_index(table, size, "ambient")) >= 0)
		if(sscanf(table[idx].value, "%lf", &config->ambient) != 1)
			fatal("invalid format for configuration  parameter ambient\n");
	if ((idx = get_str_index(table, size, "init_file")) >= 0)
		if(sscanf(table[idx].value, "%s", config->init_file) != 1)
			fatal("invalid format for configuration  parameter init_file\n");
	if ((idx = get_str_index(table, size, "init_temp")) >= 0)
		if(sscanf(table[idx].value, "%lf", &config->init_temp) != 1)
			fatal("invalid format for configuration  parameter init_temp\n");
	if ((idx = get_str_index(table, size, "steady_file")) >= 0)
		if(sscanf(table[idx].value, "%s", config->steady_file) != 1)
			fatal("invalid format for configuration  parameter steady_file\n");
  if ((idx = get_str_index(table, size, "grid_transient_file")) >= 0)
		if(sscanf(table[idx].value, "%s", config->grid_transient_file) != 1)
			fatal("invalid format for configuration  parameter grid_transient_file\n");
	if ((idx = get_str_index(table, size, "sampling_intvl")) >= 0)
		if(sscanf(table[idx].value, "%lf", &config->sampling_intvl) != 1)
			fatal("invalid format for configuration  parameter sampling_intvl\n");
	if ((idx = get_str_index(table, size, "base_proc_freq")) >= 0)
		if(sscanf(table[idx].value, "%lf", &config->base_proc_freq) != 1)
			fatal("invalid format for configuration  parameter base_proc_freq\n");
	if ((idx = get_str_index(table, size, "dtm_used")) >= 0)
		if(sscanf(table[idx].value, "%d", &config->dtm_used) != 1)
			fatal("invalid format for configuration  parameter dtm_used\n");
	if ((idx = get_str_index(table, size, "model_type")) >= 0)
		if(sscanf(table[idx].value, "%s", config->model_type) != 1)
			fatal("invalid format for configuration  parameter model_type\n");
	if ((idx = get_str_index(table, size, "leakage_used")) >= 0)
		if(sscanf(table[idx].value, "%d", &config->leakage_used) != 1)
			fatal("invalid format for configuration  parameter leakage_used\n");
	if ((idx = get_str_index(table, size, "leakage_mode")) >= 0)
		if(sscanf(table[idx].value, "%d", &config->leakage_mode) != 1)
			fatal("invalid format for configuration  parameter leakage_mode\n");
	if ((idx = get_str_index(table, size, "package_model_used")) >= 0)
		if(sscanf(table[idx].value, "%d", &config->package_model_used) != 1)
			fatal("invalid format for configuration  parameter package_model_used\n");
	if ((idx = get_str_index(table, size, "package_config_file")) >= 0)
		if(sscanf(table[idx].value, "%s", config->package_config_file) != 1)
			fatal("invalid format for configuration  parameter package_config_file\n");
	if ((idx = get_str_index(table, size, "block_omit_lateral")) >= 0)
		if(sscanf(table[idx].value, "%d", &config->block_omit_lateral) != 1)
			fatal("invalid format for configuration  parameter block_omit_lateral\n");
	if ((idx = get_str_index(table, size, "grid_rows")) >= 0)
		if(sscanf(table[idx].value, "%d", &config->grid_rows) != 1)
			fatal("invalid format for configuration  parameter grid_rows\n");
	if ((idx = get_str_index(table, size, "grid_cols")) >= 0)
		if(sscanf(table[idx].value, "%d", &config->grid_cols) != 1)
			fatal("invalid format for configuration  parameter grid_cols\n");
	if ((idx = get_str_index(table, size, "grid_layer_file")) >= 0)
		if(sscanf(table[idx].value, "%s", config->grid_layer_file) != 1)
			fatal("invalid format for configuration  parameter grid_layer_file\n");
	if ((idx = get_str_index(table, size, "grid_steady_file")) >= 0)
		if(sscanf(table[idx].value, "%s", config->grid_steady_file) != 1)
			fatal("invalid format for configuration  parameter grid_steady_file\n");
	if ((idx = get_str_index(table, size, "grid_map_mode")) >= 0)
		if(sscanf(table[idx].value, "%s", config->grid_map_mode) != 1)
			fatal("invalid format for configuration  parameter grid_map_mode\n");

	if ((config->t_chip <= 0) || (config->s_sink <= 0) || (config->t_sink <= 0) ||
		(config->s_spreader <= 0) || (config->t_spreader <= 0) ||
		(config->t_interface <= 0))
		fatal("chip and package dimensions should be greater than zero\n");
	if ((config->t_metal <= 0) || (config->n_metal <= 0) || (config->t_c4 <= 0) ||
		(config->s_c4 <= 0) || (config->n_c4 <= 0) || (config->s_sub <= 0) || (config->t_sub <= 0) ||
		(config->s_solder <= 0) || (config->t_solder <= 0) || (config->s_pcb <= 0) ||
		(config->t_solder <= 0) || (config->r_convec_sec <= 0) || (config->c_convec_sec <= 0))
		fatal("secondary heat tranfer layer dimensions should be greater than zero\n");
	/* leakage iteration is not supported in transient mode in this release */
	if (config->leakage_used == 1) {
		printf("Warning: transient leakage iteration is not supported in this release...\n");
		printf(" ...all transient results are without thermal-leakage loop.\n");
	}
	if ((config->model_secondary == 1) && (!strcasecmp(config->model_type, BLOCK_MODEL_STR)))
		fatal("secondary heat tranfer path is supported only in the grid mode\n");
	if ((config->thermal_threshold < 0) || (config->c_convec < 0) ||
		(config->r_convec < 0) || (config->ambient < 0) ||
		(config->base_proc_freq <= 0) || (config->sampling_intvl <= 0))
		fatal("invalid thermal simulation parameters\n");
	if (strcasecmp(config->model_type, BLOCK_MODEL_STR) &&
		strcasecmp(config->model_type, GRID_MODEL_STR))
		fatal("invalid model type. use 'block' or 'grid'\n");
	if(config->grid_rows <= 0 || config->grid_cols <= 0)
		fatal("grid rows and columns should both be greater than zero\n");
	if (strcasecmp(config->grid_map_mode, GRID_AVG_STR) &&
		strcasecmp(config->grid_map_mode, GRID_MIN_STR) &&
		strcasecmp(config->grid_map_mode, GRID_MAX_STR) &&
		strcasecmp(config->grid_map_mode, GRID_CENTER_STR))
		fatal("invalid mapping mode. use 'avg', 'min', 'max' or 'center'\n");

  if ((idx = get_str_index(table, size, "material_chip")) >= 0) {
    char material_name[STR_SIZE];
    if(sscanf(table[idx].value, "%s", material_name) != 1)
      fatal("invalid format for configuration parameter material_chip\n");
    config->k_chip = get_material_thermal_conductivity(materials_list, material_name);
    config->p_chip = get_material_volumetric_heat_capacity(materials_list, material_name);
    if(config->k_chip < 0 || config->p_chip < 0)
      fatal("material name specified in configuration parameter material_chip not found\n");
  }
  if ((idx = get_str_index(table, size, "material_sink")) >= 0) {
    char material_name[STR_SIZE];
    if(sscanf(table[idx].value, "%s", material_name) != 1)
      fatal("invalid format for configuration parameter material_sink\n");
    config->k_sink = get_material_thermal_conductivity(materials_list, material_name);
    config->p_sink = get_material_volumetric_heat_capacity(materials_list, material_name);

    if(config->k_sink < 0 || config->p_sink < 0)
      fatal("material name specified in configuration parameter material_sink not found\n");
  }
  if ((idx = get_str_index(table, size, "material_spreader")) >= 0) {
    char material_name[STR_SIZE];
    if(sscanf(table[idx].value, "%s", material_name) != 1)
      fatal("invalid format for configuration parameter material_spreader\n");
    config->k_spreader = get_material_thermal_conductivity(materials_list, material_name);
    config->p_spreader = get_material_volumetric_heat_capacity(materials_list, material_name);

    if(config->k_spreader < 0 || config->p_spreader < 0)
      fatal("material name specified in configuration parameter material_spreader not found\n");
  }
  if ((idx = get_str_index(table, size, "material_interface")) >= 0) {
    char material_name[STR_SIZE];
    if(sscanf(table[idx].value, "%s", material_name) != 1)
      fatal("invalid format for configuration parameter material_interface\n");
    config->k_interface = get_material_thermal_conductivity(materials_list, material_name);
    config->p_interface = get_material_volumetric_heat_capacity(materials_list, material_name);

    if(config->k_interface < 0 || config->p_interface < 0)
      fatal("material name specified in configuration parameter material_interface not found\n");
  }
}

/*
 * convert config into a table of name-value pairs. returns the no.
 * of parameters converted
 */
int thermal_config_to_strs(thermal_config_t *config, str_pair *table, int max_entries)
{
	if (max_entries < 51)
		fatal("not enough entries in table\n");

	sprintf(table[0].name, "t_chip");
	sprintf(table[1].name, "k_chip");
	sprintf(table[2].name, "p_chip");
	sprintf(table[3].name, "thermal_threshold");
	sprintf(table[4].name, "c_convec");
	sprintf(table[5].name, "r_convec");
	sprintf(table[6].name, "s_sink");
	sprintf(table[7].name, "t_sink");
	sprintf(table[8].name, "k_sink");
	sprintf(table[9].name, "p_sink");
	sprintf(table[10].name, "s_spreader");
	sprintf(table[11].name, "t_spreader");
	sprintf(table[12].name, "k_spreader");
	sprintf(table[13].name, "p_spreader");
	sprintf(table[14].name, "t_interface");
	sprintf(table[15].name, "k_interface");
	sprintf(table[16].name, "p_interface");
	sprintf(table[17].name, "model_secondary");
	sprintf(table[18].name, "r_convec_sec");
	sprintf(table[19].name, "c_convec_sec");
	sprintf(table[20].name, "n_metal");
	sprintf(table[21].name, "t_metal");
	sprintf(table[22].name, "t_c4");
	sprintf(table[23].name, "s_c4");
	sprintf(table[24].name, "n_c4");
	sprintf(table[25].name, "s_sub");
	sprintf(table[26].name, "t_sub");
	sprintf(table[27].name, "s_solder");
	sprintf(table[28].name, "t_solder");
	sprintf(table[29].name, "s_pcb");
	sprintf(table[30].name, "t_pcb");
	sprintf(table[31].name, "ambient");
	sprintf(table[32].name, "init_file");
	sprintf(table[33].name, "init_temp");
	sprintf(table[34].name, "steady_file");
	sprintf(table[35].name, "sampling_intvl");
	sprintf(table[36].name, "base_proc_freq");
	sprintf(table[37].name, "dtm_used");
	sprintf(table[38].name, "model_type");
	sprintf(table[39].name, "leakage_used");
	sprintf(table[40].name, "leakage_mode");
	sprintf(table[41].name, "package_model_used");
	sprintf(table[42].name, "package_config_file");
	sprintf(table[43].name, "block_omit_lateral");
	sprintf(table[44].name, "grid_rows");
	sprintf(table[45].name, "grid_cols");
	sprintf(table[46].name, "grid_layer_file");
	sprintf(table[47].name, "grid_steady_file");
	sprintf(table[48].name, "grid_map_mode");
  sprintf(table[49].name, "grid_transient_file");
  sprintf(table[50].name, "detailed_3D_used");

	sprintf(table[0].value, "%lg", config->t_chip);
	sprintf(table[1].value, "%lg", config->k_chip);
	sprintf(table[2].value, "%lg", config->p_chip);
	sprintf(table[3].value, "%lg", config->thermal_threshold);
	sprintf(table[4].value, "%lg", config->c_convec);
	sprintf(table[5].value, "%lg", config->r_convec);
	sprintf(table[6].value, "%lg", config->s_sink);
	sprintf(table[7].value, "%lg", config->t_sink);
	sprintf(table[8].value, "%lg", config->k_sink);
	sprintf(table[9].value, "%lg", config->p_sink);
	sprintf(table[10].value, "%lg", config->s_spreader);
	sprintf(table[11].value, "%lg", config->t_spreader);
	sprintf(table[12].value, "%lg", config->k_spreader);
	sprintf(table[13].value, "%lg", config->p_spreader);
	sprintf(table[14].value, "%lg", config->t_interface);
	sprintf(table[15].value, "%lg", config->k_interface);
	sprintf(table[16].value, "%lg", config->p_interface);
	sprintf(table[17].value, "%d", config->model_secondary);
	sprintf(table[18].value, "%lg", config->r_convec_sec);
	sprintf(table[19].value, "%lg", config->c_convec_sec);
	sprintf(table[20].value, "%d", config->n_metal);
	sprintf(table[21].value, "%lg", config->t_metal);
	sprintf(table[22].value, "%lg", config->t_c4);
	sprintf(table[23].value, "%lg", config->s_c4);
	sprintf(table[24].value, "%d", config->n_c4);
	sprintf(table[25].value, "%lg", config->s_sub);
	sprintf(table[26].value, "%lg", config->t_sub);
	sprintf(table[27].value, "%lg", config->s_solder);
	sprintf(table[28].value, "%lg", config->t_solder);
	sprintf(table[29].value, "%lg", config->s_pcb);
	sprintf(table[30].value, "%lg", config->t_pcb);
	sprintf(table[31].value, "%lg", config->ambient);
	sprintf(table[32].value, "%s", config->init_file);
	sprintf(table[33].value, "%lg", config->init_temp);
	sprintf(table[34].value, "%s", config->steady_file);
	sprintf(table[35].value, "%lg", config->sampling_intvl);
	sprintf(table[36].value, "%lg", config->base_proc_freq);
	sprintf(table[37].value, "%d", config->dtm_used);
	sprintf(table[38].value, "%s", config->model_type);
	sprintf(table[39].value, "%d", config->leakage_used);
	sprintf(table[40].value, "%d", config->leakage_mode);
	sprintf(table[41].value, "%d", config->package_model_used);
	sprintf(table[42].value, "%s", config->package_config_file);
	sprintf(table[43].value, "%d", config->block_omit_lateral);
	sprintf(table[44].value, "%d", config->grid_rows);
	sprintf(table[45].value, "%d", config->grid_cols);
	sprintf(table[46].value, "%s", config->grid_layer_file);
	sprintf(table[47].value, "%s", config->grid_steady_file);
	sprintf(table[48].value, "%s", config->grid_map_mode);
  sprintf(table[49].value, "%s", config->grid_transient_file);
  sprintf(table[50].value, "%d", config->detailed_3D_used);

	return 51;
}

/* package parameter routines	*/
void populate_package_R(package_RC_t *p, thermal_config_t *config, double width, double height)
{
	double s_spreader = config->s_spreader;
	double t_spreader = config->t_spreader;
	double s_sink = config->s_sink;
	double t_sink = config->t_sink;
	double r_convec = config->r_convec;

	double s_sub = config->s_sub;
	double t_sub = config->t_sub;
	double s_solder = config->s_solder;
	double t_solder = config->t_solder;
	double s_pcb = config->s_pcb;
	double t_pcb = config->t_pcb;
	double r_convec_sec = config->r_convec_sec;

	double k_sink = config->k_sink;
	double k_spreader = config->k_spreader;


	/* lateral R's of spreader and sink */
	p->r_sp1_x = getr(k_spreader, (s_spreader-width)/4.0, (s_spreader+3*height)/4.0 * t_spreader);
	p->r_sp1_y = getr(k_spreader, (s_spreader-height)/4.0, (s_spreader+3*width)/4.0 * t_spreader);
	p->r_hs1_x = getr(k_sink, (s_spreader-width)/4.0, (s_spreader+3*height)/4.0 * t_sink);
	p->r_hs1_y = getr(k_sink, (s_spreader-height)/4.0, (s_spreader+3*width)/4.0 * t_sink);
	p->r_hs2_x = getr(k_sink, (s_spreader-width)/4.0, (3*s_spreader+height)/4.0 * t_sink);
	p->r_hs2_y = getr(k_sink, (s_spreader-height)/4.0, (3*s_spreader+width)/4.0 * t_sink);
	p->r_hs = getr(k_sink, (s_sink-s_spreader)/4.0, (s_sink+3*s_spreader)/4.0 * t_sink);

	/* vertical R's of spreader and sink */
	p->r_sp_per_x = getr(k_spreader, t_spreader, (s_spreader+height) * (s_spreader-width) / 4.0);
	p->r_sp_per_y = getr(k_spreader, t_spreader, (s_spreader+width) * (s_spreader-height) / 4.0);
	p->r_hs_c_per_x = getr(k_sink, t_sink, (s_spreader+height) * (s_spreader-width) / 4.0);
	p->r_hs_c_per_y = getr(k_sink, t_sink, (s_spreader+width) * (s_spreader-height) / 4.0);
	p->r_hs_per = getr(k_sink, t_sink, (s_sink*s_sink - s_spreader*s_spreader) / 4.0);

	/* vertical R's to ambient (divide r_convec proportional to area) */
	p->r_amb_c_per_x = r_convec * (s_sink * s_sink) / ((s_spreader+height) * (s_spreader-width) / 4.0);
	p->r_amb_c_per_y = r_convec * (s_sink * s_sink) / ((s_spreader+width) * (s_spreader-height) / 4.0);
	p->r_amb_per = r_convec * (s_sink * s_sink) / ((s_sink*s_sink - s_spreader*s_spreader) / 4.0);

	/* lateral R's of package substrate, solder and PCB */
	p->r_sub1_x = getr(K_SUB, (s_sub-width)/4.0, (s_sub+3*height)/4.0 * t_sub);
	p->r_sub1_y = getr(K_SUB, (s_sub-height)/4.0, (s_sub+3*width)/4.0 * t_sub);
	p->r_solder1_x = getr(K_SOLDER, (s_solder-width)/4.0, (s_solder+3*height)/4.0 * t_solder);
	p->r_solder1_y = getr(K_SOLDER, (s_solder-height)/4.0, (s_solder+3*width)/4.0 * t_solder);
	p->r_pcb1_x = getr(K_PCB, (s_solder-width)/4.0, (s_solder+3*height)/4.0 * t_pcb);
	p->r_pcb1_y = getr(K_PCB, (s_solder-height)/4.0, (s_solder+3*width)/4.0 * t_pcb);
	p->r_pcb2_x = getr(K_PCB, (s_solder-width)/4.0, (3*s_solder+height)/4.0 * t_pcb);
	p->r_pcb2_y = getr(K_PCB, (s_solder-height)/4.0, (3*s_solder+width)/4.0 * t_pcb);
	p->r_pcb = getr(K_PCB, (s_pcb-s_solder)/4.0, (s_pcb+3*s_solder)/4.0 * t_pcb);

	/* vertical R's of package substrate, solder balls and PCB */
	p->r_sub_per_x = getr(K_SUB, t_sub, (s_sub+height) * (s_sub-width) / 4.0);
	p->r_sub_per_y = getr(K_SUB, t_sub, (s_sub+width) * (s_sub-height) / 4.0);
	p->r_solder_per_x = getr(K_SOLDER, t_solder, (s_solder+height) * (s_solder-width) / 4.0);
	p->r_solder_per_y = getr(K_SOLDER, t_solder, (s_solder+width) * (s_solder-height) / 4.0);
	p->r_pcb_c_per_x = getr(K_PCB, t_pcb, (s_solder+height) * (s_solder-width) / 4.0);
	p->r_pcb_c_per_y = getr(K_PCB, t_pcb, (s_solder+width) * (s_solder-height) / 4.0);
	p->r_pcb_per = getr(K_PCB, t_pcb, (s_pcb*s_pcb - s_solder*s_solder) / 4.0);

	/* vertical R's to ambient at PCB (divide r_convec_sec proportional to area) */
	p->r_amb_sec_c_per_x = r_convec_sec * (s_pcb * s_pcb) / ((s_solder+height) * (s_solder-width) / 4.0);
	p->r_amb_sec_c_per_y = r_convec_sec * (s_pcb * s_pcb) / ((s_solder+width) * (s_solder-height) / 4.0);
	p->r_amb_sec_per = r_convec_sec * (s_pcb * s_pcb) / ((s_pcb*s_pcb - s_solder*s_solder) / 4.0);
}

void populate_package_C(package_RC_t *p, thermal_config_t *config, double width, double height)
{
	double s_spreader = config->s_spreader;
	double t_spreader = config->t_spreader;
	double s_sink = config->s_sink;
	double t_sink = config->t_sink;
	double c_convec = config->c_convec;

	double s_sub = config->s_sub;
	double t_sub = config->t_sub;
	double s_solder = config->s_solder;
	double t_solder = config->t_solder;
	double s_pcb = config->s_pcb;
	double t_pcb = config->t_pcb;
	double c_convec_sec = config->c_convec_sec;

	double p_sink = config->p_sink;
	double p_spreader = config->p_spreader;

	/* vertical C's of spreader and sink */
	p->c_sp_per_x = getcap(p_spreader, t_spreader, (s_spreader+height) * (s_spreader-width) / 4.0);
	p->c_sp_per_y = getcap(p_spreader, t_spreader, (s_spreader+width) * (s_spreader-height) / 4.0);
	p->c_hs_c_per_x = getcap(p_sink, t_sink, (s_spreader+height) * (s_spreader-width) / 4.0);
	p->c_hs_c_per_y = getcap(p_sink, t_sink, (s_spreader+width) * (s_spreader-height) / 4.0);
	p->c_hs_per = getcap(p_sink, t_sink, (s_sink*s_sink - s_spreader*s_spreader) / 4.0);

	/* vertical C's to ambient (divide c_convec proportional to area) */
	p->c_amb_c_per_x = C_FACTOR * c_convec / (s_sink * s_sink) * ((s_spreader+height) * (s_spreader-width) / 4.0);
	p->c_amb_c_per_y = C_FACTOR * c_convec / (s_sink * s_sink) * ((s_spreader+width) * (s_spreader-height) / 4.0);
	p->c_amb_per = C_FACTOR * c_convec / (s_sink * s_sink) * ((s_sink*s_sink - s_spreader*s_spreader) / 4.0);

	/* vertical C's of package substrate, solder balls, and PCB */
	p->c_sub_per_x = getcap(SPEC_HEAT_SUB, t_sub, (s_sub+height) * (s_sub-width) / 4.0);
	p->c_sub_per_y = getcap(SPEC_HEAT_SUB, t_sub, (s_sub+width) * (s_sub-height) / 4.0);
	p->c_solder_per_x = getcap(SPEC_HEAT_SOLDER, t_solder, (s_solder+height) * (s_solder-width) / 4.0);
	p->c_solder_per_y = getcap(SPEC_HEAT_SOLDER, t_solder, (s_solder+width) * (s_solder-height) / 4.0);
	p->c_pcb_c_per_x = getcap(SPEC_HEAT_PCB, t_pcb, (s_solder+height) * (s_solder-width) / 4.0);
	p->c_pcb_c_per_y = getcap(SPEC_HEAT_PCB, t_pcb, (s_solder+width) * (s_solder-height) / 4.0);
	p->c_pcb_per = getcap(SPEC_HEAT_PCB, t_pcb, (s_pcb*s_pcb - s_solder*s_solder) / 4.0);

	/* vertical C's to ambient at PCB (divide c_convec_sec proportional to area) */
	p->c_amb_sec_c_per_x = C_FACTOR * c_convec_sec / (s_pcb * s_pcb) * ((s_solder+height) * (s_solder-width) / 4.0);
	p->c_amb_sec_c_per_y = C_FACTOR * c_convec_sec / (s_pcb * s_pcb) * ((s_solder+width) * (s_solder-height) / 4.0);
	p->c_amb_sec_per = C_FACTOR * c_convec_sec / (s_pcb * s_pcb) * ((s_pcb*s_pcb - s_solder*s_solder) / 4.0);
}

/* debug print	*/
void debug_print_package_RC(package_RC_t *p)
{
	fprintf(stdout, "printing package RC information...\n");
	fprintf(stdout, "r_sp1_x: %f\tr_sp1_y: %f\n", p->r_sp1_x, p->r_sp1_y);
	fprintf(stdout, "r_sp_per_x: %f\tr_sp_per_y: %f\n", p->r_sp_per_x, p->r_sp_per_y);
	fprintf(stdout, "c_sp_per_x: %f\tc_sp_per_y: %f\n", p->c_sp_per_x, p->c_sp_per_y);
	fprintf(stdout, "r_hs1_x: %f\tr_hs1_y: %f\n", p->r_hs1_x, p->r_hs1_y);
	fprintf(stdout, "r_hs2_x: %f\tr_hs2_y: %f\n", p->r_hs2_x, p->r_hs2_y);
	fprintf(stdout, "r_hs_c_per_x: %f\tr_hs_c_per_y: %f\n", p->r_hs_c_per_x, p->r_hs_c_per_y);
	fprintf(stdout, "c_hs_c_per_x: %f\tc_hs_c_per_y: %f\n", p->c_hs_c_per_x, p->c_hs_c_per_y);
	fprintf(stdout, "r_hs: %f\tr_hs_per: %f\n", p->r_hs, p->r_hs_per);
	fprintf(stdout, "c_hs_per: %f\n", p->c_hs_per);
	fprintf(stdout, "r_amb_c_per_x: %f\tr_amb_c_per_y: %f\n", p->r_amb_c_per_x, p->r_amb_c_per_y);
	fprintf(stdout, "c_amb_c_per_x: %f\tc_amb_c_per_y: %f\n", p->c_amb_c_per_x, p->c_amb_c_per_y);
	fprintf(stdout, "r_amb_per: %f\n", p->r_amb_per);
	fprintf(stdout, "c_amb_per: %f\n", p->c_amb_per);
	fprintf(stdout, "r_sub1_x: %f\tr_sub1_y: %f\n", p->r_sub1_x, p->r_sub1_y);
	fprintf(stdout, "r_sub_per_x: %f\tr_sub_per_y: %f\n", p->r_sub_per_x, p->r_sub_per_y);
	fprintf(stdout, "c_sub_per_x: %f\tc_sub_per_y: %f\n", p->c_sub_per_x, p->c_sub_per_y);
	fprintf(stdout, "r_solder1_x: %f\tr_solder1_y: %f\n", p->r_solder1_x, p->r_solder1_y);
	fprintf(stdout, "r_solder_per_x: %f\tr_solder_per_y: %f\n", p->r_solder_per_x, p->r_solder_per_y);
	fprintf(stdout, "c_solder_per_x: %f\tc_solder_per_y: %f\n", p->c_solder_per_x, p->c_solder_per_y);
	fprintf(stdout, "r_pcb1_x: %f\tr_pcb1_y: %f\n", p->r_pcb1_x, p->r_pcb1_y);
	fprintf(stdout, "r_pcb2_x: %f\tr_pcb2_y: %f\n", p->r_pcb2_x, p->r_pcb2_y);
	fprintf(stdout, "r_pcb_c_per_x: %f\tr_pcb_c_per_y: %f\n", p->r_pcb_c_per_x, p->r_pcb_c_per_y);
	fprintf(stdout, "c_pcb_c_per_x: %f\tc_pcb_c_per_y: %f\n", p->c_pcb_c_per_x, p->c_pcb_c_per_y);
	fprintf(stdout, "r_pcb: %f\tr_pcb_per: %f\n", p->r_pcb, p->r_pcb_per);
	fprintf(stdout, "c_pcb_per: %f\n", p->c_pcb_per);
	fprintf(stdout, "r_amb_sec_c_per_x: %f\tr_amb_sec_c_per_y: %f\n", p->r_amb_sec_c_per_x, p->r_amb_sec_c_per_y);
	fprintf(stdout, "c_amb_sec_c_per_x: %f\tc_amb_sec_c_per_y: %f\n", p->c_amb_sec_c_per_x, p->c_amb_sec_c_per_y);
	fprintf(stdout, "r_amb_sec_per: %f\n", p->r_amb_sec_per);
	fprintf(stdout, "c_amb_sec_per: %f\n", p->c_amb_sec_per);
}

/*
 * wrapper routines interfacing with those of the corresponding
 * thermal model (block or grid)
 */

/*
 * allocate memory for the matrices. for the block model, placeholder
 * can be an empty floorplan frame with only the names of the functional
 * units. for the grid model, it is the default floorplan
 */
RC_model_t *alloc_RC_model(thermal_config_t *config, flp_t *placeholder, microchannel_config_t *microchannel_config, materials_list_t *materials_list,
	                         int do_detailed_3D, int use_microchannels) //BU_3D: do_detailed_3D option added.
{
	RC_model_t *model= (RC_model_t *) calloc (1, sizeof(RC_model_t));
	if (!model)
		fatal("memory allocation error\n");
	if(!(strcasecmp(config->model_type, BLOCK_MODEL_STR))) {
		model->type = BLOCK_MODEL;
		model->block = alloc_block_model(config, placeholder);
		model->config = &model->block->config;
	} else if(!(strcasecmp(config->model_type, GRID_MODEL_STR))) {
		model->type = GRID_MODEL;
		model->grid = alloc_grid_model(config, placeholder, microchannel_config, materials_list, do_detailed_3D, use_microchannels);
		model->config = &model->grid->config;
	} else
		fatal("unknown model type\n");
	return model;
}

/* populate the thermal restistance values */
void populate_R_model(RC_model_t *model, flp_t *flp)
{
	if (model->type == BLOCK_MODEL)
		populate_R_model_block(model->block, flp);
	else if (model->type == GRID_MODEL)
		populate_R_model_grid(model->grid, flp);
	else fatal("unknown model type\n");
}

/* populate the thermal capacitance values */
void populate_C_model(RC_model_t *model, flp_t *flp)
{
	if (model->type == BLOCK_MODEL)
		populate_C_model_block(model->block, flp);
	else if (model->type == GRID_MODEL)
		populate_C_model_grid(model->grid, flp);
	else fatal("unknown model type\n");
}

/* steady state temperature	*/
void steady_state_temp(RC_model_t *model, double *power, double *temp)
{
	int leak_convg_true = 0;
	int leak_iter = 0;
	int n, base=0;
	//int idx=0;
	double blk_height, blk_width;
	int i, j, k;

	double *d_temp = NULL;
	double *temp_old = NULL;
	double *power_new = NULL;
	double d_max=0.0;

	if (model->type == BLOCK_MODEL) {
		n = model->block->flp->n_units;
		if (model->config->leakage_used) { // if considering leakage-temperature loop
			d_temp = hotspot_vector(model);
			temp_old = hotspot_vector(model);
			power_new = hotspot_vector(model);
			for (leak_iter=0;(!leak_convg_true)&&(leak_iter<=LEAKAGE_MAX_ITER);leak_iter++){
				for(i=0; i < n; i++) {
					blk_height = model->block->flp->units[i].height;
					blk_width = model->block->flp->units[i].width;
					power_new[i] = power[i] + calc_leakage(model->config->leakage_mode,blk_height,blk_width,temp[i]);
					temp_old[i] = temp[i]; //copy temp before update
				}
				steady_state_temp_block(model->block, power_new, temp); // update temperature
				d_max = 0.0;
				for(i=0; i < n; i++) {
					d_temp[i] = temp[i] - temp_old[i]; //temperature increase due to leakage
					if (d_temp[i]>d_max) {
						d_max = d_temp[i];
					}
				}
				if (d_max < LEAK_TOL) {// check convergence
					leak_convg_true = 1;
				}
				if (d_max > TEMP_HIGH && leak_iter > 1) {// check to make sure d_max is not "nan" (esp. in natural convection)
					fatal("temperature is too high, possible thermal runaway. Double-check power inputs and package settings.\n");
				}
			}
			free(d_temp);
			free(temp_old);
			free(power_new);
			/* if no convergence after max number of iterations, thermal runaway */
			if (!leak_convg_true)
				fatal("too many iterations before temperature-leakage convergence -- possible thermal runaway\n");
		} else // if leakage-temperature loop is not considered
			steady_state_temp_block(model->block, power, temp);
	}
	else if (model->type == GRID_MODEL)	{
		if (model->config->leakage_used) { // if considering leakage-temperature loop
			d_temp = hotspot_vector(model);
			temp_old = hotspot_vector(model);
			power_new = hotspot_vector(model);
			for (leak_iter=0;(!leak_convg_true)&&(leak_iter<=LEAKAGE_MAX_ITER);leak_iter++){
				for(k=0, base=0; k < model->grid->n_layers; k++) {
					if(model->grid->layers[k].has_power)
						for(j=0; j < model->grid->layers[k].flp->n_units; j++) {
							blk_height = model->grid->layers[k].flp->units[j].height;
							blk_width = model->grid->layers[k].flp->units[j].width;
							power_new[base+j] = power[base+j] + calc_leakage(model->config->leakage_mode,blk_height,blk_width,temp[base+j]);
							temp_old[base+j] = temp[base+j]; //copy temp before update
						}
					base += model->grid->layers[k].flp->n_units;
				}
				steady_state_temp_grid(model->grid, power_new, temp);
				d_max = 0.0;
				for(k=0, base=0; k < model->grid->n_layers; k++) {
					if(model->grid->layers[k].has_power)
						for(j=0; j < model->grid->layers[k].flp->n_units; j++) {
							d_temp[base+j] = temp[base+j] - temp_old[base+j]; //temperature increase due to leakage
							if (d_temp[base+j]>d_max)
								d_max = d_temp[base+j];
						}
					base += model->grid->layers[k].flp->n_units;
				}
				if (d_max < LEAK_TOL) {// check convergence
					leak_convg_true = 1;
				}
				if (d_max > TEMP_HIGH && leak_iter > 0) {// check to make sure d_max is not "nan" (esp. in natural convection)
					fatal("temperature is too high, possible thermal runaway. Double-check power inputs and package settings.\n");
				}
			}
			free(d_temp);
			free(temp_old);
			free(power_new);
			/* if no convergence after max number of iterations, thermal runaway */
			if (!leak_convg_true)
				fatal("too many iterations before temperature-leakage convergence -- possible thermal runaway\n");
		} else // if leakage-temperature loop is not considered
			steady_state_temp_grid(model->grid, power, temp);
	}
	else fatal("unknown model type\n");
}

/* transient (instantaneous) temperature	*/
void compute_temp(RC_model_t *model, double *power, double *temp, double time_elapsed)
{
	if (model->type == BLOCK_MODEL)
		compute_temp_block(model->block, power, temp, time_elapsed);
	else if (model->type == GRID_MODEL)
		compute_temp_grid(model->grid, power, temp, time_elapsed);
	else fatal("unknown model type\n");
}

/* differs from 'dvector()' in that memory for internal nodes is also allocated	*/
double *hotspot_vector(RC_model_t *model)
{
	if (model->type == BLOCK_MODEL)
		return hotspot_vector_block(model->block);
	else if (model->type == GRID_MODEL)
		return hotspot_vector_grid(model->grid);
	else fatal("unknown model type\n");
	return NULL;
}

/* copy 'src' to 'dst' except for a window of 'size'
 * elements starting at 'at'. useful in floorplan
 * compaction
 */
void trim_hotspot_vector(RC_model_t *model, double *dst, double *src,
						 int at, int size)
{
	if (model->type == BLOCK_MODEL)
		trim_hotspot_vector_block(model->block, dst, src, at, size);
	else if (model->type == GRID_MODEL)
		trim_hotspot_vector_grid(model->grid, dst, src, at, size);
	else fatal("unknown model type\n");
}

/* update the model's node count	*/
void resize_thermal_model(RC_model_t *model, int n_units)
{
	if (model->type == BLOCK_MODEL)
		resize_thermal_model_block(model->block, n_units);
	else if (model->type == GRID_MODEL)
		resize_thermal_model_grid(model->grid, n_units);
	else fatal("unknown model type\n");
}

/* sets the temperature of a vector 'temp' allocated using 'hotspot_vector'	*/
void set_temp(RC_model_t *model, double *temp, double val)
{
	if (model->type == BLOCK_MODEL)
		set_temp_block(model->block, temp, val);
	else if (model->type == GRID_MODEL)
		set_temp_grid(model->grid, temp, val);
	else fatal("unknown model type\n");
}

/* dump temperature vector alloced using 'hotspot_vector' to 'file' */
void dump_temp(RC_model_t *model, double *temp, char *file)
{
	if (model->type == BLOCK_MODEL)
		dump_temp_block(model->block, temp, file);
	else if (model->type == GRID_MODEL)
		dump_temp_grid(model->grid, temp, file);
	else fatal("unknown model type\n");
}

/* calculate average heatsink temperature for natural convection package */
double calc_sink_temp(RC_model_t *model, double *temp)
{
	if (model->type == BLOCK_MODEL)
		return calc_sink_temp_block(model->block, temp, model->config);
	else if (model->type == GRID_MODEL)
		return calc_sink_temp_grid(model->grid, temp, model->config);
	else fatal("unknown model type\n");
	return 0.0;
}

/* copy temperature vector from src to dst */
void copy_temp(RC_model_t *model, double *dst, double *src)
{
	if (model->type == BLOCK_MODEL)
		copy_temp_block(model->block, dst, src);
	else if (model->type == GRID_MODEL)
		copy_temp_grid(model->grid, dst, src);
	else fatal("unknown model type\n");
}

/*
 * read temperature vector alloced using 'hotspot_vector' from 'file'
 * which was dumped using 'dump_temp'. values are clipped to thermal
 * threshold based on 'clip'
 */
void read_temp(RC_model_t *model, double *temp, char *file, int clip)
{
	if (model->type == BLOCK_MODEL)
		read_temp_block(model->block, temp, file, clip);
	else if (model->type == GRID_MODEL)
		read_temp_grid(model->grid, temp, file, clip);
	else fatal("unknown model type\n");
}

/* dump power numbers to file	*/
void dump_power(RC_model_t *model, double *power, char *file)
{
	if (model->type == BLOCK_MODEL)
		dump_power_block(model->block, power, file);
	else if (model->type == GRID_MODEL)
		dump_power_grid(model->grid, power, file);
	else fatal("unknown model type\n");
}

/*
 * read power vector alloced using 'hotspot_vector' from 'file'
 * which was dumped using 'dump_power'.
 */
void read_power(RC_model_t *model, double *power, char *file)
{
	if (model->type == BLOCK_MODEL)
		read_power_block(model->block, power, file);
	else if (model->type == GRID_MODEL)
		read_power_grid(model->grid, power, file);
	else fatal("unknown model type\n");
}

/* peak temperature on chip	*/
double find_max_temp(RC_model_t *model, double *temp)
{
	if (model->type == BLOCK_MODEL)
		return find_max_temp_block(model->block, temp);
	else if (model->type == GRID_MODEL)
		return find_max_temp_grid(model->grid, temp);
	else fatal("unknown model type\n");
	return 0.0;
}

/* average temperature on chip	*/
double find_avg_temp(RC_model_t *model, double *temp)
{
	if (model->type == BLOCK_MODEL)
		return find_avg_temp_block(model->block, temp);
	else if (model->type == GRID_MODEL)
		return find_avg_temp_grid(model->grid, temp);
	else fatal("unknown model type\n");
	return 0.0;
}

/* debug print	*/
void debug_print_model(RC_model_t *model)
{
	if (model->type == BLOCK_MODEL)
		debug_print_block(model->block);
	else if (model->type == GRID_MODEL)
		debug_print_grid(model->grid);
	else fatal("unknown model type\n");
}

/* calculate temperature-dependent leakage power */
/* will support HotLeakage in future releases */
double calc_leakage(int mode, double h, double w, double temp)
{
	/* a simple leakage model.
	 * Be aware -- this model may not be accurate in some cases.
	 * You may want to use your own temperature-dependent leakage model here.
	 */
	double leak_alpha = 1.5e+4;
	double leak_beta = 0.036;
	double leak_Tbase = 383.15; /* 110C according to the above paper */

	double leakage_power;

	if (mode)
		fatal("HotLeakage currently is not implemented in this release of HotSpot, please check back later.\n");

	leakage_power = leak_alpha*h*w*exp(leak_beta*(temp-leak_Tbase));
	return leakage_power;
}

/* destructor */
void delete_RC_model(RC_model_t *model)
{
	if (model->type == BLOCK_MODEL)
		delete_block_model(model->block);
  else if (model->type == GRID_MODEL)
		delete_grid_model(model->grid);
	else fatal("unknown model type\n");
	free(model);
}

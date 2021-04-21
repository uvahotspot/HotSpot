#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "microchannel.h"

#if SUPERLU > 0
#include "slu_ddefs.h"
#endif

#define MAX_LINE_SIZE 4096
#define DEBUG 0

// How many extra nodes we need to include in the pressure circuit
int extra_pressure_nodes;

// default microchannel configuration parameters
microchannel_config_t default_microchannel_config(void)
{
  microchannel_config_t config;

  config.cell_width        = 100e-6;
  config.cell_height       = 100e-6;
  config.cell_thickness    = 100e-6;
  config.pumping_pressure  = 5000;
  config.pump_internal_res = 0;            // ideal pump
  config.inlet_temperature = 300;
  config.coolant_capac     = 4172638;      // water
  config.coolant_res       = 1.647717911;  // water
  config.coolant_visc      = 0.000889;     // water
  config.wall_capac        = 1635660;      // silicon
  config.wall_res          = 0.0076923077; // silicon
  config.htc               = 27132;
  config.num_rows          = -1;
  config.num_columns       = -1;
  config.n_fluid_cells     = -1;
  config.cell_types        = NULL;
  config.mapping           = NULL;
  config.A                 = NULL;
  config.b                 = NULL;
  config.nnz               = 0;

  return config;
}

/*
 * parse a table of name-value string pairs and add the configuration
 * parameters to 'config'
 */
void microchannel_config_add_from_strs(microchannel_config_t *config, materials_list_t *materials_list, str_pair *table, int size)
{
  int idx;

  if ((idx = get_str_index(table, size, "pumping_pressure")) >= 0)
    if(sscanf(table[idx].value, "%lf", &config->pumping_pressure) != 1)
      fatal("invalid format for configuration  parameter pumping_pressure\n");
  if ((idx = get_str_index(table, size, "pump_internal_res")) >= 0)
    if(sscanf(table[idx].value, "%lf", &config->pump_internal_res) != 1)
      fatal("invalid format for configuration  parameter pumping internal resistance\n");
  if ((idx = get_str_index(table, size, "inlet_temperature")) >= 0)
    if(sscanf(table[idx].value, "%lf", &config->inlet_temperature) != 1)
      fatal("invalid format for configuration  parameter inlet_temperature\n");
  if ((idx = get_str_index(table, size, "coolant_capac")) >= 0)
    if(sscanf(table[idx].value, "%lf", &config->coolant_capac) != 1)
      fatal("invalid format for configuration  parameter coolant_capacity\n");
  if ((idx = get_str_index(table, size, "coolant_res")) >= 0)
    if(sscanf(table[idx].value, "%lf", &config->coolant_res) != 1)
      fatal("invalid format for configuration  parameter coolant resistivity\n");
  if ((idx = get_str_index(table, size, "coolant_visc")) >= 0)
    if(sscanf(table[idx].value, "%lf", &config->coolant_visc) != 1)
      fatal("invalid format for configuration  parameter coolant dynamic viscosity\n");
  if ((idx = get_str_index(table, size, "wall_capac")) >= 0)
    if(sscanf(table[idx].value, "%lf", &config->wall_capac) != 1)
      fatal("invalid format for configuration  parameter wall_capacity\n");
  if ((idx = get_str_index(table, size, "wall_res")) >= 0)
    if(sscanf(table[idx].value, "%lf", &config->wall_res) != 1)
      fatal("invalid format for configuration  parameter wall_resistivity\n");
  if ((idx = get_str_index(table, size, "htc")) >= 0)
    if(sscanf(table[idx].value, "%lf", &config->htc) != 1)
      fatal("invalid format for configuration  parameter heat transfer coefficient\n");
  if ((idx = get_str_index(table, size, "network_file")) >= 0)
    if(sscanf(table[idx].value, "%s", config->network_file) != 1)
      fatal("invalid format for configuration  parameter network_file\n");

  if ((idx = get_str_index(table, size, "wall_material")) >= 0) {
    char material_name[STR_SIZE];
    if(sscanf(table[idx].value, "%s", material_name) != 1)
      fatal("invalid format for configuration parameter wall_material\n");

    config->wall_res = 1.0 / get_material_thermal_conductivity(materials_list, material_name);
    config->wall_capac = get_material_volumetric_heat_capacity(materials_list, material_name);

    if(config->wall_res < 0 || config->wall_capac < 0)
      fatal("material name specified in configuration parameter wall_material not found\n");
  }

  if ((idx = get_str_index(table, size, "coolant_material")) >= 0) {
    char material_name[STR_SIZE];
    if(sscanf(table[idx].value, "%s", material_name) != 1)
      fatal("invalid format for configuration parameter coolant_material\n");

    config->coolant_res = 1.0 / get_material_thermal_conductivity(materials_list, material_name);
    config->coolant_capac = get_material_volumetric_heat_capacity(materials_list, material_name);
    config->coolant_visc = get_material_dynamic_viscosity(materials_list, material_name);

    if(config->coolant_res < 0 || config->coolant_capac < 0 || config->coolant_visc < 0)
      fatal("material name specified in configuration parameter coolant_material not found\n");
  }
}

/*
 * convert config into a table of name-value pairs. returns the no.
 * of parameters converted
 */
int microchannel_config_to_strs(microchannel_config_t *config, str_pair *table, int max_entries)
{
  if (max_entries < 17)
    fatal("not enough entries in table\n");

  sprintf(table[0].name, "cell_width");
  sprintf(table[1].name, "cell_height");
  sprintf(table[2].name, "cell_thickness");
  sprintf(table[3].name, "pumping_pressure");
  sprintf(table[4].name, "pump_internal_res");
  sprintf(table[5].name, "inlet_temperature");
  sprintf(table[6].name, "coolant_capac");
  sprintf(table[7].name, "coolant_res");
  sprintf(table[8].name, "coolant_visc");
  sprintf(table[9].name, "wall_capac");
  sprintf(table[10].name, "wall_res");
  sprintf(table[11].name, "htc");
  sprintf(table[12].name, "network_file");
  sprintf(table[13].name, "floorplan_file");
  sprintf(table[14].name, "num_rows");
  sprintf(table[15].name, "num_columns");
  sprintf(table[16].name, "n_fluid_cells");

  sprintf(table[0].value, "%e", config->cell_width);
  sprintf(table[1].value, "%e", config->cell_height);
  sprintf(table[2].value, "%e", config->cell_thickness);
  sprintf(table[3].value, "%e", config->pumping_pressure);
  sprintf(table[4].value, "%e", config->pump_internal_res);
  sprintf(table[5].value, "%e", config->inlet_temperature);
  sprintf(table[6].value, "%e", config->coolant_capac);
  sprintf(table[7].value, "%e", config->coolant_res);
  sprintf(table[8].value, "%e", config->coolant_visc);
  sprintf(table[9].value, "%e", config->wall_capac);
  sprintf(table[10].value, "%e", config->wall_res);
  sprintf(table[11].value, "%e", config->htc);
  sprintf(table[12].value, "%s", config->network_file);
  sprintf(table[13].value, "%s", config->floorplan_file);
  sprintf(table[14].value, "%d", config->num_rows);
  sprintf(table[15].value, "%d", config->num_columns);
  sprintf(table[16].value, "%d", config->n_fluid_cells);

  return 16;
}

void solve_pressure_circuit(microchannel_config_t *config) {
#if SUPERLU > 0
  SuperMatrix A, L, U, B;
  double *a, *rhs;
  int *asub, *xa;
  int *perm_r;
  int *perm_c;
  int nnz, nrhs, info, i, m, n, perc_spec;
  int j;
  superlu_options_t options;
  SuperLUStat_t stat;

  m = n = config->n_fluid_cells + extra_pressure_nodes;
  nnz = config->nnz;
  if( !(a = doubleMalloc(nnz)) ) ABORT("malloc failed");
  if( !(asub = intMalloc(nnz)) ) ABORT("malloc failed");
  if( !(xa = intMalloc(n + 1)) ) ABORT("malloc failed");

  int v = 0;
  int x = 1;
  xa[0] = 0;
  for(i = 0; i < m; i++) {
    for(j = 0; j < n; j++) {
      if(config->A[i][j] != 0) {
        a[v] = config->A[i][j];
        asub[v++] = j;
      }
    }
    xa[x++] = v;
  }

  if(DEBUG) {
    fprintf(stderr, "\n");
    for(i = 0; i < nnz; i++)
      fprintf(stderr, "a[%d] = %.15lf\n", i, a[i]);

    for(i = 0; i < nnz; i++)
      fprintf(stderr, "asub[%d] = %d\n", i, asub[i]);

    for(i = 0; i < n+1; i++)
      fprintf(stderr, "xa[%d] = %d\n", i, xa[i]);
  }

  dCreate_CompCol_Matrix(&A, m, n, nnz, a, asub, xa, SLU_NR, SLU_D, SLU_GE);

  nrhs = 1;
  if( !(rhs = doubleMalloc(m * nrhs)) ) ABORT("malloc failed");

  for(i = 0; i < m; i++) {
    rhs[i] = config->b[i];
  }

  dCreate_Dense_Matrix(&B, m, nrhs, rhs, m, SLU_DN, SLU_D, SLU_GE);

  if( !(perm_r = intMalloc(m)) ) ABORT("malloc failed");
  if( !(perm_c = intMalloc(n)) ) ABORT("malloc failed");

  set_default_options(&options);
  options.ColPerm = NATURAL;

  StatInit(&stat);

  dgssv(&options, &A, perm_c, perm_r, &L, &U, &B, &stat, &info);

  if(DEBUG) {
    //dPrint_CompCol_Matrix("A", &A);
    //dPrint_Dense_Matrix("B", &B);
  }

  DNformat *Astore = (DNformat *) B.Store;
  double *dp = (double *) Astore->nzval;

    for(i = 0; i < n; i++) {
      config->b[i] = dp[i];
      if(DEBUG)
        fprintf(stderr, "config->b[%d] = %e\n", i, config->b[i]);
    }

  SUPERLU_FREE(rhs);
  SUPERLU_FREE(perm_r);
  SUPERLU_FREE(perm_c);
  Destroy_CompCol_Matrix(&A);
  Destroy_SuperMatrix_Store(&B);
  Destroy_SuperNode_Matrix(&L);
  Destroy_CompCol_Matrix(&U);
  StatFree(&stat);
#else

  gaussj(config->A, config->n_fluid_cells + extra_pressure_nodes, config->b);

  if(DEBUG) {
    int i;
    for(i = 0; i < config->n_fluid_cells + extra_pressure_nodes; i++)
      fprintf(stderr, "config->b[%d] = %e\n", i, config->b[i]);
  }

#endif
}

// Parse config's CSV file to build internal array of microchannel network
void microchannel_build_network(microchannel_config_t *config) {
  char line[MAX_LINE_SIZE], str[STR_SIZE];
  char *cell;
  int cell_type, i = 0, j = 0;
  FILE *fp = fopen(config->network_file, "r");

  if(DEBUG)
    fprintf(stderr, "network_file = %s\n", config->network_file);

  if(!fp) {
    strncpy(str, "Unable to open ", STR_SIZE);
    strncat(str, config->network_file, STR_SIZE - strlen(str));
    strncat(str, "\n", STR_SIZE - strlen(str));
    fatal(str);
  }

  int nr = config->num_rows;
  int nc = config->num_columns;

  if(DEBUG)
    fprintf(stderr, "num_rows: %d, num_cols: %d\n", nr, nc);

  config->cell_types = calloc(nr, sizeof(int *));

  if (config->cell_types == NULL) {
    fprintf(stderr, "ERROR: Couldn't allocate space for cell_types array");
  }

  for(i = 0; i < nr; i++) {
    config->cell_types[i] = calloc(nc, sizeof(int));

    if (config->cell_types[i] == NULL) {
      fprintf(stderr, "ERROR: Couldn't allocate space for cell_types[%d] array", i);
    }
  }

  // parse network file to build cell_types array
  i = j = 0;
  fgets(line, MAX_LINE_SIZE, fp);
  while(!feof(fp)) {
    cell = strtok(line, " ,\n");
    while(cell != NULL) {
      cell_type = atoi(cell);
      config->cell_types[i][j] = cell_type;
      if(DEBUG)
        fprintf(stderr, "Adding cell %d, %d; type = %d\n", i, j, config->cell_types[i][j]);
      j++;


      if(j > nc) {
        fclose(fp);
        sprintf(str, "Microchannel row %d has more cells than num_columns(%d)\n", i, nc);
        fatal(str);
      }

      cell = strtok(NULL, " ,\n");
    }

    i++; j = 0;

    if(i > nr) {
      fclose(fp);
      sprintf(str, "Microchannel has more rows than num_rows(%d)\n", nr);
      fatal(str);
    }

    fgets(line, MAX_LINE_SIZE, fp);
  }

  fclose(fp);

  // Create floorplan file for microchannel
  strcpy(config->floorplan_file, config->network_file);
  char *ext = strstr(config->floorplan_file, NETWORK_EXTENSION);
  strcpy(ext, FLOORPLAN_EXTENSION);
  fp = fopen(config->floorplan_file, "w");

  fprintf(fp, "# Name\tw\th\tx\ty\tc_v\tp\n");
  for(i = 0; i < nr; i++) {
    for(j = 0; j < nc; j++) {
      if(IS_FLUID_CELL(config, i, j)) {
        // For floorplan files, x = y = 0 is the bottom left corner, but for our parsing i = j = 0 is
        // the top left cell, so we have to account for that
        fprintf(fp, "Cell_%d_%d\t%e\t%e\t%e\t%e\t%e\t%e\n", i, j, config->cell_width, config->cell_height,
                j*config->cell_width, (config->num_rows - i - 1)*config->cell_height, config->coolant_capac,
                config->coolant_res);
      }
      else {
        fprintf(fp, "Cell_%d_%d\t%e\t%e\t%e\t%e\t%e\t%e\n", i, j, config->cell_width, config->cell_height,
                j*config->cell_width, (config->num_rows - i - 1)*config->cell_height, config->wall_capac,
                config->wall_res);
      }
    }
  }

  fclose(fp);

  printf("Creating pressure circuit...\n");
  build_pressure_matrix(config);
  printf("Solving pressure circuit...\n");
  solve_pressure_circuit(config);
}

double hydroC(microchannel_config_t *config) {
  double h = config->cell_thickness;
  double w = config->cell_width;
  double L = config->cell_height;
  double viscosity = config->coolant_visc;
  double ret_val;

  if(h == w)
    ret_val = (0.42229 * pow(h, 4)) / (12 * viscosity * L);
  else if(h > w)
    ret_val = ((1 - 0.63*(w / h)) * (pow(w, 3)) * (h)) / (12 * viscosity * L);
  else
    ret_val = ((1 - 0.63*(h / w)) * (pow(h, 3)) * (w)) / (12 * viscosity * L);

  return ret_val;
}

void build_pressure_matrix(microchannel_config_t *config) {
  int i, j;
  int nr = config->num_rows;
  int nc = config->num_columns;
  int **mapping;
  int n = 0;
  config->nnz = 0;

  mapping = calloc(nr, sizeof(int *));

  if(mapping == NULL)
    fatal("Unable to allocate pressure circuit mapping\n");

  for(i = 0; i < nr; i++) {
    mapping[i] = calloc(nc, sizeof(int));

    if(mapping[i] == NULL) {
      fatal("Unable to allocate pressure circuit mapping\n");
    }
  }

  // Assign unique number to each fluid cell
  for(i = 0; i < nr; i++) {
    for(j = 0; j < nc; j++) {
      if(IS_FLUID_CELL(config, i, j)) {
        mapping[i][j] = n++;
      }
      else
        mapping[i][j] = -1;
    }
  }
  config->n_fluid_cells = n;
  config->mapping = mapping;

  // If we're modeling a non-ideal pump, we include one extra node for the pump
  if(config->pump_internal_res == 0) {
    extra_pressure_nodes = 0;
  }
  else {
    extra_pressure_nodes = 1;
  }

  config->A = calloc(config->n_fluid_cells + extra_pressure_nodes, sizeof(double *));

  if(!config->A)
    fatal("Unable to allocate matrix A for pressure circuit\n");
  for(i = 0; i < config->n_fluid_cells + extra_pressure_nodes; i++) {
    config->A[i] = calloc(config->n_fluid_cells, sizeof(double));

    if(!config->A[i])
      fatal("Unable to allocate matrix A for pressure circuit\n");
  }

  config->b = calloc(config->n_fluid_cells + extra_pressure_nodes, sizeof(double));

  if(!config->b)
    fatal("Unable to allocate matrix b for pressure circuit\n");


  if(DEBUG) {
    fprintf(stderr, "Mapping number: %d\n", n);
    for(i = 0; i < nr; i++) {
      for(j = 0; j < nc; j++) {
        fprintf(stderr, "mapping[%d][%d] = %d\n", i, j, mapping[i][j]);
      }
    }
  }

  // Iterate through all cells
  double diagonal_val = 0;
  double hydro_conductance = -hydroC(config);
  for(i = 0; i < nr; i++) {
    for(j = 0; j < nc; j++) {
      if(config->cell_types[i][j] == FLUID ||
         (config->cell_types[i][j] == INLET && config->pump_internal_res != 0)) {
        // northern cell
          if(i > 0 && IS_FLUID_CELL(config, i-1, j)) {
            if(DEBUG)
              fprintf(stderr, "[%d, %d]: Northern cell. Setting A[%d][%d] = %.15lf\n", i, j, mapping[i][j], mapping[i-1][j], hydro_conductance);

            config->A[mapping[i][j]][mapping[i-1][j]] = hydro_conductance;
            diagonal_val += hydro_conductance;
            config->nnz++;
          }

        // southern cell
          if(i < nr - 1 && IS_FLUID_CELL(config, i+1, j)) {
            if(DEBUG)
              fprintf(stderr, "[%d, %d]: Southern cell. Setting A[%d][%d] = %.15lf\n", i, j, mapping[i][j], mapping[i+1][j], hydro_conductance);

            config->A[mapping[i][j]][mapping[i+1][j]] = hydro_conductance;
            diagonal_val += hydro_conductance;
            config->nnz++;
          }

        // western cell
          if(j > 0 && IS_FLUID_CELL(config, i, j-1)) {
            if(DEBUG)
              fprintf(stderr, "[%d, %d]: Western Cell. Setting A[%d][%d] = %.15lf\n", i, j, mapping[i][j], mapping[i][j-1], hydro_conductance);

            config->A[mapping[i][j]][mapping[i][j-1]] = hydro_conductance;
            diagonal_val += hydro_conductance;
            config->nnz++;
          }

        // eastern cell
          if(j < nc - 1 && IS_FLUID_CELL(config, i, j+1)) {
            if(DEBUG)
              fprintf(stderr, "[%d, %d]: Eastern Cell. Setting A[%d][%d] = %.15lf\n", i, j, mapping[i][j], mapping[i][j+1], hydro_conductance);

            config->A[mapping[i][j]][mapping[i][j+1]] = hydro_conductance;
            diagonal_val += hydro_conductance;
            config->nnz++;
          }

        // diagonal
        if(DEBUG)
          fprintf(stderr, "[%d, %d]: Diagonal. Setting A[%d][%d] = %.15lf\n", i, j, mapping[i][j], mapping[i][j], -diagonal_val);

        config->A[mapping[i][j]][mapping[i][j]] = -diagonal_val;
        config->nnz++;
        diagonal_val = 0;
      }

      if(IS_INLET_CELL(config, i, j)) {
        // Non-ideal pump
        if(config->pump_internal_res != 0) {
          if(DEBUG) {
            fprintf(stderr, "[%d, %d]: Inlet. Setting A[%d][%d] = %e\n", i, j, mapping[i][j], config->n_fluid_cells, -1.0 / config->pump_internal_res);
            fprintf(stderr, "[%d, %d]: Inlet. Setting A[%d][%d] = %e\n", i, j, mapping[i][j], mapping[i][j], 1.0 / config->pump_internal_res);
          }

          // Inlet cells are connected to the pump through the pump's internal
          // resistance
          config->A[mapping[i][j]][config->n_fluid_cells] = -1.0 / config->pump_internal_res;
          config->A[mapping[i][j]][mapping[i][j]] += 1.0 / config->pump_internal_res;
          config->nnz++; // diagonal val already exists, so we're only adding one nonzero value
        }

        // Ideal pump
        else {
          if(DEBUG) {
            fprintf(stderr, "[%d, %d]: Inlet. Setting A[%d][%d] = %e\n", i, j, mapping[i][j], mapping[i][j], 1.0);
            fprintf(stderr, "[%d, %d]: Inlet. Setting b[%d] = %e\n", i, j, mapping[i][j], config->pumping_pressure);
          }
          config->A[mapping[i][j]][mapping[i][j]] = 1.0;
          config->b[mapping[i][j]] = config->pumping_pressure;
          config->nnz++;
        }
      }

      else if(IS_OUTLET_CELL(config, i, j)) {
        if(DEBUG)
          fprintf(stderr, "[%d, %d]: OUTLET. Setting A[%d][%d] = 1\n", i, j, mapping[i][j], mapping[i][j]);

        config->A[mapping[i][j]][mapping[i][j]] = 1;
        config->nnz++;
      }
    }
  }

  if(config->pump_internal_res != 0) {
    if(DEBUG) {
      fprintf(stderr, "Pump Node. Setting A[%d][%d] = %e\n", config->n_fluid_cells, config->n_fluid_cells, 1.0);
      fprintf(stderr, "Pump Node. Setting b[%d] = %e\n", config->n_fluid_cells, config->pumping_pressure);
    }

    // Handle extra node representing pump
    config->A[config->n_fluid_cells][config->n_fluid_cells] = 1;
    config->b[config->n_fluid_cells] = config->pumping_pressure;
    config->nnz++;
  }

  if(DEBUG) {
    fprintf(stderr, "Nonzero values (%d total):\n", config->nnz);
     for(i = 0; i < config->n_fluid_cells + extra_pressure_nodes; i++) {
      for(j = 0; j < config->n_fluid_cells + extra_pressure_nodes; j++) {
        if(config->A[i][j] != 0)
          fprintf(stderr, "A[%d][%d] = %e\n", i, j, config->A[i][j]);
      }
    }

    if(DEBUG)
      for(i = 0; i < config->n_fluid_cells + extra_pressure_nodes; i++)
        fprintf(stderr, "b[%d] = %e\n", i, config->b[i]);

  }
}

double flow_rate(microchannel_config_t * config, int cell1_i, int cell1_j, int cell2_i, int cell2_j) {
  double *pressure = config->b;
  int **mapping = config->mapping;
  return (pressure[mapping[cell1_i][cell1_j]] - pressure[mapping[cell2_i][cell2_j]]) * hydroC(config);
}

// Copy user-defined parameters from one microchannel config to another
void copy_microchannel(microchannel_config_t *dst, microchannel_config_t *src) {
  dst->pumping_pressure  = src->pumping_pressure;
  dst->pump_internal_res = src->pump_internal_res;
  dst->inlet_temperature = src->inlet_temperature;
  dst->coolant_capac     = src->coolant_capac;
  dst->coolant_res       = src->coolant_res;
  dst->coolant_visc      = src->coolant_visc;
  dst->wall_capac        = src->wall_capac;
  dst->wall_res          = src->wall_res;
  dst->htc               = src->htc;
}

void free_microchannel(microchannel_config_t *config) {
  int i;
  if(config) {
    if(config->cell_types) {
      for(i = 0; i < config->num_rows; i++) {
        free(config->cell_types[i]);
      }
      free(config->cell_types);
    }

    if(config->A) {
      for(i = 0; i < config->n_fluid_cells; i++) {
        free(config->A[i]);
      }
      free(config->A);
    }

    if(config->b) {
      free(config->b);
    }

    if(config->mapping) {
      for(i = 0; i < config->num_rows; i++) {
        free(config->mapping[i]);
      }
      free(config->mapping);
    }

    free(config);
  }
}

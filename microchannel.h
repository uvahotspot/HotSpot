#ifndef __MICROCHANNEL_H
#define __MICROCHANNEL_H

#include "util.h"
#include "materials.h"

// file extensions used for microchannels
#define NETWORK_EXTENSION ".csv"
#define FLOORPLAN_EXTENSION ".flp"

// Different types of cells in microchannel layer
#define TSV -1
#define SOLID 0
#define FLUID 1
#define INLET 2
#define OUTLET 3

// One extra node for the pump
#define EXTRA_PRESSURE_NODES 1

#define IS_FLUID_CELL(uconf, i, j)  (uconf->cell_types[i][j] == FLUID || \
                                     uconf->cell_types[i][j] == INLET || \
                                     uconf->cell_types[i][j] == OUTLET)

#define IS_INLET_CELL(uconf, i, j)  (uconf->cell_types[i][j] == INLET)

#define IS_OUTLET_CELL(uconf, i, j) (uconf->cell_types[i][j] == OUTLET)

typedef struct microchannel_config_t_st
{
  // width of an individual cell in m
  double cell_width;

  // height of an individual cell in m
  double cell_height;

  // thickness of an individual cell in m
  double cell_thickness;

  // pumping pressure between inlet(s) and outlet(s) in Pa
  double pumping_pressure;

  // pump internal resistance in K/W
  double pump_internal_res;

  // temperature of coolant at inlet in K
  double inlet_temperature;

  // volumetric heat capacity of coolant in J/m^3-K
  double coolant_capac;

  // resistivity of coolant in m-K/W
  double coolant_res;

  // dynamic viscosity of coolant in Pa-s/m^3
  double coolant_visc;

  // volumetric heat capacity of channel walls in J/m^3-K
  double wall_capac;

  // resistivity of channel walls in m-K/W
  double wall_res;

  // Heat Transfer Coefficient from coolant to channel walls in W/m^2-K
  double htc;

  // csv file containing description of microchannel network
  char network_file[STR_SIZE];

  // floorplan file created from network_file
  char floorplan_file[STR_SIZE+3];

  // rows and columns in microchannel network
  int num_rows;
  int num_columns;

  // number of fluid cells in microchannel network
  int n_fluid_cells;

  // array of cell types
  int **cell_types;

  // mapping from cell indices to pressure circuit node indices
  int **mapping;

  // For SuperLU
  double **A;
  double *b;
  int nnz;
} microchannel_config_t;

microchannel_config_t default_microchannel_config(void);
void microchannel_config_add_from_strs(microchannel_config_t *config, materials_list_t *materials_list, str_pair *table, int size);
int microchannel_config_to_strs(microchannel_config_t *config, str_pair *table, int max_entries);
void microchannel_build_network(microchannel_config_t *config);
void build_pressure_matrix(microchannel_config_t *config);
double flow_rate(microchannel_config_t *config, int cell1_i, int cell1_j, int cell2_i, int cell2_j);
void copy_microchannel(microchannel_config_t *src, microchannel_config_t *dst);
void free_microchannel(microchannel_config_t * config);

#endif

#ifndef MATERIALS_H
#define MATERIALS_H

#include "util.h"

#define SOLID_MATERIAL 0
#define FLUID_MATERIAL 1

typedef struct material_t_st
{
  int material_type; // solid or fluid
  double thermal_conductivity;
  double volumetric_heat_capacity;
  double dynamic_viscosity;
} material_t;

typedef struct materials_list_t_st
{
  // number of materials entries
  int size;

  // names of all materials
  char **names;

  // properties of each material
  material_t *material_properties;
} materials_list_t;

void default_materials(materials_list_t *materials_list);
void materials_add_from_file(materials_list_t *materials_list, char *materials_filename);
void free_materials(materials_list_t *materials_list);
material_t get_material_properties(materials_list_t *materials_list, char *name);
double get_material_thermal_conductivity(materials_list_t *materials_list, char *name);
double get_material_volumetric_heat_capacity(materials_list_t *materials_list, char *name);
double get_material_dynamic_viscosity(materials_list_t *materials_list, char *name);

#endif

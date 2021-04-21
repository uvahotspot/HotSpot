#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "materials.h"
#include "util.h"

#define MATERIAL_NAME 0
#define MATERIAL_TYPE 1
#define MATERIAL_THERMAL_CONDUCTIVITY 2
#define MATERIAL_VOLUMETRIC_HEAT_CAPACITY 3
#define MATERIAL_DYNAMIC_VISCOSITY 4

void default_materials(materials_list_t *materials_list) {
  materials_list->size = 0;
}

void materials_add_from_file(materials_list_t *materials_list, char *materials_filename) {
  FILE *fp = fopen(materials_filename, "r");

  if(!fp) {
    char err_message[STR_SIZE];
    sprintf(err_message, "Unable to open file %s\n", materials_filename);
    fatal(err_message);
  }
  char line[LINE_SIZE], current_type[LINE_SIZE], sval[STR_SIZE];
  int field = MATERIAL_NAME, number_materials = 0;
  double dval;

  // first pass: count number of entries
  while(!feof(fp)) {
    fgets(line, LINE_SIZE, fp);
    if(feof(fp))
      break;

    // skip comments and empty lines
    char *ptr = strtok(line, " \r\t\n");
    if (!ptr || ptr[0] == '#')
      continue;

    switch(field) {
      case MATERIAL_NAME:
        field = MATERIAL_TYPE;
        break;
      case MATERIAL_TYPE:
        field = MATERIAL_THERMAL_CONDUCTIVITY;
        if(sscanf(ptr, "%s", current_type) != 1)
          fatal("invalid material type: must be solid or fluid\n");
        break;
      case MATERIAL_THERMAL_CONDUCTIVITY:
        field = MATERIAL_VOLUMETRIC_HEAT_CAPACITY;
        break;
      case MATERIAL_VOLUMETRIC_HEAT_CAPACITY:
        if(!strncmp(current_type, "solid", STR_SIZE)) {
          number_materials++;
          field = MATERIAL_NAME;
        }
        else if(!strncmp(current_type, "fluid", STR_SIZE)) {
          field = MATERIAL_DYNAMIC_VISCOSITY;
        }
        else {
          fatal("invalid material type: must be solid or fluid\n");
        }
        break;
      case MATERIAL_DYNAMIC_VISCOSITY:
        number_materials++;
        field = MATERIAL_NAME;
        break;
      default:
        fatal("ill-formed materials file\n");
    }
  }

  materials_list->size = number_materials;
  materials_list->names = calloc(number_materials, sizeof(char *));
  materials_list->material_properties = calloc(number_materials, sizeof(material_t));
  fseek(fp, 0, SEEK_SET);

  // second pass: fill in material properties
  field = MATERIAL_NAME;
  int i = 0;
  while(!feof(fp)) {
    fgets(line, LINE_SIZE, fp);
    if(feof(fp))
      break;

    // skip comments and empty lines
    char *ptr = strtok(line, " \r\t\n");
    if (!ptr || ptr[0] == '#')
      continue;

    switch(field) {
      case MATERIAL_NAME:
        if(sscanf(ptr, "%s", sval) != 1)
          fatal("invalid material name\n");
        materials_list->names[i] = calloc(1, STR_SIZE);
        strncpy(materials_list->names[i], sval, STR_SIZE);
        field = MATERIAL_TYPE;
        break;
      case MATERIAL_TYPE:
        if(sscanf(ptr, "%s", current_type) != 1)
          fatal("invalid material type: must be solid or fluid\n");
        if(!strncmp(current_type, "solid", STR_SIZE))
          materials_list->material_properties[i].material_type = SOLID_MATERIAL;
        else if(!strncmp(current_type, "fluid", STR_SIZE))
          materials_list->material_properties[i].material_type = FLUID_MATERIAL;
        else
          fatal("invalid material type: must be solid or fluid\n");
        field = MATERIAL_THERMAL_CONDUCTIVITY;
        break;
      case MATERIAL_THERMAL_CONDUCTIVITY:
        if(sscanf(ptr, "%lf", &dval) != 1)
          fatal("invalid material thermal conductivity\n");
        if(dval < 0)
          fatal("thermal conductivities must be nonnegative\n");
        materials_list->material_properties[i].thermal_conductivity = dval;
        field = MATERIAL_VOLUMETRIC_HEAT_CAPACITY;
        break;
      case MATERIAL_VOLUMETRIC_HEAT_CAPACITY:
        if(sscanf(ptr, "%lf", &dval) != 1)
          fatal("invalid material volumetric heat capacity\n");
        if(dval < 0)
          fatal("volumetric heat capacities must be nonnegative\n");
        materials_list->material_properties[i].volumetric_heat_capacity = dval;
        if(!strncmp(current_type, "solid", STR_SIZE)) {
          materials_list->material_properties[i].dynamic_viscosity = 0.0;
          i++;
          field = MATERIAL_NAME;
        }
        else if(!strncmp(current_type, "fluid", STR_SIZE)) {
          field = MATERIAL_DYNAMIC_VISCOSITY;
        }
        else {
          fatal("invalid material type: must be solid or fluid\n");
        }
        break;
      case MATERIAL_DYNAMIC_VISCOSITY:
        if(sscanf(ptr, "%lf", &dval) != 1)
          fatal("invalid material dynamic viscosity\n");
        if(dval < 0)
          fatal("dynamic viscosities must be nonnegative\n");
        materials_list->material_properties[i].dynamic_viscosity = dval;
        i++;
        field = MATERIAL_NAME;
        break;
      default:
        fatal("ill-formed materials file\n");
    }
  }
  fclose(fp);
}

void free_materials(materials_list_t *materials_list) {
  int i;
  for(i = 0; i < materials_list->size; i++) {
    free(materials_list->names[i]);
  }

  free(materials_list->names);
  free(materials_list->material_properties);
}

material_t get_material_properties(materials_list_t *materials_list, char *name) {
  int i;
  for(i = 0; i < materials_list->size; i++) {
    if(!strncmp(materials_list->names[i], name, STR_SIZE)) {
      return materials_list->material_properties[i];
    }
  }

  char err_message[STR_SIZE];
  strncpy(err_message, "Unable to find material properties for ", STR_SIZE);
  strncat(err_message, name, STR_SIZE - strlen(err_message));
  strncat(err_message, "\n", STR_SIZE - strlen(err_message));
  fatal(err_message);
}

double get_material_thermal_conductivity(materials_list_t *materials_list, char *name) {
  int i;
  for(i = 0; i < materials_list->size; i++) {
    if(!strncmp(materials_list->names[i], name, STR_SIZE)) {
      return materials_list->material_properties[i].thermal_conductivity;
    }
  }
  return -1;
}

double get_material_volumetric_heat_capacity(materials_list_t *materials_list, char *name) {
  int i;
  for(i = 0; i < materials_list->size; i++) {
    if(!strncmp(materials_list->names[i], name, STR_SIZE)) {
      return materials_list->material_properties[i].volumetric_heat_capacity;
    }
  }
  return -1;
}

double get_material_dynamic_viscosity(materials_list_t *materials_list, char *name) {
  int i;
  for(i = 0; i < materials_list->size; i++) {
    if(!strncmp(materials_list->names[i], name, STR_SIZE)) {
      return materials_list->material_properties[i].dynamic_viscosity;
    }
  }
  return -1;
}

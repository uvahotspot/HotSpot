#ifndef __SHAPE_H_
#define __SHAPE_H_

#include "flp.h"
#include "npe.h"

#define	MAX_STACK		128

/* piecewise linear shape curve	*/
typedef struct shape_t_st
{
  double *x;	/* width	*/
  double *y;	/* height	*/
  /* 
   * if this is not a leaf in an flp,
   * the position of the orientations
   * in the contributing shape curve
   */
  int *left_pos;
  int *right_pos;
  /*
   * the position of the demarcation
   * in the combined shape curve
   */
  double *median;
  int size;
}shape_t;

/* slicing tree node	*/
typedef struct tree_node_t_st
{
	shape_t *curve;
	union {
		int cut_type;
		int unit;
	}label;	
	struct tree_node_t_st *left;
	struct tree_node_t_st *right;
}tree_node_t;

/* tree node stack	*/
typedef struct tree_node_stack_t_st
{
	tree_node_t *array[MAX_STACK];
	int top;
}tree_node_stack_t;

/* 
 * make a shape curve from area and aspect ratios
 * - allocates memory internally
 */
shape_t *shape_from_aspect(double area, double min, 
                           double max, int rotable,
						   int n_orients);
/* shape curve arithmetic	*/
shape_t *shape_add(shape_t *shape1, shape_t *shape2, int cut_type);
/* uninitialization	*/
void free_shape(shape_t *shape); 
/* position in the shape curve with minimum area	*/
int min_area_pos(shape_t *curve);
/* debug print	*/
void print_shape(shape_t *shape);

/* stack operations	*/
tree_node_stack_t *new_tree_node_stack(void);
void tree_node_stack_push(tree_node_stack_t *stack, tree_node_t *node);
tree_node_t *tree_node_stack_pop(tree_node_stack_t *stack);
int tree_node_stack_isfull(tree_node_stack_t *stack);
int tree_node_stack_isempty(tree_node_stack_t *stack);
void free_tree_node_stack(tree_node_stack_t *stack);
void tree_node_stack_clear(tree_node_stack_t *stack);

/* tree operations	*/
/* construct floorplan slicing tree from NPE	*/
tree_node_t *tree_from_NPE(flp_desc_t *flp_desc, 
						   tree_node_stack_t *stack, 
						   NPE_t *expr);
void free_tree(tree_node_t *root);
void print_tree(tree_node_t *root, flp_desc_t *flp_desc);
/* 
 * print only the portion of the shape curves 
 * corresponding to the `pos'th entry of root->curve
 */
void print_tree_relevant(tree_node_t *root, int pos, flp_desc_t *flp_desc);
/* 
 * convert slicing tree into actual floorplan
 * returns the number of dead blocks compacted
 */
int tree_to_flp(tree_node_t *root, flp_t *flp, int compact_dead, 
				double compact_ratio);
#endif

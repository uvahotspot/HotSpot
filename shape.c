#include <stdio.h>
#ifdef _MSC_VER
#define strcasecmp    _stricmp
#define strncasecmp   _strnicmp
#else
#include <strings.h>
#endif
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "shape.h"
#include "flp.h"
#include "util.h"

/* make a shape curve from area and aspect ratios	*/
shape_t *shape_from_aspect(double area, double min, 
                           double max, int rotable,
						   int n_orients)
{
	shape_t *shape;
	double r=1, minx, maxx, tmin, tmax;
	int i, overlap = 0;

	/* rotatable blocks have 2n no. of orients	*/
	if (n_orients <= 1 || (n_orients & 1))
		fatal("n_orients should be an even number greater than 1\n");

	shape = (shape_t *) calloc(1, sizeof(shape_t));
	if (!shape)
		fatal("memory allocation error\n");
	if (min == max)	
		shape->size = 1 + (!!rotable);
	else
		shape->size = n_orients;

	shape->x = (double *) calloc(shape->size, sizeof(double));
	shape->y = (double *) calloc(shape->size, sizeof(double));
	if (!shape->x || !shape->y)
		fatal("memory allocation error\n");

	/* overlapping regions of aspect ratios */
	if (rotable && min <= 1.0 && max >= 1.0) {
		overlap = 1;
		tmin = MIN4(min, max, 1.0/min, 1.0/max);
		tmax = MAX4(min, max, 1.0/min, 1.0/max);
		min = tmin;
		max = tmax;
	}

	if (!rotable || overlap) {
		minx = sqrt(area * min);
		maxx = sqrt(area * max);
		if (shape->size > 1)
			r = pow((maxx / minx) , 1.0/(shape->size-1));
		for (i = 0; i < shape->size; i++) {
			shape->x[i] = minx * pow(r, i);
			shape->y[i] = area / shape->x[i];
		}
	/* rotable but no overlap, hence two sets of orientations	*/		
	} else {
		int n = shape->size / 2;
		/* orientations with aspect ratios < 1	*/
		tmin = MIN(min, 1.0/min);
		tmax = MIN(max, 1.0/max);
		minx = sqrt(area * MIN(tmin, tmax));
		maxx = sqrt(area * MAX(tmin, tmax));
		if ( n > 1)
			r = pow((maxx / minx) , 1.0/(n-1));
		for (i = 0; i < n; i++) {
			shape->x[i] = minx * pow(r, i);
			shape->y[i] = area / shape->x[i];
		}
		/* orientations with aspect ratios > 1	*/
		tmin = MAX(min, 1.0/min);
		tmax = MAX(max, 1.0/max);
		minx = sqrt(area * MIN(tmin, tmax));
		maxx = sqrt(area * MAX(tmin, tmax));
		if ( n > 1)
			r = pow((maxx / minx) , 1.0/(n-1));
		for (i = 0; i < n; i++) {
			shape->x[n+i] = minx * pow(r, i);
			shape->y[n+i] = area / shape->x[n+i];
		}
	}
	return shape;
}

void free_shape(shape_t *shape)
{
	free (shape->x);
	free (shape->y);
	if(shape->left_pos) {
		free(shape->left_pos);
		free(shape->right_pos);
		free(shape->median);
	}	
	free(shape);
}

void print_shape_entry(shape_t *shape, int i)
{
	fprintf(stdout, "%g\t%g", shape->x[i], shape->y[i]);
	if (shape->left_pos)
		fprintf(stdout, "\t%d\t%d\t%g", shape->left_pos[i], 
		        shape->right_pos[i], shape->median[i]);
	fprintf(stdout, "\n");
}

/* debug print  */
void print_shape(shape_t *shape)
{
	int i;
	if (!shape) {
		fprintf(stdout, "printing shape curve: NULL\n");
		return;
	}
	fprintf(stdout, "printing shape curve with %d elements\n", shape->size);
	for (i=0; i < shape->size; i++) 
		print_shape_entry(shape, i);
	fprintf(stdout, "\n");	
}

/* shape curve arithmetic	*/
shape_t *shape_add(shape_t *shape1, shape_t *shape2, int cut_type)
{
	int i=0, j=0, k=0, total=0, m, n;
	shape_t *sum;

	sum = (shape_t *) calloc(1, sizeof(shape_t));
	if (!sum)
		fatal("memory allocation error\n");

	/* shortcuts	*/	
	m = shape1->size;
	n = shape2->size;

	/* determine result size	*/
	while(i < m && j < n) {
		if (cut_type == CUT_VERTICAL) {
			if (shape1->y[i] >= shape2->y[j])
				i++;
			else
				j++;
		} else {
			if (shape1->x[m-1-i] >= shape2->x[n-1-j])
				i++;
			else
				j++;
		}
		total++;
	}

	sum->x = (double *) calloc(total, sizeof(double));
	sum->y = (double *) calloc(total, sizeof(double));
	sum->left_pos = (int *) calloc(total, sizeof(int));
	sum->right_pos = (int *) calloc(total, sizeof(int));
	sum->median = (double *) calloc(total, sizeof(double));
	if (!sum->x || !sum->y || !sum->left_pos || 
	    !sum->right_pos || !sum->median)
		fatal("memory allocation error\n");
	sum->size = total;	

	i=j=0;
	while(i < m && j < n) {
		/* vertical add	*/
		if (cut_type == CUT_VERTICAL) {
			sum->x[k] = shape1->x[i] + shape2->x[j];
			sum->y[k] = MAX(shape1->y[i], shape2->y[j]);
			sum->left_pos[k] = i;
			sum->right_pos[k] = j;
			sum->median[k] = shape1->x[i];
			if (shape1->y[i] >= shape2->y[j])
				i++;
			else
				j++;
		/* horizontal add */
		} else {
			/* 
			 * direction of increasing 'y' is the reverse of the
			 * regular direction of the shape curve. hence reverse
			 * the curve before adding
			 */
			sum->x[total-1-k] = MAX(shape1->x[m-1-i], shape2->x[n-1-j]);
			sum->y[total-1-k] =  shape1->y[m-1-i] + shape2->y[n-1-j];
			sum->left_pos[total-1-k] = m-1-i;
			sum->right_pos[total-1-k] = n-1-j;
			sum->median[total-1-k] = shape1->y[m-1-i];
			if (shape1->x[m-1-i] >= shape2->x[n-1-j])
				i++;
			else
				j++;
		}
		k++;
	}

	return sum;
}

/* copy a shape curve	*/
shape_t *shape_duplicate(shape_t *shape)
{
	shape_t *copy;
	int i;
	copy = (shape_t *) calloc(1, sizeof(shape_t));
	if (!copy)
		fatal("memory allocation error\n");
	copy->size = shape->size;
	copy->x = (double *) calloc(copy->size, sizeof(double));
	copy->y = (double *) calloc(copy->size, sizeof(double));
	if (!copy->x || !copy->y)
		fatal("memory allocation error\n");
	if (shape->left_pos) {
		copy->left_pos = (int *) calloc(copy->size, sizeof(int));
		copy->right_pos = (int *) calloc(copy->size, sizeof(int));
		copy->median = (double *) calloc(copy->size, sizeof(double));
		if(!copy->left_pos || !copy->right_pos || !copy->median)
			fatal("memory allocation error\n");
	}
	for(i=0; i < shape->size; i++) {
		copy->x[i] = shape->x[i];
		copy->y[i] = shape->y[i];
		if (shape->left_pos) {
			copy->left_pos[i] = shape->left_pos[i];
			copy->right_pos[i] = shape->right_pos[i];
			copy->median[i] = shape->median[i];
		}
	}
	return copy;
}

/* tree node stack operations	*/

/* constructor	*/
tree_node_stack_t *new_tree_node_stack(void)
{
	tree_node_stack_t *stack;
	stack = (tree_node_stack_t *)calloc(1, sizeof(tree_node_stack_t));
	if (!stack)
		fatal("memory allocation error\n");
	/* first free location	*/	
	stack->top = 0;
	return stack;
}

/* destructor	*/
void free_tree_node_stack(tree_node_stack_t *stack)
{
	free(stack);
}

/* is empty?	*/
int tree_node_stack_isempty(tree_node_stack_t *stack)
{
	if (stack->top <= 0)
		return TRUE;
	return FALSE;
}

/* is full?	*/
int tree_node_stack_isfull(tree_node_stack_t *stack)
{
	if (stack->top >= MAX_STACK)
		return TRUE;
	return FALSE;
}

/* push	*/
void tree_node_stack_push(tree_node_stack_t *stack, tree_node_t *node)
{
	if (!tree_node_stack_isfull(stack)) {
		stack->array[stack->top] = node;
		stack->top++;
	} else
		fatal("attempting to push into an already full stack\n");
}

/* pop	*/
tree_node_t *tree_node_stack_pop(tree_node_stack_t *stack)
{
	if (!tree_node_stack_isempty(stack)) {
		stack->top--;
		return stack->array[stack->top];
	} else {
		fatal("attempting to pop from an already empty stack\n");
		return NULL;
	}	
}

/* clear	*/
void tree_node_stack_clear(tree_node_stack_t *stack)
{
	stack->top = 0;
}

/* slicing tree routines	*/

/* construct floorplan slicing tree from NPE	*/
tree_node_t *tree_from_NPE(flp_desc_t *flp_desc,
						   tree_node_stack_t *stack,
						   NPE_t *expr)
{
	int i;
	tree_node_t *node = NULL, *left, *right;
	
	if (!tree_node_stack_isempty(stack))
		fatal("stack not empty\n");

	for (i=0; i < expr->size; i++) {
		node = (tree_node_t *) calloc(1, sizeof(tree_node_t));
		if (!node)
			fatal("memory allocation error\n");
		/* leaf */
		if (expr->elements[i] >= 0) {
			node->curve = shape_duplicate(flp_desc->units[expr->elements[i]].shape);
			node->left = node->right = NULL;
			node->label.unit = expr->elements[i];
		/*	internal node denoting a cut	*/
		} else {
			right = tree_node_stack_pop(stack);
			left = tree_node_stack_pop(stack);
			node->curve = shape_add(left->curve, right->curve, expr->elements[i]);
			node->left = left;
			node->right = right;
			node->label.cut_type = expr->elements[i];
		}
		tree_node_stack_push(stack, node);
	}

	tree_node_stack_clear(stack);
	return node;
}

void free_tree(tree_node_t *root)
{
	if (root->left != NULL)
		free_tree(root->left);
	if (root->right != NULL)
		free_tree(root->right);
	free_shape(root->curve);
	free(root);
}

/* debug print	*/
void print_tree(tree_node_t *root, flp_desc_t *flp_desc)
{
	if(root->left != NULL)
		print_tree(root->left, flp_desc);
	if (root->right != NULL)
		print_tree(root->right, flp_desc);
		
	if (root->label.unit >= 0)
		fprintf(stdout, "printing shape curve for %s\n", flp_desc->units[root->label.unit].name);
	else if (root->label.cut_type == CUT_VERTICAL)	
		fprintf(stdout, "printing shape curve for VERTICAL CUT\n");
	else if (root->label.cut_type == CUT_HORIZONTAL)
		fprintf(stdout, "printing shape curve for HORIZONTAL CUT\n");
	else
		fprintf(stdout, "unknown cut type\n");
	
	print_shape(root->curve);
}

/* 
 * print only the portion of the shape curves 
 * corresponding to the `pos'th entry of root->curve
 */
void print_tree_relevant(tree_node_t *root, int pos, flp_desc_t *flp_desc)
{
	if(root->left != NULL)
		print_tree_relevant(root->left, root->curve->left_pos[pos], flp_desc);
	if (root->right != NULL)
		print_tree_relevant(root->right, root->curve->right_pos[pos], flp_desc);
		
	if (root->label.unit >= 0)
		fprintf(stdout, "printing orientation for %s:\t", flp_desc->units[root->label.unit].name);
	else if (root->label.cut_type == CUT_VERTICAL)	
		fprintf(stdout, "printing orientation for VERTICAL CUT:\t");
	else if (root->label.cut_type == CUT_HORIZONTAL)
		fprintf(stdout, "printing orientation for HORIZONTAL CUT:\t");
	else
		fprintf(stdout, "unknown cut type\n");
	
	print_shape_entry(root->curve, pos);
}

int min_area_pos(shape_t *curve)
{
	int i = 1, pos = 0;
	double min = curve->x[0] * curve->y[0];
	for (; i < curve->size; i++)
		if (min > curve->x[i] * curve->y[i]) {
			min = curve->x[i] * curve->y[i];
			pos = i;
		}
	return pos;	
}

/* 
 * recursive sizing - given the slicing tree containing
 * the added up shape curves. 'pos' denotes the current
 * added up orientation. 'leftx' & 'bottomy' denote the
 * left and bottom ends of the current bounding rectangle
 */
int recursive_sizing (tree_node_t *node, int pos, 
					   double leftx, double bottomy,
					   int dead_count, int compact_dead,
					   double compact_ratio,
#if VERBOSE > 1					   
					   double *compacted_area, 
#endif					   
					   flp_t *flp)
{
	/* shortcut	*/
	shape_t *self = node->curve;

	/* leaf node. fill the placeholder	*/
	if (node->label.unit >= 0) {
		flp->units[node->label.unit].width = self->x[pos];
		flp->units[node->label.unit].height = self->y[pos];
		flp->units[node->label.unit].leftx = leftx;
		flp->units[node->label.unit].bottomy = bottomy;
	} else {
		/* shortcuts */
		int idx;
		double x1, x2, y1, y2;
		shape_t *left = node->left->curve;
		shape_t *right = node->right->curve;

		/* location of the first dead block	+ offset */
		idx = (flp->n_units + 1) / 2 + dead_count;

		x1 = left->x[self->left_pos[pos]];
		x2 = right->x[self->right_pos[pos]];
		y1 = left->y[self->left_pos[pos]];
		y2 = right->y[self->right_pos[pos]];

		/* add a dead block - possibly of zero area	*/
		if (node->label.unit == CUT_VERTICAL) {
			/* 
			 * if a dead block has been previously compacted away from this
			 * bounding rectangle, absorb that area into the child also
			 */
			if(self->y[pos] > MAX(y1, y2)) {
				left->y[self->left_pos[pos]] += (self->y[pos] - MAX(y1, y2));
				right->y[self->right_pos[pos]] += (self->y[pos] - MAX(y1, y2));
				y1 = left->y[self->left_pos[pos]];
				y2 = right->y[self->right_pos[pos]];
			}	
			if(self->x[pos] > (x1+x2)) {
				right->x[self->right_pos[pos]] += self->x[pos]-(x1+x2);
				x2 = right->x[self->right_pos[pos]];
			}	

			flp->units[idx].width = (y2 >= y1) ? x1 : x2;
			flp->units[idx].height = fabs(y2 - y1);
			flp->units[idx].leftx = leftx + ((y2 >= y1) ? 0 : x1);
			flp->units[idx].bottomy = bottomy + MIN(y1, y2);

			/* 
			 * ignore dead blocks smaller than compact_ratio times the area
			 * of the smaller of the constituent rectangles. instead, increase
			 * the size of the rectangle by that amount
			 */
			if (compact_dead && fabs(y2-y1) / MIN(y1, y2) <= compact_ratio) {
				#if VERBOSE > 1
				*compacted_area += (flp->units[idx].width * flp->units[idx].height);
				#endif
				if (y2 >= y1) 
					left->y[self->left_pos[pos]] = y2;
				else
					right->y[self->right_pos[pos]] = y1;
			} else {
				dead_count++;
			}

			/* left and bottom don't change for the left child	*/
			dead_count = recursive_sizing(node->left, self->left_pos[pos],
										 leftx, bottomy, dead_count, compact_dead, 
										 compact_ratio,
			#if VERBOSE > 1
										 compacted_area,
			#endif
										 flp);
			dead_count = recursive_sizing(node->right, self->right_pos[pos],
										 leftx + self->median[pos], bottomy, 
										 dead_count, compact_dead, compact_ratio,
			#if VERBOSE > 1
										 compacted_area,
			#endif
										 flp);
		} else {
			if(self->x[pos] > MAX(x1, x2)) {
				left->x[self->left_pos[pos]] += (self->x[pos] - MAX(x1, x2));
				right->x[self->right_pos[pos]] += (self->x[pos] - MAX(x1, x2));
				x1 = left->x[self->left_pos[pos]];
				x2 = right->x[self->right_pos[pos]];
			}	
			if(self->y[pos] > (y1+y2)) {
				right->y[self->right_pos[pos]] += self->y[pos]-(y1+y2);
				y2 = right->y[self->right_pos[pos]];
			}	

			flp->units[idx].width = fabs(x2 - x1);
			flp->units[idx].height = (x2 >= x1) ? y1 : y2;
			flp->units[idx].leftx = leftx + MIN(x1, x2);
			flp->units[idx].bottomy = bottomy + ((x2 >= x1) ? 0 : y1);

			if (compact_dead && fabs(x2-x1) / MIN(x1, x2) <= compact_ratio) {
				#if VERBOSE > 1
				*compacted_area += (flp->units[idx].width * flp->units[idx].height);
				#endif
				if (x2 >= x1) 
					left->x[self->left_pos[pos]] = x2;
				else
					right->x[self->right_pos[pos]] = x1;
			} else {
				dead_count++;
			}

			/* left and bottom don't change for the left child	*/
			dead_count = recursive_sizing(node->left, self->left_pos[pos],
										 leftx, bottomy, dead_count, compact_dead, 
										 compact_ratio,
			#if VERBOSE > 1
										 compacted_area,
			#endif
										 flp);
			dead_count = recursive_sizing(node->right, self->right_pos[pos],
							 			 leftx, bottomy + self->median[pos], 
										 dead_count, compact_dead, compact_ratio,
			#if VERBOSE > 1
										 compacted_area,
			#endif
										 flp);
		}
	}
	return dead_count;
}

/* 
 * convert slicing tree into actual floorplan
 * returns the number of dead blocks compacted
 */
int tree_to_flp(tree_node_t *root, flp_t *flp, int compact_dead, 
				double compact_ratio)
{
	/* for now, only choose the floorplan with
	 * the minimum area regardless of the overall
	 * aspect ratio
	 */
	int pos = min_area_pos(root->curve);
	#if VERBOSE > 1									  
	double compacted_area = 0.0;
	#endif								  
	int dead_count = recursive_sizing(root, pos, 0.0, 0.0, 0, 
					 				  compact_dead, compact_ratio,
	#if VERBOSE > 1									  
									  &compacted_area,
	#endif								  
									  flp);
					 
	int compacted = (flp->n_units - 1) / 2 - dead_count;
	flp->n_units -= compacted;
	#if VERBOSE > 1
	fprintf(stdout, "%d dead blocks, %.2f%% of the core compacted\n", compacted, 
			compacted_area / (get_total_area(flp)-compacted_area) * 100);
	#endif
	return compacted;
}

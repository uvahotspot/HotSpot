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

#include "npe.h"
#include "flp.h"
#include "util.h"

void fill_unit_pos(NPE_t *expr)
{
	int i, j=0; 
	for (i=0; i < expr->size; i++)
		if (expr->elements[i] >= 0) {
			expr->unit_pos[j] = i;
			j++;
		}
	expr->n_units = j;		
}

void fill_flip_pos(NPE_t *expr)
{
	int i, j=0; 
	for (i=0; i < expr->size - 1; i++)
		if ((expr->elements[i] < 0 && expr->elements[i+1] >= 0) ||
			(expr->elements[i] >= 0 && expr->elements[i+1] < 0)) {
			expr->flip_pos[j] = i;
			j++;
		}
	expr->n_flips = j;		
}

void fill_chain_pos(NPE_t *expr)
{
	int i=0, j=0, prev; 

	while (i < expr->size) {
		if (expr->elements[i] < 0) {
			expr->chain_pos[j] = i;
			j++;

			/* skip this chain of successive cuts	*/
			prev = expr->elements[i];
			i++;
			while(i < expr->size && expr->elements[i] < 0) {
				if (expr->elements[i] == prev)
					fatal("NPE not normalized\n");
				prev = expr->elements[i];	
				i++;
			}
		} else
			i++;
	}
	expr->n_chains = j;		
}

void fill_ballot_count(NPE_t *expr)
{
	int i, ballot_count = 0;

	for (i=0; i < expr->size; i++) {
		if (expr->elements[i] < 0)
			ballot_count++;
		expr->ballot_count[i] = ballot_count;	
	}
}

/* the starting solution for simulated annealing	*/
NPE_t *NPE_get_initial(flp_desc_t *flp_desc)
{
	int i;
	NPE_t *expr = (NPE_t *) calloc(1, sizeof(NPE_t));
	if (!expr)
		fatal("memory allocation error\n");
	expr->size = 2 * flp_desc->n_units - 1;
	expr->elements = (int *) calloc(expr->size, sizeof(int));
	expr->unit_pos = (int *) calloc(flp_desc->n_units, sizeof(int));
	expr->flip_pos = (int *) calloc(expr->size, sizeof(int));
	expr->chain_pos = (int *) calloc(flp_desc->n_units-1, sizeof(int));
	expr->ballot_count = (int *) calloc(expr->size, sizeof(int));

	if(!expr->elements || !expr->unit_pos || !expr->flip_pos 
	   || !expr->chain_pos || !expr->ballot_count)
	   fatal("memory allocation error\n");

	/* starting solution - 0, 1, V, 2, V, ..., n-1, V	*/
	expr->elements[0] = 0;
	for (i=1; i < expr->size; i+=2) {
		expr->elements[i] = (i+1) / 2;
		expr->elements[i+1] = CUT_VERTICAL;
	}

	fill_unit_pos(expr);
	fill_flip_pos(expr);
	fill_chain_pos(expr);
	fill_ballot_count(expr);

	return expr;
}

void free_NPE(NPE_t *expr)
{
	free(expr->elements);
	free(expr->unit_pos);
	free(expr->flip_pos);
	free(expr->chain_pos);
	free(expr->ballot_count);
	free(expr);
}

/* debug print	*/
void print_NPE(NPE_t *expr, flp_desc_t *flp_desc)
{
	int i;

	fprintf(stdout, "printing normalized polish expression ");
	fprintf(stdout, "of size %d\n", expr->size);
	fprintf(stdout, "%s", flp_desc->units[expr->elements[0]].name);

	for(i=1; i < expr->size; i++) {
		if (expr->elements[i] >= 0)
			fprintf(stdout, ", %s", flp_desc->units[expr->elements[i]].name);
		else if (expr->elements[i] == CUT_VERTICAL)
			fprintf(stdout, ", V");
		else if (expr->elements[i] == CUT_HORIZONTAL)
			fprintf(stdout, ", H");
		else
			fprintf(stdout, ", X");
	}
	fprintf(stdout, "\n");

	fprintf(stdout, "unit_pos:\n");
	for(i=0; i < expr->n_units; i++)
		fprintf(stdout, "%d\t", expr->unit_pos[i]);
	fprintf(stdout, "\nflip_pos:\n");
	for(i=0; i < expr->n_flips; i++)
		fprintf(stdout, "%d\t", expr->flip_pos[i]);
	fprintf(stdout, "\nchain_pos:\n");
	for(i=0; i < expr->n_chains; i++)
		fprintf(stdout, "%d\t", expr->chain_pos[i]);
	fprintf(stdout, "\nballot_count:\n");
	for(i=0; i < expr->size; i++)
		fprintf(stdout, "%d\t", expr->ballot_count[i]);
	fprintf(stdout, "\n");
}

/* 
 * move M1 of the floorplan paper 
 * swap two units adjacent in the NPE	
 */
void NPE_swap_units(NPE_t *expr, int pos)
{
	int i, t;
	
	/* find adjacent unit	*/
	for (i=pos+1; i < expr->size; i++)
		if (expr->elements[i] >= 0)
			break;
	if (i >= expr->size)
		fatal("unable to find adjacent unit\n");

	/* swap	*/
	t = expr->elements[pos];
	expr->elements[pos] = expr->elements[i]; 
	expr->elements[i] = t;
}

/* move M2 - invert a chain of cut_types in the NPE	*/
void NPE_invert_chain(NPE_t *expr, int pos)
{
	int i = pos+1, prev = expr->elements[pos];

	if (expr->elements[pos] == CUT_VERTICAL)
		expr->elements[pos] = CUT_HORIZONTAL;
	else if (expr->elements[pos] == CUT_HORIZONTAL)
		expr->elements[pos] = CUT_VERTICAL;
	else
		fatal("invalid NPE in invert_chain\n");

	while(i < expr->size && expr->elements[i] < 0) {
		if (expr->elements[i] == prev)
			fatal("NPE not normalized\n");
		prev = expr->elements[i];
		if (expr->elements[i] == CUT_VERTICAL)
			expr->elements[i] = CUT_HORIZONTAL;
		else if (expr->elements[i] == CUT_HORIZONTAL)
			expr->elements[i] = CUT_VERTICAL;
		else
			fatal("unknown cut type\n");
		i++;	
	}
}

/* binary search and increment the unit position by delta	*/
int update_unit_pos(NPE_t *expr, int pos, int delta,
					int start, int end)
{
	int mid;

	if (start > end)
		return FALSE;

	mid = (start + end) / 2;
	
	if (expr->unit_pos[mid] == pos) {
		expr->unit_pos[mid] += delta;
		return TRUE;
	} else if (expr->unit_pos[mid] > pos)
		return update_unit_pos(expr, pos, delta, start, mid-1);
	else
		return update_unit_pos(expr, pos, delta, mid+1, end);
}

/* 
 * move M3 - swap adjacent cut_type and unit in the NPE	
 * - could result in a non-allowable move. hence returns
 * if the move is legal or not
 */
int NPE_swap_cut_unit(NPE_t *expr, int pos)
{
	int t;

	if (pos <= 0 || pos >= expr->size -1)
		fatal("invalid position in NPE_swap_cut_unit\n");

	/* unit, cut_type swap	*/
	if (expr->elements[pos] >= 0) {
		/* swap leads to consecutive cut_types that are identical?	*/
		if (expr->elements[pos-1] ==  expr->elements[pos+1])
			return FALSE;
		/* move should not violate the balloting property	*/
		if (2 * expr->ballot_count[pos+1] >= pos+1)
			return FALSE;
		/* unit's position is advanced by 1	*/
		if (!update_unit_pos(expr, pos, 1, 0, expr->n_units-1))
			fatal("unit position not found\n");
		expr->ballot_count[pos]++;
	} else {	/* cut_type, unit swap	*/
		/* swap leads to consecutive cut_types that are identical?	*/
		if ((pos < expr->size - 2) && (expr->elements[pos] ==  expr->elements[pos+2]))
			return FALSE;
		/* unit's position is reduced by 1	*/
		if (!update_unit_pos(expr, pos+1, -1, 0, expr->n_units-1))
			fatal("unit position not found\n");
		expr->ballot_count[pos]--;
	}
	
	/* swap O.K	*/
	t = expr->elements[pos];
	expr->elements[pos] = expr->elements[pos+1];
	expr->elements[pos+1] = t;

	/* flip and chain positions altered. recompute them	*/
	fill_flip_pos(expr);
	fill_chain_pos(expr);

	return TRUE;
}

/* make a random move out of the above	*/
NPE_t *make_random_move(NPE_t *expr)
{
	int i, move, count = 0, done = FALSE, m3_count;
	NPE_t *copy = NPE_duplicate(expr);

	while (!done && count < MAX_MOVES) {
		/* choose one of three moves	*/
		move = rand_upto(3);
		switch(move) {
			case 0:	/* swap adjacent units	*/
				/* leave the unit last in the NPE	*/
				i = rand_upto(expr->n_units-1);
				#if VERBOSE > 2
				fprintf(stdout, "making M1 at %d\n", expr->unit_pos[i]);
				#endif
				NPE_swap_units(copy, expr->unit_pos[i]);
				done = TRUE;
				break;

			case 1:	/* invert an arbitrary chain	*/
				i = rand_upto(expr->n_chains);
				#if VERBOSE > 2
				fprintf(stdout, "making M2 at %d\n", expr->chain_pos[i]);
				#endif
				NPE_invert_chain(copy, expr->chain_pos[i]);
				done = TRUE;
				break;

			case 2:	/* swap a unit and an adjacent cut_type	*/
				m3_count = 0; 
				while (!done && m3_count < MAX_MOVES) {
					i = rand_upto(expr->n_flips);
					#if VERBOSE > 2
					fprintf(stdout, "making M3 at %d\n", expr->flip_pos[i]);
					#endif
					done = NPE_swap_cut_unit(copy, expr->flip_pos[i]);
					m3_count++;
				}
				break;

			default:
				fatal("unknown move type\n");
				break;
		}
		count++;
	}

	if (count == MAX_MOVES) {
		char msg[STR_SIZE];
		sprintf(msg, "tried %d moves, now giving up\n", MAX_MOVES); 
		fatal(msg);
	}

	return copy;
}

/* make a copy of this NPE	*/
NPE_t *NPE_duplicate(NPE_t *expr)
{
	int i;
	NPE_t *copy = (NPE_t *) calloc(1, sizeof(NPE_t));
	if (!copy)
		fatal("memory allocation error\n");
	copy->elements = (int *) calloc(expr->size, sizeof(int));
	copy->unit_pos = (int *) calloc(expr->n_units, sizeof(int));
	copy->flip_pos = (int *) calloc(expr->size, sizeof(int));
	copy->chain_pos = (int *) calloc(expr->n_units-1, sizeof(int));
	copy->ballot_count = (int *) calloc(expr->size, sizeof(int));

	if(!copy->elements || !copy->unit_pos || !copy->flip_pos 
	   || !copy->chain_pos || !copy->ballot_count)
	   fatal("memory allocation error\n");

	copy->size = expr->size;
	for (i=0; i < expr->size; i++)
		copy->elements[i] = expr->elements[i];
	copy->n_units = expr->n_units;
	for (i=0; i < expr->n_units; i++)
		copy->unit_pos[i] = expr->unit_pos[i];
	copy->n_flips = expr->n_flips;
	for (i=0; i < expr->n_flips; i++)
		copy->flip_pos[i] = expr->flip_pos[i];
	copy->n_chains = expr->n_chains;
	for (i=0; i < expr->n_chains; i++)
		copy->chain_pos[i] = expr->chain_pos[i];
	for (i=0; i < expr->size; i++)
		copy->ballot_count[i] = expr->ballot_count[i];
	
	return copy;
}


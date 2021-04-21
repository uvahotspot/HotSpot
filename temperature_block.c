#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef _MSC_VER
#define strcasecmp    _stricmp
#define strncasecmp   _strnicmp
#else
#include <strings.h>
#endif

#include "temperature_block.h"
#include "flp.h"
#include "util.h"

/*
 * allocate memory for the matrices. placeholder can be an empty
 * floorplan frame with only the names of the functional units
 */
block_model_t *alloc_block_model(thermal_config_t *config, flp_t *placeholder)
{
	/* shortcuts	*/
	int n = placeholder->n_units;
	int m = NL*n+EXTRA;

	block_model_t *model = (block_model_t *) calloc (1, sizeof(block_model_t));
	if (!model)
		fatal("memory allocation error\n");
	model->config = *config;
	model->n_units = model->base_n_units = n;
	model->n_nodes = m;

	model->border = imatrix(n, 4);
	model->len = dmatrix(n, n);	/* len[i][j] = length of shared edge bet. i & j	*/
	model->g = dmatrix(m, m);	/* g[i][j] = conductance bet. nodes i & j */
	model->gx = dvector(n);		/* lumped conductances in x direction	*/
	model->gy = dvector(n);		/* lumped conductances in y direction	*/
	model->gx_int = dvector(n);	/* lateral conductances in the interface layer	*/
	model->gy_int = dvector(n);
	model->gx_sp = dvector(n);	/* lateral conductances in the spreader	layer */
	model->gy_sp = dvector(n);
	model->gx_hs = dvector(n);	/* lateral conductances in the heatsink	layer */
	model->gy_hs = dvector(n);
	/* vertical conductances to ambient	*/
	model->g_amb = dvector(n+EXTRA);
	model->t_vector = dvector(m);/* scratch pad	*/
	model->p = ivector(m);		/* permutation vector for b's LUP decomposition	*/

	model->a = dvector(m);		/* vertical Cs - diagonal matrix stored as a 1-d vector	*/
	model->inva = dvector(m);	/* inverse of the above 	*/
	/* B, C and LU are (NL*n+EXTRA)x(NL*n+EXTRA) matrices	*/
	model->b = dmatrix(m, m);
	model->c = dmatrix(m, m);
	model->lu = dmatrix(m, m);

	model->flp = placeholder;
	return model;
}

/* creates matrices  B and invB: BT = Power in the steady state.
 * NOTE: EXTRA nodes: 4 heat spreader peripheral nodes, 4 heat
 * sink inner peripheral nodes, 4 heat sink outer peripheral
 * nodes(north, south, east and west) and 1 ambient node.
 */
void populate_R_model_block(block_model_t *model, flp_t *flp)
{
	/*	shortcuts	*/
	double **b = model->b;
	double *gx = model->gx, *gy = model->gy;
	double *gx_int = model->gx_int, *gy_int = model->gy_int;
	double *gx_sp = model->gx_sp, *gy_sp = model->gy_sp;
	double *gx_hs = model->gx_hs, *gy_hs = model->gy_hs;
	double *g_amb = model->g_amb;
	double **len = model->len, **g = model->g, **lu = model->lu;
	int **border = model->border;
	int *p = model->p;
	double t_chip = model->config.t_chip;
	double r_convec = model->config.r_convec;
	double s_sink = model->config.s_sink;
	double t_sink = model->config.t_sink;
	double s_spreader = model->config.s_spreader;
	double t_spreader = model->config.t_spreader;
	double t_interface = model->config.t_interface;
	double k_chip = model->config.k_chip;
	double k_sink = model->config.k_sink;
	double k_spreader = model->config.k_spreader;
	double k_interface = model->config.k_interface;

	int i, j, n = flp->n_units;
	double gn_sp=0, gs_sp=0, ge_sp=0, gw_sp=0;
	double gn_hs=0, gs_hs=0, ge_hs=0, gw_hs=0;
	double r_amb;

	double w_chip = get_total_width (flp);	/* x-axis	*/
	double l_chip = get_total_height (flp);	/* y-axis	*/

	/* sanity check on floorplan sizes	*/
	if (w_chip > s_sink || l_chip > s_sink ||
		w_chip > s_spreader || l_chip > s_spreader) {
		print_flp(flp, FALSE);
		print_flp_fig(flp);
		fatal("inordinate floorplan size!\n");
	}
	if(model->flp != flp || model->n_units != flp->n_units ||
	   model->n_nodes != NL * flp->n_units + EXTRA)
	   fatal("mismatch between the floorplan and the thermal model\n");

	/* gx's and gy's of blocks	*/
	for (i = 0; i < n; i++) {
		/* at the silicon layer	*/
		if (model->config.block_omit_lateral) {
			gx[i] = gy[i] = 0;
		}
		else {
			gx[i] = 1.0/getr(k_chip, flp->units[i].width / 2.0, flp->units[i].height * t_chip);
			gy[i] = 1.0/getr(k_chip, flp->units[i].height / 2.0, flp->units[i].width * t_chip);
		}

		/* at the interface layer	*/
		gx_int[i] = 1.0/getr(k_interface, flp->units[i].width / 2.0, flp->units[i].height * t_interface);
		gy_int[i] = 1.0/getr(k_interface, flp->units[i].height / 2.0, flp->units[i].width * t_interface);

		/* at the spreader layer	*/
		gx_sp[i] = 1.0/getr(k_spreader, flp->units[i].width / 2.0, flp->units[i].height * t_spreader);
		gy_sp[i] = 1.0/getr(k_spreader, flp->units[i].height / 2.0, flp->units[i].width * t_spreader);

		/* at the heatsink layer	*/
		gx_hs[i] = 1.0/getr(k_sink, flp->units[i].width / 2.0, flp->units[i].height * t_sink);
		gy_hs[i] = 1.0/getr(k_sink, flp->units[i].height / 2.0, flp->units[i].width * t_sink);
	}

	/* shared lengths between blocks	*/
	for (i = 0; i < n; i++)
		for (j = i; j < n; j++)
			len[i][j] = len[j][i] = get_shared_len(flp, i, j);

	/* package R's	*/
	populate_package_R(&model->pack, &model->config, w_chip, l_chip);

	/* short the R's from block centers to a particular chip edge	*/
	for (i = 0; i < n; i++) {
		if (eq(flp->units[i].bottomy + flp->units[i].height, l_chip)) {
			gn_sp += gy_sp[i];
			gn_hs += gy_hs[i];
			border[i][2] = 1;	/* block is on northern border 	*/
		} else
			border[i][2] = 0;

		if (eq(flp->units[i].bottomy, 0)) {
			gs_sp += gy_sp[i];
			gs_hs += gy_hs[i];
			border[i][3] = 1;	/* block is on southern border	*/
		} else
			border[i][3] = 0;

		if (eq(flp->units[i].leftx + flp->units[i].width, w_chip)) {
			ge_sp += gx_sp[i];
			ge_hs += gx_hs[i];
			border[i][1] = 1;	/* block is on eastern border	*/
		} else
			border[i][1] = 0;

		if (eq(flp->units[i].leftx, 0)) {
			gw_sp += gx_sp[i];
			gw_hs += gx_hs[i];
			border[i][0] = 1;	/* block is on western border	*/
		} else
			border[i][0] = 0;
	}

	/* initialize g	*/
	zero_dmatrix(g, NL*n+EXTRA, NL*n+EXTRA);
	zero_dvector(g_amb, n+EXTRA);

	/* overall Rs between nodes */
	for (i = 0; i < n; i++) {
		double area = (flp->units[i].height * flp->units[i].width);
		/* amongst functional units	in the various layers	*/
		for (j = 0; j < n; j++) {
			double part = 0, part_int = 0, part_sp = 0, part_hs = 0;
			if (is_horiz_adj(flp, i, j)) {
				part = gx[i] / flp->units[i].height;
				part_int = gx_int[i] / flp->units[i].height;
				part_sp = gx_sp[i] / flp->units[i].height;
				part_hs = gx_hs[i] / flp->units[i].height;
			}
			else if (is_vert_adj(flp, i,j))  {
				part = gy[i] / flp->units[i].width;
				part_int = gy_int[i] / flp->units[i].width;
				part_sp = gy_sp[i] / flp->units[i].width;
				part_hs = gy_hs[i] / flp->units[i].width;
			}
			g[i][j] = part * len[i][j];
			g[IFACE*n+i][IFACE*n+j] = part_int * len[i][j];
			g[HSP*n+i][HSP*n+j] = part_sp * len[i][j];
			g[HSINK*n+i][HSINK*n+j] = part_hs * len[i][j];
		}
		/* the 2.0 factor in the following equations is
		 * explained during the calculation of the B matrix
		 */
 		/* vertical g's in the silicon layer	*/
		g[i][IFACE*n+i]=g[IFACE*n+i][i]=2.0/getr(k_chip, t_chip, area);
 		/* vertical g's in the interface layer	*/
		g[IFACE*n+i][HSP*n+i]=g[HSP*n+i][IFACE*n+i]=2.0/getr(k_interface, t_interface, area);
		/* vertical g's in the spreader layer	*/
		g[HSP*n+i][HSINK*n+i]=g[HSINK*n+i][HSP*n+i]=2.0/getr(k_spreader, t_spreader, area);
		/* vertical g's in the heatsink core layer	*/
		/* vertical R to ambient: divide r_convec proportional to area	*/
		r_amb = r_convec * (s_sink * s_sink) / area;
		g_amb[i] = 1.0 / (getr(k_sink, t_sink, area) + r_amb);

		/* lateral g's from block center (spreader layer) to peripheral (n,s,e,w) spreader nodes	*/
		g[HSP*n+i][NL*n+SP_N]=g[NL*n+SP_N][HSP*n+i]=2.0*border[i][2] /
							  ((1.0/gy_sp[i])+model->pack.r_sp1_y*gn_sp/gy_sp[i]);
		g[HSP*n+i][NL*n+SP_S]=g[NL*n+SP_S][HSP*n+i]=2.0*border[i][3] /
							  ((1.0/gy_sp[i])+model->pack.r_sp1_y*gs_sp/gy_sp[i]);
		g[HSP*n+i][NL*n+SP_E]=g[NL*n+SP_E][HSP*n+i]=2.0*border[i][1] /
							  ((1.0/gx_sp[i])+model->pack.r_sp1_x*ge_sp/gx_sp[i]);
		g[HSP*n+i][NL*n+SP_W]=g[NL*n+SP_W][HSP*n+i]=2.0*border[i][0] /
							  ((1.0/gx_sp[i])+model->pack.r_sp1_x*gw_sp/gx_sp[i]);

		/* lateral g's from block center (heatsink layer) to peripheral (n,s,e,w) heatsink nodes	*/
		g[HSINK*n+i][NL*n+SINK_C_N]=g[NL*n+SINK_C_N][HSINK*n+i]=2.0*border[i][2] /
									((1.0/gy_hs[i])+model->pack.r_hs1_y*gn_hs/gy_hs[i]);
		g[HSINK*n+i][NL*n+SINK_C_S]=g[NL*n+SINK_C_S][HSINK*n+i]=2.0*border[i][3] /
									((1.0/gy_hs[i])+model->pack.r_hs1_y*gs_hs/gy_hs[i]);
		g[HSINK*n+i][NL*n+SINK_C_E]=g[NL*n+SINK_C_E][HSINK*n+i]=2.0*border[i][1] /
									((1.0/gx_hs[i])+model->pack.r_hs1_x*ge_hs/gx_hs[i]);
		g[HSINK*n+i][NL*n+SINK_C_W]=g[NL*n+SINK_C_W][HSINK*n+i]=2.0*border[i][0] /
									((1.0/gx_hs[i])+model->pack.r_hs1_x*gw_hs/gx_hs[i]);
	}

	/* g's from peripheral(n,s,e,w) nodes	*/
	/* vertical g's between peripheral spreader nodes and center peripheral heatsink nodes */
	g[NL*n+SP_N][NL*n+SINK_C_N]=g[NL*n+SINK_C_N][NL*n+SP_N]=2.0/model->pack.r_sp_per_y;
	g[NL*n+SP_S][NL*n+SINK_C_S]=g[NL*n+SINK_C_S][NL*n+SP_S]=2.0/model->pack.r_sp_per_y;
	g[NL*n+SP_E][NL*n+SINK_C_E]=g[NL*n+SINK_C_E][NL*n+SP_E]=2.0/model->pack.r_sp_per_x;
	g[NL*n+SP_W][NL*n+SINK_C_W]=g[NL*n+SINK_C_W][NL*n+SP_W]=2.0/model->pack.r_sp_per_x;
	/* lateral g's between peripheral outer sink nodes and center peripheral sink nodes	*/
	g[NL*n+SINK_C_N][NL*n+SINK_N]=g[NL*n+SINK_N][NL*n+SINK_C_N]=2.0/(model->pack.r_hs + model->pack.r_hs2_y);
	g[NL*n+SINK_C_S][NL*n+SINK_S]=g[NL*n+SINK_S][NL*n+SINK_C_S]=2.0/(model->pack.r_hs + model->pack.r_hs2_y);
	g[NL*n+SINK_C_E][NL*n+SINK_E]=g[NL*n+SINK_E][NL*n+SINK_C_E]=2.0/(model->pack.r_hs + model->pack.r_hs2_x);
	g[NL*n+SINK_C_W][NL*n+SINK_W]=g[NL*n+SINK_W][NL*n+SINK_C_W]=2.0/(model->pack.r_hs + model->pack.r_hs2_x);
	/* vertical g's between inner peripheral sink nodes and ambient	*/
	g_amb[n+SINK_C_N] = g_amb[n+SINK_C_S] = 1.0 / (model->pack.r_hs_c_per_y+model->pack.r_amb_c_per_y);
	g_amb[n+SINK_C_E] = g_amb[n+SINK_C_W] = 1.0 / (model->pack.r_hs_c_per_x+model->pack.r_amb_c_per_x);
	/* vertical g's between outer peripheral sink nodes and ambient	*/
	g_amb[n+SINK_N] = g_amb[n+SINK_S] = g_amb[n+SINK_E] =
					  g_amb[n+SINK_W] = 1.0 / (model->pack.r_hs_per+model->pack.r_amb_per);

	/* calculate matrix B such that BT = POWER in steady state */
	/* non-diagonal elements	*/
	for (i = 0; i < NL*n+EXTRA; i++)
		for (j = 0; j < i; j++)
			if ((g[i][j] == 0.0) || (g[j][i] == 0.0))
				b[i][j] = b[j][i] = 0.0;
			else
				/* here is why the 2.0 factor comes when calculating g[][]	*/
				b[i][j] = b[j][i] = -1.0/((1.0/g[i][j])+(1.0/g[j][i]));
	/* diagonal elements	*/
	for (i = 0; i < NL*n+EXTRA; i++) {
		/* functional blocks in the heat sink layer	*/
		if (i >= HSINK*n && i < NL*n)
			b[i][i] = g_amb[i%n];
		/* heat sink peripheral nodes	*/
		else if (i >= NL*n+SINK_C_W)
			b[i][i] = g_amb[n+i-NL*n];
		/* all other nodes that are not connected to the ambient	*/
		else
			b[i][i] = 0.0;
		/* sum up the conductances	*/
		for(j=0; j < NL*n+EXTRA; j++)
			if (i != j)
				b[i][i] -= b[i][j];
	}

	/* compute the LUP decomposition of B and store it too	*/
	copy_dmatrix(lu, b, NL*n+EXTRA, NL*n+EXTRA);
	/*
	 * B is a symmetric positive definite matrix. It is
	 * symmetric because if a node A is connected to B,
	 * then B is also connected to A with the same R value.
	 * It is positive definite because of the following
	 * informal argument from Professor Lieven Vandenberghe's
	 * lecture slides for the spring 2004-2005 EE 103 class
	 * at UCLA: http://www.ee.ucla.edu/~vandenbe/103/chol.pdf
	 * x^T*B*x = voltage^T * (B*x) = voltage^T * current
	 * = total power dissipated in the resistors > 0
	 * for x != 0.
	 */
	lupdcmp(lu, NL*n+EXTRA, p, 1);

	/* done	*/
	model->flp = flp;
	model->r_ready = TRUE;
}

/* creates 2 matrices: invA, C: dT + A^-1*BT = A^-1*Power,
 * C = A^-1 * B. note that A is a diagonal matrix (no lateral
 * capacitances. all capacitances are to ground). also note that
 * it is stored as a 1-d vector. so, for computing the inverse,
 * inva[i] = 1/a[i] is just enough.
 */
void populate_C_model_block(block_model_t *model, flp_t *flp)
{
	/*	shortcuts	*/
	double *inva = model->inva, **c = model->c;
	double **b = model->b;
	double *a = model->a;
	double t_chip = model->config.t_chip;
	double c_convec = model->config.c_convec;
	double s_sink = model->config.s_sink;
	double t_sink = model->config.t_sink;
	double t_spreader = model->config.t_spreader;
	double t_interface = model->config.t_interface;
	double p_chip = model->config.p_chip;
	double p_sink = model->config.p_sink;
	double p_spreader = model->config.p_spreader;
	double p_interface = model->config.p_interface;
	double c_amb;
	double w_chip, l_chip;

	int i, n = flp->n_units;

	if (!model->r_ready)
		fatal("R model not ready\n");
	if (model->flp != flp || model->n_units != flp->n_units ||
		model->n_nodes != NL * flp->n_units + EXTRA)
		fatal("different floorplans for R and C models!\n");

	w_chip = get_total_width (flp);	/* x-axis	*/
	l_chip = get_total_height (flp);	/* y-axis	*/

	/* package C's	*/
	populate_package_C(&model->pack, &model->config, w_chip, l_chip);

	/* functional block C's */
	for (i = 0; i < n; i++) {
		double area = (flp->units[i].height * flp->units[i].width);
		/* C's from functional units to ground	*/
		a[i] = getcap(p_chip, t_chip, area);
		/* C's from interface portion of the functional units to ground	*/
		a[IFACE*n+i] = getcap(p_interface, t_interface, area);
		/* C's from spreader portion of the functional units to ground	*/
		a[HSP*n+i] = getcap(p_spreader, t_spreader, area);
		/* C's from heatsink portion of the functional units to ground	*/
		/* vertical C to ambient: divide c_convec proportional to area	*/
		c_amb = C_FACTOR * c_convec / (s_sink * s_sink) * area;
		a[HSINK*n+i] = getcap(p_sink, t_sink, area) + c_amb;
	}

	/* C's from peripheral(n,s,e,w) nodes	*/
 	/* from peripheral spreader nodes to ground	*/
	a[NL*n+SP_N] = a[NL*n+SP_S] = model->pack.c_sp_per_y;
	a[NL*n+SP_E] = a[NL*n+SP_W] = model->pack.c_sp_per_x;
 	/* from center peripheral sink nodes to ground
	 * NOTE: this treatment of capacitances (and
	 * the corresponding treatment of resistances
	 * in populate_R_model) as parallel (series)
	 * is only approximate and is done in order
	 * to avoid creating an extra layer of nodes
	 */
	a[NL*n+SINK_C_N] = a[NL*n+SINK_C_S] = model->pack.c_hs_c_per_y +
										  model->pack.c_amb_c_per_y;
	a[NL*n+SINK_C_E] = a[NL*n+SINK_C_W] = model->pack.c_hs_c_per_x +
										  model->pack.c_amb_c_per_x;
	/* from outer peripheral sink nodes to ground	*/
	a[NL*n+SINK_N] = a[NL*n+SINK_S] = a[NL*n+SINK_E] = a[NL*n+SINK_W] =
					 model->pack.c_hs_per + model->pack.c_amb_per;

	/* calculate A^-1 (for diagonal matrix A) such that A(dT) + BT = POWER */
	for (i = 0; i < NL*n+EXTRA; i++)
		inva[i] = 1.0/a[i];

	/* we are always going to use the eqn dT + A^-1 * B T = A^-1 * POWER. so, store  C = A^-1 * B	*/
	diagmatmult(c, inva, b, NL*n+EXTRA);

	/*	done	*/
	model->c_ready = TRUE;
}

/* setting package nodes' power numbers	*/
void set_internal_power_block(block_model_t *model, double *power)
{
	int i;
	zero_dvector(&power[IFACE*model->n_units], model->n_units);
	zero_dvector(&power[HSP*model->n_units], model->n_units);
	for(i=0; i < model->n_units + EXTRA; i++)
		power[HSINK*model->n_units+i] = model->config.ambient * model->g_amb[i];
}

/* power and temp should both be alloced using hotspot_vector.
 * 'b' is the 'thermal conductance' matrix. i.e, b * temp = power
 *  => temp = invb * power. instead of computing invb, we have
 * stored the LUP decomposition of B in 'lu' and 'p'. Using
 * forward and backward substitution, we can then solve the
 * equation b * temp = power.
 */
void steady_state_temp_block(block_model_t *model, double *power, double *temp)
{
	if (!model->r_ready)
		fatal("R model not ready\n");

	/* set power numbers for the virtual nodes */
	set_internal_power_block(model, power);

	/*
	 * find temperatures (spd flag is set to 1 by the same argument
	 * as mentioned in the populate_R_model_block function)
	 */
	lusolve(model->lu, model->n_nodes, model->p, power, temp, 1);
}

/* compute the slope vector dy for the transient equation
 * dy + cy = p. useful in the transient solver
 */
void slope_fn_block(block_model_t *model, double *y, double *p, double *dy)
{
	/* shortcuts	*/
	int n = model->n_nodes;
	double **c = model->c;

	/* for our equation, dy = p - cy */
	#if (MATHACCEL == MA_INTEL || MATHACCEL == MA_APPLE)
	/* dy = p	*/
	cblas_dcopy(n, p, 1, dy, 1);
	/* dy = dy - c*y = p - c*y */
	cblas_dgemv(CblasRowMajor, CblasNoTrans, n, n, -1, c[0],
				n, y, 1, 1, dy, 1);
	#elif (MATHACCEL == MA_AMD || MATHACCEL == MA_SUN)
	/* dy = p	*/
	dcopy(n, p, 1, dy, 1);
	/* dy = dy - c*y = p - c*y */
	dgemv('T', n, n, -1, c[0], n, y, 1, 1, dy, 1);
	#else
	int i;
	double *t = dvector(n);
	matvectmult(t, c, y, n);
	for (i = 0; i < n; i++)
		dy[i] = p[i]-t[i];
	free_dvector(t);
	#endif
}

/* compute_temp: solve for temperature from the equation dT + CT = inv_A * Power
 * Given the temperature (temp) at time t, the power dissipation per cycle during the
 * last interval (time_elapsed), find the new temperature at time t+time_elapsed.
 * power and temp should both be alloced using hotspot_vector
 */
void compute_temp_block(block_model_t *model, double *power, double *temp, double time_elapsed)
{
	double t, h, new_h;

	#if VERBOSE > 1
	unsigned int i = 0;
	#endif

	if (!model->r_ready || !model->c_ready)
		fatal("block model not ready\n");
	if (temp == model->t_vector)
		fatal("output same as scratch pad\n");

	/* set power numbers for the virtual nodes */
	set_internal_power_block(model, power);

	/* use the scratch pad vector to find (inv_A)*POWER */
	diagmatvectmult(model->t_vector, model->inva, power, model->n_nodes);

	/* Obtain temp at time (t+time_elapsed).
	 * Instead of getting the temperature at t+time_elapsed directly, we do it
	 * in multiple steps with the correct step size at each time
	 * provided by rk4
	 */
	for (t = 0, new_h = MIN_STEP; t < time_elapsed && new_h >= MIN_STEP*DELTA; t+=h) {
		h = new_h;
		new_h = rk4(model, temp, model->t_vector, model->n_nodes, &h,
		/* the slope function callback is typecast accordingly */
					temp, (slope_fn_ptr) slope_fn_block);
		new_h = MIN(new_h, time_elapsed-t-h);
		#if VERBOSE > 1
		i++;
		#endif
	}
	#if VERBOSE > 1
	fprintf(stdout, "no. of rk4 calls during compute_temp: %d\n", i+1);
	#endif
}

/* differs from 'dvector()' in that memory for internal nodes is also allocated	*/
double *hotspot_vector_block(block_model_t *model)
{
	return dvector(model->n_nodes);
}

/* copy 'src' to 'dst' except for a window of 'size'
 * elements starting at 'at'. useful in floorplan
 * compaction
 */
void trim_hotspot_vector_block(block_model_t *model, double *dst, double *src,
						 	   int at, int size)
{
	int i;

	for (i=0; i < at && i < model->n_nodes; i++)
		dst[i] = src[i];
	for(i=at+size; i < model->n_nodes; i++)
		dst[i-size] = src[i];
}

/* update the model's node count	*/
void resize_thermal_model_block(block_model_t *model, int n_units)
{
	if (n_units > model->base_n_units)
		fatal("resizing block model to more than the allocated space\n");
	model->n_units = n_units;
	model->n_nodes = NL * n_units + EXTRA;
	/* resize the 2-d matrices whose no. of columns changes	*/
	resize_dmatrix(model->len, model->n_units, model->n_units);
	resize_dmatrix(model->g, model->n_nodes, model->n_nodes);
	resize_dmatrix(model->b, model->n_nodes, model->n_nodes);
	resize_dmatrix(model->c, model->n_nodes, model->n_nodes);
	resize_dmatrix(model->lu, model->n_nodes, model->n_nodes);
}

/* sets the temperature of a vector 'temp' allocated using 'hotspot_vector'	*/
void set_temp_block(block_model_t *model, double *temp, double val)
{
	int i;
	for(i=0; i < model->n_nodes; i++)
		temp[i] = val;
}

/* dump temperature vector alloced using 'hotspot_vector' to 'file' */
void dump_temp_block(block_model_t *model, double *temp, char *file)
{
	flp_t *flp = model->flp;
	int i;
	char str[STR_SIZE];
	FILE *fp;

	if (!strcasecmp(file, "stdout"))
		fp = stdout;
	else if (!strcasecmp(file, "stderr"))
		fp = stderr;
	else
		fp = fopen (file, "w");

	if (!fp) {
		sprintf (str,"error: %s could not be opened for writing\n", file);
		fatal(str);
	}
	/* on chip temperatures	*/
	for (i=0; i < flp->n_units; i++)
		fprintf(fp, "%s\t%.2f\n", flp->units[i].name, temp[i]);

	/* interface temperatures	*/
	for (i=0; i < flp->n_units; i++)
		fprintf(fp, "iface_%s\t%.2f\n", flp->units[i].name, temp[IFACE*flp->n_units+i]);

	/* spreader temperatures	*/
	for (i=0; i < flp->n_units; i++)
		fprintf(fp, "hsp_%s\t%.2f\n", flp->units[i].name, temp[HSP*flp->n_units+i]);

	/* heatsink temperatures	*/
	for (i=0; i < flp->n_units; i++)
		fprintf(fp, "hsink_%s\t%.2f\n", flp->units[i].name, temp[HSINK*flp->n_units+i]);

	/* internal node temperatures	*/
	for (i=0; i < EXTRA; i++) {
		sprintf(str, "inode_%d", i);
		fprintf(fp, "%s\t%.2f\n", str, temp[i+NL*flp->n_units]);
	}

	if(fp != stdout && fp != stderr)
		fclose(fp);
}

void copy_temp_block(block_model_t *model, double *dst, double *src)
{
	copy_dvector(dst, src, NL*model->flp->n_units+EXTRA);
}

/*
 * read temperature vector alloced using 'hotspot_vector' from 'file'
 * which was dumped using 'dump_temp'. values are clipped to thermal
 * threshold based on 'clip'
 */
void read_temp_block(block_model_t *model, double *temp, char *file, int clip)
{
	/*	shortcuts	*/
	flp_t *flp = model->flp;
	double thermal_threshold = model->config.thermal_threshold;
	double ambient = model->config.ambient;

	int i, n, idx;
	double max=0, val;
	char *ptr, str1[LINE_SIZE], str2[LINE_SIZE];
	char name[STR_SIZE], format[STR_SIZE];
	FILE *fp;

	if (!strcasecmp(file, "stdin"))
		fp = stdin;
	else
		fp = fopen (file, "r");

	if (!fp) {
		sprintf (str1,"error: %s could not be opened for reading\n", file);
		fatal(str1);
	}

	/* temperatures of the different layers	*/
	for (n=0; n < NL; n++) {
		switch(n)
		{
			case 0:
					strcpy(format,"%s%lf");
					break;
			case IFACE:
					strcpy(format,"iface_%s%lf");
					break;
			case HSP:
					strcpy(format,"hsp_%s%lf");
					break;
			case HSINK:
					strcpy(format,"hsink_%s%lf");
					break;
			default:
					fatal("unknown layer\n");
					break;
		}
		for (i=0; i < flp->n_units; i++) {
			fgets(str1, LINE_SIZE, fp);
			if (feof(fp))
				fatal("not enough lines in temperature file\n");
			strcpy(str2, str1);

			/* ignore comments and empty lines	*/
			ptr = strtok(str1, " \r\t\n");
			if (!ptr || ptr[0] == '#') {
				i--;
				continue;
			}

			if (sscanf(str2, format, name, &val) != 2)
				fatal("invalid temperature file format\n");
			idx = get_blk_index(flp, name);
			if (idx >= 0)
				temp[idx + n*flp->n_units] = val;
			else	/* since get_blk_index calls fatal, the line below cannot be reached	*/
				fatal ("unit in temperature file not found in floorplan\n");

			/* find max temp on the chip	*/
			if (n == 0 && temp[idx] > max)
				max = temp[idx];
		}
	}

	/* internal node temperatures	*/
	for (i=0; i < EXTRA; i++) {
		fgets(str1, LINE_SIZE, fp);
		if (feof(fp))
			fatal("not enough lines in temperature file\n");
		strcpy(str2, str1);
		/* ignore comments and empty lines	*/
		ptr = strtok(str1, " \r\t\n");
		if (!ptr || ptr[0] == '#') {
			i--;
			continue;
		}
		if (sscanf(str2, "%s%lf", name, &val) != 2)
			fatal("invalid temperature file format\n");
		sprintf(str1, "inode_%d", i);
		if (strcasecmp(str1, name))
			fatal("invalid temperature file format\n");
		temp[i+NL*flp->n_units] = val;
	}

	fgets(str1, LINE_SIZE, fp);
	if (!feof(fp))
		fatal("too many lines in temperature file\n");

	if(fp != stdin)
		fclose(fp);

	/* clipping	*/
	if (clip && (max > thermal_threshold)) {
		/* if max has to be brought down to thermal_threshold,
		 * (w.r.t the ambient) what is the scale down factor?
		 */
		double factor = (thermal_threshold - ambient) / (max - ambient);

		/* scale down all temperature differences (from ambient) by the same factor	*/
		for (i=0; i < NL*flp->n_units + EXTRA; i++)
			temp[i] = (temp[i]-ambient)*factor + ambient;
	}
}

/* dump power numbers to file	*/
void dump_power_block(block_model_t *model, double *power, char *file)
{
	flp_t *flp = model->flp;
	int i;
	char str[STR_SIZE];
	FILE *fp;

	if (!strcasecmp(file, "stdout"))
		fp = stdout;
	else if (!strcasecmp(file, "stderr"))
		fp = stderr;
	else
		fp = fopen (file, "w");

	if (!fp) {
		sprintf (str,"error: %s could not be opened for writing\n", file);
		fatal(str);
	}
	for (i=0; i < flp->n_units; i++)
		fprintf(fp, "%s\t%.6f\n", flp->units[i].name, power[i]);
	if(fp != stdout && fp != stderr)
		fclose(fp);
}

/*
 * read power vector alloced using 'hotspot_vector' from 'file'
 * which was dumped using 'dump_power'.
 */
void read_power_block (block_model_t *model, double *power, char *file)
{
	flp_t *flp = model->flp;
	int idx;
	double val;
	char *ptr, str1[LINE_SIZE], str2[LINE_SIZE], name[STR_SIZE];
	FILE *fp;

	if (!strcasecmp(file, "stdin"))
		fp = stdin;
	else
		fp = fopen (file, "r");

	if (!fp) {
		sprintf (str1,"error: %s could not be opened for reading\n", file);
		fatal(str1);
	}
	while(!feof(fp)) {
		fgets(str1, LINE_SIZE, fp);
		if (feof(fp))
			break;
		strcpy(str2, str1);

		/* ignore comments and empty lines	*/
		ptr = strtok(str1, " \r\t\n");
		if (!ptr || ptr[0] == '#')
			continue;

		if (sscanf(str2, "%s%lf", name, &val) != 2)
			fatal("invalid power file format\n");
		idx = get_blk_index(flp, name);
		if (idx >= 0)
			power[idx] = val;
		else	/* since get_blk_index calls fatal, the line below cannot be reached	*/
			fatal ("unit in power file not found in floorplan\n");
	}
	if(fp != stdin)
		fclose(fp);
}

double find_max_temp_block(block_model_t *model, double *temp)
{
	int i;
	double max = 0.0;
	for(i=0; i < model->n_units; i++) {
		if (temp[i] < 0)
			fatal("negative temperature!\n");
		else if (max < temp[i])
			max = temp[i];
	}

	return max;
}

double find_avg_temp_block(block_model_t *model, double *temp)
{
	int i;
	double sum = 0.0;
	for(i=0; i < model->n_units; i++) {
		if (temp[i] < 0)
			fatal("negative temperature!\n");
		else
			sum += temp[i];
	}

	return (sum / model->n_units);
}

/* calculate avg sink temp for natural convection package model */
double calc_sink_temp_block(block_model_t *model, double *temp, thermal_config_t *config)
{
	flp_t *flp = model->flp;
	int i;
	double sum = 0.0;
	double width = get_total_width(flp);
	double height = get_total_height(flp);
	double spr_size = config->s_spreader*config->s_spreader;
	double sink_size = config->s_sink*config->s_sink;

	/* heatsink temperatures	*/
	for (i=0; i < flp->n_units; i++)
		if (temp[HSINK*flp->n_units+i] < 0)
			fatal("negative temperature!\n");
		else  /* area-weighted average */
			sum += temp[HSINK*flp->n_units+i]*(flp->units[i].width*flp->units[i].height);

	for(i=SINK_C_W; i <= SINK_C_E; i++)
		if (temp[i+NL*flp->n_units] < 0)
			fatal("negative temperature!\n");
		else
			sum += temp[i+NL*flp->n_units]*0.25*(config->s_spreader+height)*(config->s_spreader-width);

	for(i=SINK_C_N; i <= SINK_C_S; i++)
		if (temp[i+NL*flp->n_units] < 0)
			fatal("negative temperature!\n");
		else
			sum += temp[i+NL*flp->n_units]*0.25*(config->s_spreader+width)*(config->s_spreader-height);

	for(i=SINK_W; i <= SINK_S; i++)
		if (temp[i+NL*flp->n_units] < 0)
			fatal("negative temperature!\n");
		else
			sum += temp[i+NL*flp->n_units]*0.25*(sink_size-spr_size);

	return (sum / sink_size);
}

void delete_block_model(block_model_t *model)
{
	free_dvector(model->a);
	free_dvector(model->inva);
	free_dmatrix(model->b);
	free_dmatrix(model->c);

	free_dvector(model->gx);
	free_dvector(model->gy);
	free_dvector(model->gx_int);
	free_dvector(model->gy_int);
	free_dvector(model->gx_sp);
	free_dvector(model->gy_sp);
	free_dvector(model->gx_hs);
	free_dvector(model->gy_hs);
	free_dvector(model->g_amb);
	free_dvector(model->t_vector);
	free_ivector(model->p);

	free_dmatrix(model->len);
	free_dmatrix(model->g);
	free_dmatrix(model->lu);

	free_imatrix(model->border);

	free(model);
}

void debug_print_block(block_model_t *model)
{
	fprintf(stdout, "printing block model information...\n");
	fprintf(stdout, "n_nodes: %d\n", model->n_nodes);
	fprintf(stdout, "n_units: %d\n", model->n_units);
	fprintf(stdout, "base_n_units: %d\n", model->base_n_units);
	fprintf(stdout, "r_ready: %d\n", model->r_ready);
	fprintf(stdout, "c_ready: %d\n", model->c_ready);

	debug_print_package_RC(&model->pack);

	fprintf(stdout, "printing matrix b:\n");
	dump_dmatrix(model->b, model->n_nodes, model->n_nodes);
	fprintf(stdout, "printing vector a:\n");
	dump_dvector(model->a, model->n_nodes);
	fprintf(stdout, "printing vector inva:\n");
	dump_dvector(model->inva, model->n_nodes);
	fprintf(stdout, "printing matrix c:\n");
	dump_dmatrix(model->c, model->n_nodes, model->n_nodes);
	fprintf(stdout, "printing vector g_amb:\n");
	dump_dvector(model->g_amb, model->n_units+EXTRA);
}

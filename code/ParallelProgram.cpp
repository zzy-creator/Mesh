#include "Mesh.h"
#include "ParallelProgram.h"
#include <mpi.h>

ParallelProgram::ParallelProgram(Mesh& ms) {
	s = ms;
}
ParallelProgram::~ParallelProgram() {

}

void ParallelProgram::initflow() {
	printf("initialising flow field \n");

	gam = 1.4f;
	gm1 = gam - 1.0f;
	cfl = 0.9f;
	eps = 0.05f;

	double mach = 0.4f;
	double alpha = 3.0f * atan(1.0f) / 45.0f;
	double p = 1.0f;
	double r = 1.0f;
	double u = sqrt(gam * p / r) * mach;
	double e = p / (r * gm1) + 0.5f * u * u;

	qinf[0] = r;
	qinf[1] = r * u;
	qinf[2] = 0.0f;
	qinf[3] = r * e;

	return;
}


//kernels

void ParallelProgram::save(const double* q, double* qold) {
	for (int n = 0; n < 4; n++) qold[n] = q[n];//每次只存了四个
}

void ParallelProgram::area(const double* x1, const double* x2, const double* x3,
	const double* x4, const double* q, double* adt) {  //loop on cells
	double dx, dy, ri, u, v, c;

	ri = 1.0f / q[0];
	u = ri * q[1];
	v = ri * q[2];
	c = sqrt(gam * gm1 * (ri * q[3] - 0.5f * (u * u + v * v)));

	dx = x2[0] - x1[0];
	dy = x2[1] - x1[1];
	*adt = fabs(u * dy - v * dx) + c * sqrt(dx * dx + dy * dy);

	dx = x3[0] - x2[0];
	dy = x3[1] - x2[1];
	*adt += fabs(u * dy - v * dx) + c * sqrt(dx * dx + dy * dy);

	dx = x4[0] - x3[0];
	dy = x4[1] - x3[1];
	*adt += fabs(u * dy - v * dx) + c * sqrt(dx * dx + dy * dy);

	dx = x1[0] - x4[0];
	dy = x1[1] - x4[1];
	*adt += fabs(u * dy - v * dx) + c * sqrt(dx * dx + dy * dy);

	*adt = (*adt) / cfl;
}

void ParallelProgram::flux(const double* x1, const double* x2, const double* q1,
	const double* q2, const double* adt1, const double* adt2,
	double* res1, double* res2) {  //loop on edges
	double dx, dy, mu, ri, p1, vol1, p2, vol2, f;

	dx = x1[0] - x2[0];
	dy = x1[1] - x2[1];

	ri = 1.0f / q1[0];
	p1 = gm1 * (q1[3] - 0.5f * ri * (q1[1] * q1[1] + q1[2] * q1[2]));
	vol1 = ri * (q1[1] * dy - q1[2] * dx);

	ri = 1.0f / q2[0];
	p2 = gm1 * (q2[3] - 0.5f * ri * (q2[1] * q2[1] + q2[2] * q2[2]));
	vol2 = ri * (q2[1] * dy - q2[2] * dx);

	mu = 0.5f * ((*adt1) + (*adt2)) * eps;

	f = 0.5f * (vol1 * q1[0] + vol2 * q2[0]) + mu * (q1[0] - q2[0]);
	res1[0] += f;
	res2[0] -= f;
	f = 0.5f * (vol1 * q1[1] + p1 * dy + vol2 * q2[1] + p2 * dy) +
		mu * (q1[1] - q2[1]);
	res1[1] += f;
	res2[1] -= f;
	f = 0.5f * (vol1 * q1[2] - p1 * dx + vol2 * q2[2] - p2 * dx) +
		mu * (q1[2] - q2[2]);
	res1[2] += f;
	res2[2] -= f;
	f = 0.5f * (vol1 * (q1[3] + p1) + vol2 * (q2[3] + p2)) + mu * (q1[3] - q2[3]);
	res1[3] += f;
	res2[3] -= f;
}

void ParallelProgram::bcond(const double* x1, const double* x2, const double* q1,
	const double* adt1, double* res1, const int* bound) {  //loop on bedges
	double dx, dy, mu, ri, p1, vol1, p2, vol2, f;

	dx = x1[0] - x2[0];
	dy = x1[1] - x2[1];

	ri = 1.0f / q1[0];
	p1 = gm1 * (q1[3] - 0.5f * ri * (q1[1] * q1[1] + q1[2] * q1[2]));

	if (*bound == 1) {
		res1[1] += +p1 * dy;
		res1[2] += -p1 * dx;
	}
	else {
		vol1 = ri * (q1[1] * dy - q1[2] * dx);

		ri = 1.0f / qinf[0];
		p2 = gm1 * (qinf[3] - 0.5f * ri * (qinf[1] * qinf[1] + qinf[2] * qinf[2]));
		vol2 = ri * (qinf[1] * dy - qinf[2] * dx);

		mu = (*adt1) * eps;

		f = 0.5f * (vol1 * q1[0] + vol2 * qinf[0]) + mu * (q1[0] - qinf[0]);
		res1[0] += f;
		f = 0.5f * (vol1 * q1[1] + p1 * dy + vol2 * qinf[1] + p2 * dy) +
			mu * (q1[1] - qinf[1]);
		res1[1] += f;
		f = 0.5f * (vol1 * q1[2] - p1 * dx + vol2 * qinf[2] - p2 * dx) +
			mu * (q1[2] - qinf[2]);
		res1[2] += f;
		f = 0.5f * (vol1 * (q1[3] + p1) + vol2 * (qinf[3] + p2)) +
			mu * (q1[3] - qinf[3]);
		res1[3] += f;
	}

}

void ParallelProgram::update(const double* qold, double* q, double* res,
	const double* adt, double* rms) {  //loop on cells
	double del, adti;

	adti = 1.0f / (*adt);

	for (int n = 0; n < 4; n++) {
		del = adti * res[n];
		q[n] = qold[n] - del;
		res[n] = 0.0f;
		*rms += del * del;
	}

}
//
// user declared functions
//

static int compute_local_size(int global_size, int mpi_comm_size,
	int mpi_rank) {
	int local_size = global_size / mpi_comm_size;
	int remainder = (int)fmod(global_size, mpi_comm_size);

	if (mpi_rank < remainder) {
		local_size = local_size + 1;
	}
	return local_size;
}

static void scatter_double_array(double* g_array, double* l_array,
	int comm_size, int g_size, int l_size,
	int elem_size) {
	int* sendcnts = (int*)malloc(comm_size * sizeof(int));
	int* displs = (int*)malloc(comm_size * sizeof(int));
	int disp = 0;

	for (int i = 0; i < comm_size; i++) {
		sendcnts[i] = elem_size * compute_local_size(g_size, comm_size, i);
	}
	for (int i = 0; i < comm_size; i++) {
		displs[i] = disp;
		disp = disp + sendcnts[i];
	}

	MPI_Scatterv(g_array, sendcnts, displs, MPI_DOUBLE, l_array,
		l_size * elem_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	free(sendcnts);
	free(displs);
}

static void scatter_int_array(int* g_array, int* l_array, int comm_size,
	int g_size, int l_size, int elem_size) {
	int* sendcnts = (int*)malloc(comm_size * sizeof(int));
	int* displs = (int*)malloc(comm_size * sizeof(int));
	int disp = 0;

	for (int i = 0; i < comm_size; i++) {
		sendcnts[i] = elem_size * compute_local_size(g_size, comm_size, i);
	}
	for (int i = 0; i < comm_size; i++) {
		displs[i] = disp;
		disp = disp + sendcnts[i];
	}

	MPI_Scatterv(g_array, sendcnts, displs, MPI_INT, l_array, l_size * elem_size,
		MPI_INT, 0, MPI_COMM_WORLD);

	free(sendcnts);
	free(displs);
}

static void gather_double_array(double* g_array, double* l_array,
	int comm_size, int g_size, int l_size,
	int elem_size) {
	int* sendcnts = (int*)malloc(comm_size * sizeof(int));
	int* displs = (int*)malloc(comm_size * sizeof(int));
	int disp = 0;

	for (int i = 0; i < comm_size; i++) {
		sendcnts[i] = elem_size * compute_local_size(g_size, comm_size, i);
	}
	for (int i = 0; i < comm_size; i++) {
		displs[i] = disp;
		disp = disp + sendcnts[i];
	}
	MPI_Gatherv(l_array, l_size * elem_size, MPI_DOUBLE, g_array, sendcnts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	free(sendcnts);
	free(displs);
}

static void gather_int_array(int* g_array, int* l_array, int comm_size,
	int g_size, int l_size, int elem_size) {
	int* sendcnts = (int*)malloc(comm_size * sizeof(int));
	int* displs = (int*)malloc(comm_size * sizeof(int));
	int disp = 0;

	for (int i = 0; i < comm_size; i++) {
		sendcnts[i] = elem_size * compute_local_size(g_size, comm_size, i);
	}
	for (int i = 0; i < comm_size; i++) {
		displs[i] = disp;
		disp = disp + sendcnts[i];
	}
	MPI_Gatherv(l_array, l_size * elem_size, MPI_INT, g_array, sendcnts, displs, MPI_INT, 0, MPI_COMM_WORLD);
	free(sendcnts);
	free(displs);
}

static void check_scan(int items_received, int items_expected) {
	if (items_received != items_expected) {
		printf("error reading from new_grid.dat\n");
		exit(-1);
	}
}

static void compute_local_range(int global_size, int mpi_comm_size,
	int mpi_rank, int& begin, int& end) {
	begin = 0;
	for (int i = 0; i < mpi_rank; i++) {
		begin += compute_local_size(global_size, mpi_comm_size, i);
	}
	end = begin + compute_local_size(global_size, mpi_comm_size, mpi_rank);
}

int ParallelProgram::run(int argc, char* argv[])
{
	MPI_Init(&argc, &argv);

	int my_rank, namelen, comm_size;
	char processor_name[MPI_MAX_PROCESSOR_NAME];

	MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Get_processor_name(processor_name, &namelen);

	//Global Mesh 
	int g_nnode = s.element_list[0]->size;
	int g_ncell = s.element_list[3]->size;
	int g_nedge = s.element_list[1]->size;
	int g_nbedge = s.element_list[2]->size;
	int* g_becell = 0, * g_ecell = 0, * g_bound = 0, * g_bedge = 0, * g_edge = 0,
		* g_cell = 0;
	double* g_x = 0, * g_q = 0, * g_qold = 0, * g_adt = 0, * g_res = 0;

	//local mesh
	int nnode, ncell, nedge, nbedge;
	int* becell, * ecell, * bound, * bedge, * edge, * cell;
	double* x, * q, * qold, * adt, * res;
	double rms, g_rms;

	initflow();//

	g_cell = (int*)(s.map_list[4]->map);	//一个cell 4个点
	g_edge = (int*)(s.map_list[0]->map);	//一条边 2个点
	g_ecell = (int*)(s.map_list[1]->map);	//一个ecell 2个点
	g_bedge = (int*)(s.map_list[2]->map);	//一条b边2个点
	g_becell = (int*)(s.map_list[3]->map);	//一个becell 1个cell

	g_x = (double*)(s.dat_list[1]->data);//一个点三个坐标
	g_q = (double*)(s.dat_list[2]->data);
	g_qold = (double*)(s.dat_list[3]->data);
	g_res = (double*)(s.dat_list[5]->data);
	g_adt = (double*)(s.dat_list[4]->data);
	g_bound = (int*)(s.dat_list[0]->data);

	for (int n = 0; n < g_ncell; n++) { //ncell
		for (int m = 0; m < 4; m++) {
			g_q[4 * n + m] = qinf[m];
			g_res[4 * n + m] = 0.0f;
		}
	}
	//为了最后线程0聚合所有数据
	qold = (double*)malloc(4 * g_ncell * sizeof(double));
	res = (double*)malloc(4 * g_ncell * sizeof(double));
	adt = (double*)malloc(g_ncell * sizeof(double));

	MPI_Barrier(MPI_COMM_WORLD);


	nnode = compute_local_size(g_nnode, comm_size, my_rank);
	ncell = compute_local_size(g_ncell, comm_size, my_rank);
	nedge = compute_local_size(g_nedge, comm_size, my_rank);
	nbedge = compute_local_size(g_nbedge, comm_size, my_rank);

	fprintf(stdout, "Number of nodes, cells, edges, bedges on process %d = %d, %d, %d, %d\n",
		my_rank, nnode, ncell, nedge, nbedge);
	fflush(stdout);

	//partition(cell, ecell, x);  
	int node_begin, node_end, cell_begin, cell_end, edge_begin, edge_end, bedge_begin, bedge_end;
	compute_local_range(g_nnode, comm_size, my_rank, node_begin, node_end);
	compute_local_range(g_ncell, comm_size, my_rank, cell_begin, cell_end);
	compute_local_range(g_nedge, comm_size, my_rank, edge_begin, edge_end);
	compute_local_range(g_nbedge, comm_size, my_rank, bedge_begin, bedge_end);

	int niter = 100;

	rms = 0.0;
	for (int iter = 1; iter <= niter; iter++) {
		for (int i = cell_begin; i < cell_end; i++) {
			save(g_q + 4 * i, g_qold + 4 * i);
		}
		MPI_Barrier(MPI_COMM_WORLD);

		//for (int k = 0; k < 2; k++) {
		for (int i = 0; i < ncell; i++) { //////////////
			area(
				g_x + (g_cell[4 * i]) * 3,
				g_x + (g_cell[4 * i + 1]) * 3,
				g_x + (g_cell[4 * i + 2]) * 3,
				g_x + (g_cell[4 * i + 3]) * 3,
				g_q + 4 * i,
				g_adt + i);
		}
		MPI_Barrier(MPI_COMM_WORLD);

		for (int i = 0; i < nedge; i++) { //////////////
			flux(
				g_x + (g_edge[2 * i]) * 3,
				g_x + (g_edge[2 * i + 1]) * 3,
				g_q + (g_ecell[2 * i]) * 4,
				g_q + (g_ecell[2 * i + 1]) * 4,
				g_adt + (g_ecell[2 * i]),
				g_adt + (g_ecell[2 * i + 1]),
				g_res + (g_ecell[2 * i]) * 4,
				g_res + (g_ecell[2 * i + 1]) * 4);
		}
		MPI_Barrier(MPI_COMM_WORLD);

		for (int i = 0; i < nbedge; i++) {  //////////////////
			bcond(
				g_x + (g_bedge[2 * i]) * 3,
				g_x + (g_bedge[2 * i + 1]) * 3,
				g_q + (g_becell[i]) * 4,
				g_adt + (g_becell[i]),
				g_res + (g_becell[i]) * 4,
				g_bound + i);
		}
		rms = 0.0;
		MPI_Barrier(MPI_COMM_WORLD);

		for (int i = 0; i < ncell; i++) {
			update(
				g_qold + 4 * i,
				g_q + 4 * i,
				g_res + 4 * i,
				g_adt + i,
				&rms);
		}
		//}  //for

		rms = sqrt(rms / (double)ncell);
		printf("Local %d  %10.5e \n", iter, rms);

	}  //for

	MPI_Reduce(&rms, &g_rms, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	//聚合数据到线程0
	gather_double_array(qold, g_qold + 4 * cell_begin, comm_size, g_ncell, ncell, 4);
	gather_double_array(res, g_res + 4 * cell_begin, comm_size, g_ncell, ncell, 4);
	gather_double_array(adt, g_adt + cell_begin, comm_size, g_ncell, ncell, 1);

	if (my_rank == 0) {
		printf("ROOT: Total residual %10.5e \n", g_rms);
	}

	MPI_Finalize();

	
	free(qold);
	free(adt);
	free(res);

	return 0;
}

#include "DistributedProgram.h"
#include "Mesh.h"
#include <mpi.h>

void DistributedProgram::initflow() {
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

DistributedProgram::DistributedProgram(Mesh& ms) {
	s = ms;
}
DistributedProgram::~DistributedProgram() {

}

//kernels

void DistributedProgram::save(const double* q, double* qold) {
	for (int n = 0; n < 4; n++) qold[n] = q[n];//每次只存了四个
}

void DistributedProgram::area(const double* x1, const double* x2, const double* x3,
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

void DistributedProgram::flux(const double* x1, const double* x2, const double* q1,
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

void DistributedProgram::bcond(const double* x1, const double* x2, const double* q1,
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

void DistributedProgram::update(const double* qold, double* q, double* res,
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


static void compute_local_range(int global_size, int mpi_comm_size,
	int mpi_rank, int& begin, int& end) {
	begin = 0;
	for (int i = 0; i < mpi_rank; i++) {
		begin += compute_local_size(global_size, mpi_comm_size, i);
	}
	end = begin + compute_local_size(global_size, mpi_comm_size, mpi_rank);
}

int DistributedProgram::partitionTo(int* cells, int size, double* p_x, double* partition)
{
	int index = 0;//partition的下标
	for (int i = 0;i < size;i++) {//遍历cells
		int* cell = cells + 4 * i;
		for (int j = 0;j < 4;j++) {//遍历每个cell的4个点
			int node = *(cell + j);//node编号
			*(cell + j) = index;
			for (int k = 0;k < 3;k++) {//遍历每个点的xyz
				double data = *(p_x + 3 * node + k);
				partition[3*index+k] = data;
			}
			index++;
		}
	}
	return 0;
}

int DistributedProgram::partitionFrom(int* cells, int size, int size1, int* p_x)
{
	int index = 0;//partition的下标
	for (int i = 0;i < size;i++) {//遍历每个边
		int* edge = cells + 2 * i;
		for (int j = 0;j < 2;j++) {//遍历每个边的node
			int old_node = *(edge + j);//旧的node编号
			//查找新的编号,并更新
			for (int k = 0;k < size1;k++) {
				if (p_x[k] == old_node) *(edge + j) = k;
			}
		}
	}
	return 0;
}


//The following program is logically not correct, just for demo
int DistributedProgram::run(int argc, char* argv[])
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
		* g_cell = 0, * partition_edge, * partition_bedge, * partition_nodes;
	double* g_x = 0, * g_q = 0, * g_qold = 0, * g_adt = 0, * g_res = 0, * partition_x;

	//local mesh
	int nnode, ncell, nedge, nbedge;
	int* becell, * ecell, * bound, * bedge, * edge, * cell;
	double* x, * q, * qold, * adt, * res;
	double rms, g_rms;

	initflow();//

	g_ecell = (int*)(s.map_list[1]->map);	//一个ecell 2个点
	g_becell = (int*)(s.map_list[3]->map);	//一个becell 1个cell

	if (my_rank == 0) {
		g_cell = (int*)(s.map_list[4]->map);	//一个cell 4个点
		g_edge = (int*)(s.map_list[0]->map);	//一条边 2个点
		g_bedge = (int*)(s.map_list[2]->map);	//一条b边2个点

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
		//partition_x = (double*)malloc(3 * 4 * g_ncell * sizeof(double));//用于存储partition后的x，按cell对应

	}//my_rank == 0

	MPI_Barrier(MPI_COMM_WORLD);

	MPI_Bcast(g_ecell, 2 * g_nedge, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(g_becell, g_nbedge, MPI_INT, 0, MPI_COMM_WORLD);

	nnode = compute_local_size(g_nnode, comm_size, my_rank);
	ncell = compute_local_size(g_ncell, comm_size, my_rank);
	nedge = compute_local_size(g_nedge, comm_size, my_rank);
	nbedge = compute_local_size(g_nbedge, comm_size, my_rank);

	fprintf(stdout, "Number of nodes, cells, edges, bedges on process %d = %d, %d, %d, %d\n",
		my_rank, nnode, ncell, nedge, nbedge);
	fflush(stdout);

	partition_x = (double*)malloc(3 * 4 * g_ncell * sizeof(double));//用于存储partition后的x，按cell对应
	partitionTo(g_cell, g_ncell, g_x, partition_x);//将x（nodes）按照cells排列
	if (my_rank == 0) {
		partitionFrom(g_edge, g_nedge, 4 * g_ncell, g_cell);//由于x（nodes）换了新的顺序，所以原来的node编号失效了，要将旧的node编号更换成新的
		partitionFrom(g_bedge, g_nbedge, 4 * g_ncell, g_cell);//由于x（nodes）换了新的顺序，所以原来的node编号失效了，要将旧的node编号更换成新的
	}

	/*Allocate memory to hold local sets, mapping tables and data*/
	cell = (int*)malloc(4 * ncell * sizeof(int));
	edge = (int*)malloc(2 * nedge * sizeof(int));
	//ecell = (int*)malloc(2 * nedge * sizeof(int));
	bedge = (int*)malloc(2 * nbedge * sizeof(int));
	//becell = (int*)malloc(nbedge * sizeof(int));
	bound = (int*)malloc(nbedge * sizeof(int));

	x = (double*)malloc(3 * 4 * ncell * sizeof(double));
	q = (double*)malloc(4 * ncell * sizeof(double));
	qold = (double*)malloc(4 * ncell * sizeof(double));
	res = (double*)malloc(4 * ncell * sizeof(double));
	adt = (double*)malloc(ncell * sizeof(double));

	/* scatter sets, mappings and data on sets*/
	scatter_int_array(g_cell, cell, comm_size, g_ncell, ncell, 4);
	scatter_int_array(g_edge, edge, comm_size, g_nedge, nedge, 2);
	//scatter_int_array(g_ecell, ecell, comm_size, g_nedge, nedge, 2);
	scatter_int_array(g_bedge, bedge, comm_size, g_nbedge, nbedge, 2);
	//scatter_int_array(g_becell, becell, comm_size, g_nbedge, nbedge, 1);
	scatter_int_array(g_bound, bound, comm_size, g_nbedge, nbedge, 1);

	scatter_double_array(partition_x, x, comm_size, g_ncell, ncell, 3 * 4);
	scatter_double_array(g_q, q, comm_size, g_ncell, ncell, 4);
	scatter_double_array(g_qold, qold, comm_size, g_ncell, ncell, 4);
	scatter_double_array(g_res, res, comm_size, g_ncell, ncell, 4);
	scatter_double_array(g_adt, adt, comm_size, g_ncell, ncell, 1);

	//partition(cell, ecell, x);  

	int niter = 100;

	rms = 0.0;
	for (int iter = 1; iter <= niter; iter++) {
		for (int i = 0; i < ncell; i++) {
			save(q + 4 * i, qold + 4 * i);
		}
		MPI_Barrier(MPI_COMM_WORLD);

		for (int i = 0; i < ncell; i++) { //////////////
			area(
				x + 3 * 4 * i,
				x + 3 * 4 * i + 3 * 1,
				x + 3 * 4 * i + 3 * 2,
				x + 3 * 4 * i + 3 * 3,
				q + 4 * i,
				adt + i);
		}
		MPI_Barrier(MPI_COMM_WORLD);

		for (int i = 0; i < nedge; i++) { //////////////
			int ii = 0;
			for (int k = 0;k < my_rank;k++) {
				ii += compute_local_size(g_nedge, comm_size, k);
			}
			flux(
				partition_x + (edge[2 * i]) * 3,
				partition_x + (edge[2 * i + 1]) * 3,
				q + (g_ecell[2 * (ii + i)]) * 4,
				q + (g_ecell[2 * (ii + i) + 1]) * 4,
				adt + (g_ecell[2 * (ii + i)]),
				adt + (g_ecell[2 * (ii + i) + 1]),
				res + (g_ecell[2 * (ii + i)]) * 4,
				res + (g_ecell[2 * (ii + i) + 1]) * 4);
		}
		MPI_Barrier(MPI_COMM_WORLD);

		for (int i = 0; i < nbedge; i++) {  //////////////////
			int ii = 0;
			for (int k = 0;k < my_rank;k++) {
				ii += compute_local_size(g_nbedge, comm_size, k);
			}
			bcond(
				partition_x + (bedge[2 * i]) * 3,
				partition_x + (bedge[2 * i + 1]) * 3,
				q + (g_becell[ii + i]) * 4,
				adt + (g_becell[ii + i]),
				res + (g_becell[ii + i]) * 4,
				bound + i);
		}
		rms = 0.0;
		MPI_Barrier(MPI_COMM_WORLD);

		for (int i = 0; i < ncell; i++) {
			update(
				qold + 4 * i,
				q + 4 * i,
				res + 4 * i,
				adt + i,
				&rms);
		}

		rms = sqrt(rms / (double)ncell);
		printf("Local %d  %10.5e \n", iter, rms);

	}  

	MPI_Reduce(&rms, &g_rms, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);


	if (my_rank == 0) {
		printf("ROOT: Total residual %10.5e \n", g_rms);
	}

	MPI_Finalize();

	free(cell);
	free(edge);
	//free(ecell);
	free(bedge);
	//free(becell);
	free(bound);
	free(x);
	free(q);
	free(qold);
	free(res);
	free(adt);

	return 0;
}

#include "Mesh.h"
#include<map>

Mesh::Mesh() {
	element_list_index = 0;
	map_list_index = 0;
	dat_list_index = 0;
	element_list_size = 10;
	map_list_size = 10;
	dat_list_size = 10;
	element_list = (Elements*)malloc(sizeof(Elements) * element_list_size);
	map_list = (Map*)malloc(sizeof(Map) * map_list_size);
	dat_list = (Data*)malloc(sizeof(Data) * dat_list_size);
	memset(element_list, 0, sizeof(Elements) * map_list_size);
	memset(map_list, 0, sizeof(Map) * element_list_size);
	memset(dat_list, 0, sizeof(Data) * dat_list_size);
}

Mesh::Mesh(char* fname) {
	element_list_index = 0;
	map_list_index = 0;
	dat_list_index = 0;
	element_list_size = 10;
	map_list_size = 10;
	dat_list_size = 10;
	element_list = (Elements*)malloc(sizeof(Elements) * element_list_size);
	map_list = (Map*)malloc(sizeof(Map) * map_list_size);
	dat_list = (Data*)malloc(sizeof(Data) * dat_list_size);
	memset(map_list, 0, sizeof(Map) * map_list_size);
	memset(element_list, 0, sizeof(Elements) * element_list_size);
	memset(dat_list, 0, sizeof(Data) * dat_list_size);
	readfromfile(fname);
}

Mesh::~Mesh() {
	for (int i = 0;i < element_list_index;i++) {
		free(element_list[i]);
	}
	for (int i = 0;i < map_list_index;i++) {
		free(map_list[i]->map);
		free(map_list[i]);
	}
	for (int i = 0;i < dat_list_index;i++) {
		free(dat_list[i]->data);
		free(dat_list[i]);
	}
	free(element_list);
	free(map_list);
	free(dat_list);
}

////////////////////////////////////////////////
void Mesh::init() {
	//清除上一个
	for (int i = 0;i < element_list_index;i++) {
		free(element_list[i]);
	}
	for (int i = 0;i < map_list_index;i++) {
		free(map_list[i]->map);
		free(map_list[i]);
	}
	for (int i = 0;i < dat_list_index;i++) {
		free(dat_list[i]->data);
		free(dat_list[i]);
	}
	free(element_list);
	free(map_list);
	free(dat_list);
	//创建新的
	element_list_index = 0;
	map_list_index = 0;
	dat_list_index = 0;
	element_list_size = 10;
	map_list_size = 10;
	dat_list_size = 10;
	element_list = (Elements*)malloc(sizeof(Elements) * element_list_size);
	map_list = (Map*)malloc(sizeof(Map) * map_list_size);
	dat_list = (Data*)malloc(sizeof(Data) * dat_list_size);
	memset(map_list, 0, sizeof(Map) * map_list_size);
	memset(element_list, 0, sizeof(Elements) * element_list_size);
	memset(dat_list, 0, sizeof(Data) * dat_list_size);
}

Elements Mesh::makeElements(int ele_size, char const* ele_name) {
	Elements ele = (Elements)malloc(sizeof(elements));
	strcpy((char*)(ele->name), ele_name);
	ele->size = ele_size;
	ele->index = element_list_index;
	element_list[element_list_index++] = ele;
	return ele;
}
Map Mesh::makeMap(Elements map_from, Elements map_to, int map_dim, int* map_map, char const* map_name) {
	int arr_cnt = map_from->size;
	map_list[map_list_index] = (Map)malloc(sizeof(map));
	map_list[map_list_index]->index = map_list_index;
	map_list[map_list_index]->from = map_from;
	map_list[map_list_index]->to = map_to;
	map_list[map_list_index]->dim = map_dim;
	//map_list[map_list_index]->map = map_map;
	map_list[map_list_index]->map = (int*)malloc(sizeof(int) * map_dim * arr_cnt);
	for (int i = 0; i < arr_cnt * map_dim; i++)
		map_list[map_list_index]->map[i] = map_map[i];
	strcpy((char*)(map_list[map_list_index]->name), map_name);
	map_list_index++;
	return map_list[map_list_index - 1];

}
Data Mesh::makeData(Elements data_set, int data_dim, char const* data_type, char* data_data, char const* data_name) {
	Data data = (Data)malloc(sizeof(dat));
	data->index = dat_list_index;
	data->set = data_set;
	data->dim = data_dim;
	//data->size = 8;   //xxxxxxxxxxx    ///////////////////////////
	//data->data = data_data;
	if (!strcmp(data_type, "int")) {
		data->size = 4;
		int* t = (int*)malloc(sizeof(int) * data->set->size * data->dim);
		for (int j = 0; j < data->set->size * data->dim; j++)
			t[j] = ((int*)data_data)[j];
		data->data = (char*)t;
	}
	if (!strcmp(data_type, "double")) {
		data->size = 8;
		double* t = (double*)malloc(sizeof(double) * data->set->size * data->dim);
		for (int j = 0; j < data->set->size * data->dim; j++)
			t[j] = ((double*)data_data)[j];
		data->data = (char*)t;
	}
	strcpy((char*)(data->name), data_name);
	strcpy((char*)(data->type), data_type);
	dat_list[dat_list_index++] = data;
	return data;
}

/////////////////////////////////////////////////
bool Mesh::readfromfileuser(const char* fileName) {
	return 0;
}

bool Mesh::readfromfile(const char* fileName) {     /////////////////////////////
	printf("start read:%s\n", fileName);
	FILE* fp;
	if ((fp = fopen(fileName, "rb")) == NULL) {
		printf("can't open file\n");
		return 0;
	}
	init();
	readHeader(fp);
	readElements(fp);
	readMaps(fp);
	readData(fp);

	fclose(fp);
	return true;
}

bool Mesh::savetofile(const char* fileName) {    //////////////////////////
	printf("writing in grid \n");
	FILE* fp;
	if ((fp = fopen(fileName, "wb")) == NULL) {
		printf("can't open file\n");
		return 0;
	}

	writeHeader(fp);
	writeElements(fp);
	writeMaps(fp);
	writeData(fp);

	fclose(fp);
	return true;
}

bool Mesh::writeHeader(FILE* fp) {   //////////////////////
	int header[6] = { element_list_size, map_list_size, dat_list_size,
			element_list_index, map_list_index, dat_list_index };
	int count = fwrite(header, sizeof(int), 6, fp);
	if (count != 6)
		return 0;
	return 1;
}

bool Mesh::writeElements(FILE* fp) {   ////////////////////
	for (int i = 0; i < element_list_index; i++) {
		int count = fwrite(element_list[i], sizeof(elements), 1, fp);
		if (count != 1)
			return 0;
	}
	return 1;
}

bool Mesh::writeMaps(FILE* fp) {    /////////////////////////
	for (int i = 0; i < map_list_index; i++) {
		int temp[5] = { map_list[i]->index, map_list[i]->from->index,
			map_list[i]->to->index, map_list[i]->dim,strlen(map_list[i]->name) };
		int count = fwrite(temp, sizeof(int), 5, fp);
		if (count != 5)
			return 0;
		count = fwrite(map_list[i]->name, sizeof(char), temp[4], fp);
		if (count != temp[4])
			return 0;
		count = fwrite(map_list[i]->map, sizeof(int), map_list[i]->from->size * map_list[i]->dim, fp);
		if (count != map_list[i]->from->size * map_list[i]->dim)
			return 0;
	}
	return 1;
}
bool Mesh::writeData(FILE* fp) {    ////////////////
	for (int i = 0; i < dat_list_index; i++) {
		int temp[5] = { dat_list[i]->index, dat_list[i]->set->index,
			dat_list[i]->dim, dat_list[i]->size,strlen(dat_list[i]->name) };
		int count = fwrite(temp, sizeof(int), 5, fp);
		if (count != 5)
			return 0;
		count = fwrite(dat_list[i]->name, sizeof(char), temp[4], fp);
		if (count != temp[4])
			return 0;
		if (temp[3] == 8) {
			double* temp = (double*)dat_list[i]->data;
			int msize = dat_list[i]->dim * dat_list[i]->set->size;
			int count = fwrite(temp, sizeof(double), msize, fp);
			if (count != msize)
				return 0;
		}
		else if (temp[3] == 4) {
			int count = fwrite(dat_list[i]->data, sizeof(int), dat_list[i]->dim * dat_list[i]->set->size, fp);
			if (count != dat_list[i]->dim * dat_list[i]->set->size)
				return 0;
		}
	}
	return 1;
}

bool Mesh::readHeader(FILE* fp) {   ///////////////////
	int header[6];
	int count = fread(header, sizeof(int), 6, fp);
	if (count != 6)
		return 0;
	element_list_size = header[0];
	map_list_size = header[1];
	dat_list_size = header[2];
	element_list_index = header[3];
	map_list_index = header[4];
	dat_list_index = header[5];
	return 1;
}

bool Mesh::readElements(FILE* fp) {    ////////////////////
	for (int i = 0; i < element_list_index; i++) {
		element_list[i] = (Elements)malloc(sizeof(elements));
		int count = fread(element_list[i], sizeof(elements), 1, fp);
		if (count != 1)
			return 0;
	}
	return 1;
}

bool Mesh::readMaps(FILE* fp) {    ///////////////////////////
	for (int i = 0; i < map_list_index; i++) {
		int temp[5];
		int count = fread(temp, sizeof(int), 5, fp);
		if (count != 5)
			return 0;
		map_list[i] = (Map)malloc(sizeof(map));
		map_list[i]->index = temp[0];
		map_list[i]->from = element_list[temp[1]];
		map_list[i]->to = element_list[temp[2]];
		map_list[i]->dim = temp[3];
		char t[100] = { 0 };
		count = fread(t, sizeof(char), temp[4], fp);
		if (count != temp[4])
			return 0;
		strcpy((char*)(map_list[i]->name), t);
		map_list[i]->map = (int*)malloc(map_list[i]->dim * map_list[i]->from->size * sizeof(int));
		count = fread(map_list[i]->map, sizeof(int), map_list[i]->dim * map_list[i]->from->size, fp);
		if (count != map_list[i]->dim * map_list[i]->from->size)
			return 0;
	}
	return 1;
}
bool Mesh::readData(FILE* fp) {     /////////////////////
	for (int i = 0; i < dat_list_index; i++) {
		int temp[5];
		dat_list[i] = (Data)malloc(sizeof(dat));
		int count = fread(temp, sizeof(int), 5, fp);
		if (count != 5)
			return 0;
		char nam[100] = { 0 };
		count = fread(nam, sizeof(char), temp[4], fp);
		if (temp[4] != count)
			return 0;
		strcpy((char*)(dat_list[i]->name), nam);
		dat_list[i]->index = temp[0];
		dat_list[i]->set = element_list[temp[1]];
		dat_list[i]->dim = temp[2];
		dat_list[i]->size = temp[3];
		if (temp[3] == 4) {
			int* t = (int*)malloc(sizeof(int) * dat_list[i]->dim * dat_list[i]->set->size);
			count = fread(t, sizeof(int), dat_list[i]->dim * dat_list[i]->set->size, fp);
			if (count != dat_list[i]->dim * dat_list[i]->set->size)
				return 0;
			dat_list[i]->data = (char*)t;
		}
		else if (temp[3] == 8) {
			int msize = dat_list[i]->dim * dat_list[i]->set->size;
			double* t = (double*)malloc(sizeof(double) * (msize + 5));
			int count = fread(t, sizeof(double), msize, fp);
			if (count != msize)
				return 0;
			dat_list[i]->data = (char*)t;
		}
	}
	return 1;
}

/////////////////////////////////////// vtk /////////////////////////
bool Mesh::readvtk(const char* fileName) {
	printf("start read:%s\n", fileName);
	FILE* fp;
	if ((fp = fopen(fileName, "r")) == NULL) {
		printf("can't open file\n");
		exit(-1);
	}
	char meshDescription[256];
	readvtkHeader(fp, meshDescription);

	int npoint;
	double* xyz;
	readvtkPoints(fp, npoint, xyz);

	int ncell, cellsize;
	int* cells;
	readvtkCells(fp, ncell, cellsize, cells);

	int* celltypes;
	readvtkCellTypes(fp, ncell, celltypes);

	double* celldata;
	int dimcell, ctype;
	readvtkCellData(fp, ncell, ctype, dimcell, (char*&)celldata);

	int dimpoint, ptype;
	double* pointdata;
	readvtkPointData(fp, npoint, ptype, dimpoint, (char*&)pointdata);

	//获得edge、ecell
	std::map<std::pair<int, int>, int> m1, m2;//存储<边的两个point，边号>
	int en = 0;//边的编号
	for (int i = 0;i < ncell;i++) {//遍历每个cell
		for (int j = 0;j < 4;j++) {
			int n1, n2;//获取边的2个point
			n1 = cells[4 * i + j];
			n2 = cells[4 * i + ((j + 1) % 4)];
			std::pair<int, int> t1(n1, n2), t2(n2, n1);
			if (m1.find(t1) == m1.end() && m1.find(t2) == m1.end()) {//该边是第一次出现
				m1[t1] = en;
				m1[t2] = en;
				en++;
			}
		}
	}

	int nedge = en;
	int* edge = (int*)malloc(2 * nedge * sizeof(int));
	int* ecell = (int*)malloc(2 * nedge * sizeof(int));
	for (int i = 0;i < 2 * nedge;i++) {
		ecell[i] = -1;
	}

	en = 0;
	for (int i = 0;i < ncell;i++) {//遍历每个cell
		for (int j = 0;j < 4;j++) {
			int n1, n2;//获取边的2个point
			n1 = cells[4 * i + j];
			n2 = cells[4 * i + ((j + 1) % 4)];
			std::pair<int, int> t1(n1, n2), t2(n2, n1);
			if (m2.find(t1) != m2.end() || m2.find(t2) != m2.end()) {//该边已出现过，map中已存入该边
				int e = m2[t1];//获得该边的编号
				ecell[2 * e + 1] = i;
			}
			else {//该边第一次出现
				m2[t1] = en;
				m2[t2] = en;
				edge[2 * en] = n1;
				edge[2 * en + 1] = n2;
				ecell[2 * en] = i;
				en++;
			}
		}
	}
	//获得bedge、becell
	//通过：若是边界边，那么对应的ecell只有一个，该边号下第二个ecell为初始值-1，来判断出哪个边号是边界边
	int ben = 0;
	for (int i = 0;i < nedge;i++) {
		if (ecell[2 * i + 1] == -1) ben++;
	}
	int nbedge = ben;
	int* bedge = (int*)malloc(2 * nbedge * sizeof(int));
	int* becell = (int*)malloc(nbedge * sizeof(int));
	int k = 0;
	for (int i = 0;i < nedge;i++) {
		if (ecell[2 * i + 1] == -1) {
			bedge[2 * k] = edge[2 * i];
			bedge[2 * k + 1] = edge[2 * i + 1];
			becell[k] = ecell[2 * i];
			k++;
		}
	}

	Elements enodes = makeElements(npoint, "set_nodes");
	Elements e_edges = makeElements(nedge, "set_edges");
	Elements e_bedges = makeElements(nbedge, "set_bedges");
	Elements ecells = makeElements(ncell, "set_cells");

	makeMap(e_edges, enodes, 2, edge, "edge_to_node_map");
	makeMap(e_edges, ecells, 2, ecell, "edge_to_ecell_map");
	makeMap(e_bedges, enodes, 2, bedge, "bedge_to_node_map");
	makeMap(e_bedges, ecells, 1, becell, "bedge_to_becell_map");
	makeMap(ecells, enodes, 4, cells, "cell_to_node_map");

	makeData(ecells, 1, "int", (char*)celltypes, "celltypes");
	makeData(ecells, 1, "double", (char*)celldata, "data_on_cells");
	makeData(enodes, 3, "double", (char*)xyz, "nodes_xyz");
	makeData(enodes, 1, "double", (char*)pointdata, "data_on_nodes");

	fclose(fp);
	return true;
}


bool Mesh::readvtkHeader(FILE* fp, char* meshDescription) {
	char str[100];
	fscanf(fp, "%[^\n]%*c", str);
	fscanf(fp, "%[^\n]%*c", meshDescription); //
	fscanf(fp, "%[^\n]%*c", str); // ASCII
	if (strcmp(str, "ASCII") != 0) {
		exit(-1);
	}
	fscanf(fp, "%[^\n]%*c", str); //数据信息

	return 0;
}

bool Mesh::readvtkPoints(FILE* fp, int& npoint, double*& xyz) {
	char str[100];
	fscanf(fp, "%s", str);  //POINTS
	fscanf(fp, "%d %s\n", &npoint, str);
	if (strcmp(str, "double") == 0) {
		xyz = (double*)malloc(3 * npoint * sizeof(double));
	}
	for (int n = 0; n < npoint; n++) {
		fscanf(fp, "%lf %lf %lf\n", &(xyz[3 * n]), &(xyz[3 * n + 1]), &(xyz[3 * n + 2]));
	}
	return 0;
}

bool Mesh::readvtkCells(FILE* fp, int& ncell, int& cellsize, int*& cells) {
	char str[100];
	fscanf(fp, "%s", str);  //CELLS
	fscanf(fp, "%d %d\n", &ncell, &cellsize);
	cells = (int*)malloc(cellsize * sizeof(int));
	int count = 0, tt;
	while (count < cellsize) {
		fscanf(fp, "%d\n", &tt);
		int vex = tt;
		//count++;
		while (vex) {
			vex--;
			fscanf(fp, "%d\n", &cells[count++]);
		}
	}
	return 0;
}

bool Mesh::readvtkCellTypes(FILE* fp, int ncell, int*& celltypes) {
	char str[100];
	fscanf(fp, "%s", str);  //CELL_TYPES
	fscanf(fp, "%d\n", &ncell);
	celltypes = (int*)malloc(ncell * sizeof(int));
	for (int n = 0; n < ncell; n++) {
		fscanf(fp, "%d \n", &celltypes[n]);
	}
	return 0;
}

bool Mesh::readvtkPointData(FILE* fp, int npoint, int& type, int& dim, char*& pointdata) {
	char str[100];
	fscanf(fp, "%s", str);  //POINT_DATA
	fscanf(fp, "%d\n", &npoint);
	fscanf(fp, "%s", str);//SCALARS
	fscanf(fp, "%s", str);//CellEntityIds
	fscanf(fp, "%s", str);//double
	if (strcmp(str, "double") == 0) {
		type = 8;
	}
	else if (strcmp(str, "int") == 0) {
		type = 4;
	}
	else type = 1;
	fscanf(fp, "%d\n", &dim);
	pointdata = (char*)malloc(npoint * dim * type * sizeof(char));   //npoint
	//fscanf(fp, "%[^\n]%*c", str);
	//int npComp = str[strlen(str) - 1] - 48;
	fscanf(fp, "%[^\n]%*c", str);
	if (type == 8) {
		for (int n = 0; n < npoint; n++) {
			for (int d = 0; d < dim; d++)
				fscanf(fp, "%lf ", &((double*)pointdata)[n * dim + d]);
			fscanf(fp, "\n");
		}
	}
	else if (type == 4) {
		for (int n = 0; n < npoint; n++) {
			for (int d = 0; d < dim; d++)
				fscanf(fp, "%d ", &((int*)pointdata)[n * dim + d]);
			fscanf(fp, "\n");
		}
	}
	else {
		for (int n = 0; n < npoint; n++) {
			for (int d = 0; d < dim; d++)
				fscanf(fp, "%c ", &pointdata[n * dim + d]);
			fscanf(fp, "\n");
		}
	}
	return 0;
}

bool Mesh::readvtkCellData(FILE* fp, int ncell, int& type, int& dim, char*& celldata) {
	char str[100];
	fscanf(fp, "%s", str);  //CELL_DATA
	fscanf(fp, "%d\n", &ncell);
	fscanf(fp, "%s", str);//SCALARS
	fscanf(fp, "%s", str);//CellEntityIds
	fscanf(fp, "%s", str);//double
	if (strcmp(str, "double") == 0) {
		type = 8;
	}
	else if (strcmp(str, "int") == 0) {
		type = 4;
	}
	else type = 1;
	fscanf(fp, "%d\n", &dim);
	celldata = (char*)malloc(ncell * type * dim * sizeof(char));
	/*fscanf(fp, "%[^\n]%*c", str);
	int ncComp = str[strlen(str) - 1] - 48;*/
	fscanf(fp, "%[^\n]%*c", str);
	if (type == 8) {
		for (int n = 0; n < ncell; n++) {
			for (int d = 0; d < dim; d++)
				fscanf(fp, "%lf ", &((double*)celldata)[n * dim + d]);
			fscanf(fp, "\n");
		}
	}
	else if (type == 4) {
		for (int n = 0; n < ncell; n++) {
			for (int d = 0; d < dim; d++)
				fscanf(fp, "%d ", &((int*)celldata)[n * dim + d]);
			fscanf(fp, "\n");
		}
	}
	else {
		for (int n = 0; n < ncell; n++) {
			for (int d = 0; d < dim; d++)
				fscanf(fp, "%c ", &celldata[n * dim + d]);
			fscanf(fp, "\n");
		}
	}
	return 0;
}


//////////////////////////////////////////////////raw file reading
bool Mesh::readraw(const char* fileName) {
	printf("start read:%s\n", fileName);
	FILE* fp;
	if ((fp = fopen(fileName, "r")) == NULL) {
		printf("can't open file\n");
		return 0;
	}
	int nnode, ncell, nedge, nbedge;
	readrawHeader(fp,&nnode,&ncell,&nedge,&nbedge);
	int* cell = (int*)malloc(4 * ncell * sizeof(int));
	int* edge = (int*)malloc(2 * nedge * sizeof(int));
	int* ecell = (int*)malloc(2 * nedge * sizeof(int));
	int* bedge = (int*)malloc(2 * nbedge * sizeof(int));
	int* becell = (int*)malloc(nbedge * sizeof(int));
	int* bound = (int*)malloc(nbedge * sizeof(int));
	double* x = (double*)malloc(3 * nnode * sizeof(double));
	double*  q = (double*)malloc(4 * ncell * sizeof(double));
	double*  qold = (double*)malloc(4 * ncell * sizeof(double));
	double*  res = (double*)malloc(4 * ncell * sizeof(double));
	double*  adt = (double*)malloc(ncell * sizeof(double));

	readrawNode(fp,nnode,x);
	readrawCell(fp,ncell,cell);
	readrawEdge(fp,nedge,edge,ecell);
	readrawBEdge(fp,nbedge,bedge,becell,bound);
	fclose(fp);

	makeElements(nnode, "nodes"); //0
	makeElements(nedge, "edges"); //1
	makeElements(nbedge, "bedges"); //2
	makeElements(ncell, "cells"); //3
	makeMap(element_list[1], element_list[0], 2, edge, "pedge");
	makeMap(element_list[1], element_list[3], 2, ecell, "pecell");
	makeMap(element_list[2], element_list[0], 2, bedge, "pbedge");
	makeMap(element_list[2], element_list[3], 1, becell, "pbecell");
	makeMap(element_list[3], element_list[0], 4, cell, "pcell");
	makeData(element_list[2], 1, "int", (char*)bound, "p_bound");
	makeData(element_list[0], 3, "double", (char*)x, "p_x");
	makeData(element_list[3], 4, "double", (char*)q, "p_q");
	makeData(element_list[3], 4, "double", (char*)qold, "p_qold");
	makeData(element_list[3], 1, "double", (char*)adt, "p_adt");
	makeData(element_list[3], 4, "double", (char*)res, "p_res");
	
	return true;
}

bool Mesh::readrawHeader(FILE* fp, int* nnode, int* ncell, int* nedge, int* nbedge) {
	return fscanf(fp, "%d %d %d %d \n", nnode, ncell, nedge, nbedge) == 4;
}

bool Mesh::readrawNode(FILE* fp, int nnode, double* x) {
	for (int n = 0; n < nnode; n++) {
		if (fscanf(fp, "%lf %lf \n", &x[3 * n], &x[3 * n + 1]) != 2) {
			printf("error reading from 1new_grid.dat\n"); exit(-1);
		}
		x[3 * n + 2] = 0;
	}
	return 1;
}

bool Mesh::readrawCell(FILE* fp, int ncell, int* cell) {
	for (int n = 0; n < ncell; n++) {
		if (fscanf(fp, "%d %d %d %d \n", &cell[4 * n], &cell[4 * n + 1],
			&cell[4 * n + 2], &cell[4 * n + 3]) != 4) {
			printf("error reading from 2new_grid.dat\n"); exit(-1);
		}
	}
	return 1;
}

bool Mesh::readrawEdge(FILE* fp, int nedge, int* edge, int* ecell) {
	for (int n = 0; n < nedge; n++) {
		if (fscanf(fp, "%d %d %d %d \n", &edge[2 * n], &edge[2 * n + 1],
			&ecell[2 * n], &ecell[2 * n + 1]) != 4) {
			printf("error reading from 3new_grid.dat\n"); exit(-1);
		}
	}
	return 1;
}

bool Mesh::readrawBEdge(FILE* fp, int nbedge, int* bedge, int* becell, int* bound) {
	for (int n = 0; n < nbedge; n++) {
		if (fscanf(fp, "%d %d %d %d \n", &bedge[2 * n], &bedge[2 * n + 1],
			&becell[n], &bound[n]) != 4) {
			printf("error reading from new_grid.dat\n"); exit(-1);
		}
	}
	return 1;
}

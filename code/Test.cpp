#include "Mesh.h"
#include "Program.h"
#include "ParallelProgram.h"
#include "DistributedProgram.h"
#include<iostream>

int main() {
	//test vtk read
	Mesh s;
	s.readvtk("t1.vtk");

	//test emd write
	s.savetofile("new_grid.emd");

	//test emd read
	Mesh s2((char*)"new_grid.emd");

	//test raw read
	Mesh s3;
	s3.readraw("new_grid.dat");

	//test program
	Mesh s4;
	Program prog(s4);
	prog.run(0);

	//test distributed program
	DistributedProgram prog2(s3);
	int argc = 3;
	char* argv[] = { {(char *)"mpiexec"},{(char*)"-n"},{(char*)"4"} };
	prog2.run(argc,argv);

	//test parallel program
	ParallelProgram prog3(s3);
	prog3.run(argc, argv);

	return 0;
}
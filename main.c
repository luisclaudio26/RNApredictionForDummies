#include <stdio.h>
#include <stdlib.h>
#include "nussanov.h"

void fillMainDiagonal(FILE* fileStream, int nBases, MatElem* matrix);
void buildEnergyTable(double* table);

int main(int argc, char** args)
{
	if(argc < 2) return 0;

	//Read bases from file and build matrix
	FILE* f = fopen(args[1], "r");
	int nBases = 0;
	fscanf(f, "%d\n", &nBases);

	MatElem* matrix = (MatElem*)calloc( nBases * nBases , sizeof(MatElem) );
	if(nBases > 0)
		fillMainDiagonal(f, nBases, matrix);

	//Read table of pair energies
	double energies[4*4];
	buildEnergyTable(energies);	

	//do stuff
	nussanov(matrix, nBases, energies);

	//print result
	getResult(matrix, 0, nBases-1, nBases);
	printf("-1 -1");

	//Clean everything
	free(matrix);

	return 0;
}

void buildEnergyTable(double* table)
{
	table[ in(hashBase('A'),hashBase('A'),4) ] = 0.0;
	table[ in(hashBase('A'),hashBase('C'),4) ] = 0.0;
	table[ in(hashBase('A'),hashBase('G'),4) ] = 0.0;
	table[ in(hashBase('A'),hashBase('U'),4) ] = -2.0;
	table[ in(hashBase('C'),hashBase('A'),4) ] = 0.0;
	table[ in(hashBase('C'),hashBase('C'),4) ] = 0.0;
	table[ in(hashBase('C'),hashBase('G'),4) ] = -3.0;
	table[ in(hashBase('C'),hashBase('U'),4) ] = 0.0;
	table[ in(hashBase('G'),hashBase('A'),4) ] = 0.0;
	table[ in(hashBase('G'),hashBase('C'),4) ] = -3.0;
	table[ in(hashBase('G'),hashBase('G'),4) ] = 0.0;
	table[ in(hashBase('G'),hashBase('U'),4) ] = -1.0;
	table[ in(hashBase('U'),hashBase('A'),4) ] = -2.0;
	table[ in(hashBase('U'),hashBase('C'),4) ] = 0.0;
	table[ in(hashBase('U'),hashBase('G'),4) ] = -1.0;
	table[ in(hashBase('U'),hashBase('U'),4) ] = 0.0;
}

void fillMainDiagonal(FILE* fileStream, int nBases, MatElem* matrix)
{
	for(int i = 0; i < nBases; i++)
	{
		char base = 0;
		fscanf(fileStream, " %c\n", &base);

		matrix[in(i,i,nBases)].v = base;
		matrix[in(i,i,nBases)].c = BASE_PAIRED;
		matrix[in(i,i,nBases)].kval = -1;
	}
}
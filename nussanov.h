#ifndef _NUSSANOV_H_
#define _NUSSANOV_H_

#define BASE_PAIRED 0
#define I_UNPAIRED  1
#define J_UNPAIRED  2
#define BIFURCATION 3

#define in(i,j,w) 	((i)*(w)+(j))
#define inDiag(i,w) ((i)*(w)+(i))

typedef struct {
	char v;
	int c;
	int kval;
} MatElem;


int hashBase(char base);
void nussanov(MatElem* matrix, int nBases, double* energies);
void getResult(MatElem* matrix, int i, int j, int n);
void printMatrix(MatElem* m, int n);

#endif
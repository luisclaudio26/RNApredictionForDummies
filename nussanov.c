#include "nussanov.h"
#include <stdio.h>
#include <limits.h>

int hashBase(char base)
{
	switch(base)
	{
		case 'A': return 0;
		case 'C': return 1;
		case 'G': return 2;
		case 'U': return 3;
		default: return -1;
	}
}

void printMatrix(MatElem* m, int n)
{
	printf("     ");
	
	for(int i = 0; i < n; i++) 
		printf("%3d ", i);
	printf("\n");
	
	for(int i = 0; i < n+1; i++)
		printf("----");
	printf("\n");

	for(int i = 0; i < n; i++)
	{
		printf("%2d | ", i);
		for(int j = 0; j < n; j++) 
		{
			if(i == j)
				printf("%3c ", m[in(i,j,n)].v);
			else
				printf("%3d ", m[in(i,j,n)].v);
		}
		printf("\n");
	}
}


static double energyTable(double* energies, char bI, char bJ)
{
	int i = hashBase(bI);
	int j = hashBase(bJ);
	return energies[in(i,j,4)];
}

void nussanov(MatElem* matrix, int n, double* energies)
{
	//Forbidden diagonal
	for(int i = 0; i < n-3; i++) 
	{		
		matrix[ in(i, i+3, n) ].v = 0;
		matrix[ in(i, i+3, n) ].c = -1;
		matrix[ in(i, i+3, n) ].kval = -1;
	}

	//Base case
	for(int i = 0; i < n-4; i++) 
	{
		matrix[ in(i, i+4,n) ].v = energyTable(energies, matrix[inDiag(i,n)].v, matrix[inDiag(i+4,n)].v);
		matrix[ in(i, i+4,n) ].c = BASE_PAIRED;
		matrix[ in(i, i+4,n) ].kval = -1;
	}

	//Fill remaining diagonals
	for(int d = 5; d < n; d++)
	{	
		//Go through diagonal
		for(int i = 0; i < n-d; i++)
		{
			int j = i + d;

			#ifdef VERBOSE
			getc(stdin);
			printf("Calculate (%d,%d)\n", i, j);
			#endif

			//Case 1: (i,j) are paired
			int pairEnergy = energyTable(energies, matrix[inDiag(i,n)].v, matrix[inDiag(j,n)].v );
			int case1 = matrix[in(i+1,j-1,n)].v + pairEnergy;

			#ifdef VERBOSE	
			printf("Case 1: M[%d,%d] (= %d) + (%d)\n", i+1, j-1, matrix[in(i+1,j-1,n)].v, pairEnergy);
			#endif

			//Case 2: i is unpaired
			int case2 = matrix[in(i+1,j,n)].v;

			#ifdef VERBOSE
			printf("Case 2: M[%d,%d] = %d\n", i+1, j, matrix[in(i+1,j,n)].v);
			#endif

			//Case 3: j is unpaired
			int case3 = matrix[in(i,j-1,n)].v;

			#ifdef VERBOSE
			printf("Case 3: M[%d,%d] = %d\n", i, j-1, matrix[in(i,j-1,n)].v);
			#endif

			//Case 4: bifurcation
			int case4 = INT_MAX;
			if(j > i + 9)
				for(int k = i + 4; k <= j-5; k++)
				{
					int temp = matrix[in(i,k,n)].v + matrix[in(k+1,j,n)].v;

					#ifdef VERBOSE
					printf("Case 4.%d: M(%d,%d) + M(%d,%d) = %d\n", k, i, k, k+1, j, temp);
					#endif

					if(temp < case4)
					{
						case4 = temp;
						matrix[in(i,j,n)].kval = k;
					}
				}

			//Get minimun
			int min = INT_MAX;

			if(case4 < min)
			{
				min = case4;
				matrix[in(i,j,n)].c = BIFURCATION;
			}

			if(case1 < min && 
				energyTable(energies, matrix[inDiag(i,n)].v, matrix[inDiag(j,n)].v) < 0)
			{
				min = case1;
				matrix[in(i,j,n)].c = BASE_PAIRED;
			}

			if(case2 < min)
			{
				min = case2;
				matrix[in(i,j,n)].c = I_UNPAIRED;
			}

			if(case3 < min)
			{
				min = case3;
				matrix[in(i,j,n)].c = J_UNPAIRED;
			}

			matrix[in(i,j,n)].v = min;

			#ifdef VERBOSE
			printMatrix(matrix, n);
			#endif
		}
	}
	return;
}

void getResult(MatElem* matrix, int i, int j, int n)
{
	if(i+4 > j) return;

	int c = matrix[in(i,j,n)].c;
	if(c == BASE_PAIRED)
	{
		#ifdef VERBOSE
		printf("M(%d,%d) best conformation is via base pairing\n", i, j);
		#endif

		printf("%d %d\n", i, j);
		getResult(matrix, i+1,j-1,n);
	}

	if(c == I_UNPAIRED)
	{
		#ifdef VERBOSE
		printf("M(%d,%d) best conformation is via leaving %d unpaired\n", i, j, i);
		#endif

		getResult(matrix, i+1, j, n);
	}

	if(c == J_UNPAIRED)
	{
		#ifdef VERBOSE
		printf("M(%d,%d) best conformation is via leaving %d unpaired\n", i, j, j);
		#endif

		getResult(matrix, i, j-1, n);
	}

	if(c == BIFURCATION)
	{
		int k = matrix[in(i,j,n)].kval;
		
		#ifdef VERBOSE
		printf("M(%d,%d) best conformation is via bifurcation in (%d,%d) and (%d,%d) unpaired\n", i, j, i, k, k+1, j);
		#endif

		getResult(matrix, i, k, n);
		getResult(matrix, k+1, j, n);
	}

	return;
}
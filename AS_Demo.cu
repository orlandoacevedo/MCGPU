#include <stdlib.h>
#include <stdio.h>
#define N 100

//Given two indices in an array (representing atoms),
//calculate their product (potential energy),
//and store in energies array.
__global__ void calcEnergy(int *atoms, int numAtoms, int *energies, int numEnergies)
{
	int energyIdx = blockIdx.x * blockDim.x + threadIdx.x;
	
	/*if (energyIdx < numEnergies && threadIdx.x > blockIdx.x)
	{
		energies[energyIdx] = atoms[blockIdx.x] * atoms[threadIdx.x];
	}*/
	if (energyIdx < numEnergies)
		energies[energyIdx] += 1;
}

int main()
{
	int *atomsHost, *atomsDevice, *energiesHost, *energiesDevice,
		atomsSize, energiesSize,
		totalEnergy;
	
	atomsSize = N * sizeof(int);
	energiesSize = sizeof(int) * N * (N - 1) / 2;
	
	atomsHost = (int*) malloc(atomsSize);
	energiesHost = (int*) malloc(energiesSize);
	
	int i;
	for (i = 0; i < N; i++)
	{
		atomsHost[i] = i;
	}
	
	for (i = 0; i < energiesSize / sizeof(int); i++)
	{
		energiesHost[i] = 1;
	}
	
	cudaMalloc(&atomsDevice, atomsSize);
	cudaMalloc(&energiesDevice, energiesSize);
	
	cudaMemcpy(atomsDevice, atomsHost, atomsSize, cudaMemcpyHostToDevice);
	cudaMemcpy(energiesDevice, energiesHost, energiesSize, cudaMemcpyHostToDevice);
	
	//N blocks of N threads (every atom pair)
	calcEnergy<<<N, N>>>(atomsDevice, N, energiesDevice, energiesSize / sizeof(int));
	//calcEnergy<<<10, 10>>>(atomsDevice, N, energiesDevice, energiesSize / sizeof(int));	

	cudaMemcpy(energiesHost, energiesDevice, energiesSize, cudaMemcpyDeviceToHost);
	
	/*for (i = 0; i < N; i++)
	{
		printf("%d\n", atomsHost[i]);
	}*/

	totalEnergy = 0;
	for (i = 0; i < energiesSize / sizeof(int); i++)
	{
		totalEnergy += energiesHost[i];
	}
	
	printf("Total Energy for %d atoms is %d Pseudo-Joules.\n", N, totalEnergy);
	
	return 0;
}
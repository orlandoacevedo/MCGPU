#include <iostream>
#include <time.h>

#define N 50000
#define BLK_SIZE 256

using namespace std;

struct Atom
{
	int x;
	int y;
	int z;
	int a, b, c, d, e, f;
};

__global__ void AtomKernel(Atom *atoms, int *sum);
__global__ void CoalescedKernel(int *x, int *y, int *z, int *sum);

int main()
{
	int i;
	
	//host
	Atom *atoms = (Atom*) malloc(N * sizeof(Atom));
	for (i = 0; i < N; i++)
	{
		atoms[i].x = 10;
		atoms[i].y = 10;
		atoms[i].z = 10;
	}
	int *x = (int*) malloc(N * sizeof(int));
	int *y = (int*) malloc(N * sizeof(int));
	int *z = (int*) malloc(N * sizeof(int));
	memset(x, 10, N * sizeof(int));
	memset(y, 10, N * sizeof(int));
	memset(z, 10, N * sizeof(int));
	
	//device
	Atom *atomsD;
	cudaMalloc(&atomsD, N * sizeof(Atom));
	cudaMemcpy(atomsD, atoms, N * sizeof(Atom), cudaMemcpyHostToDevice);
	int *xD;
	int *yD;
	int *zD;
	int *sumD;
	cudaMalloc(&xD, N * sizeof(int));
	cudaMalloc(&yD, N * sizeof(int));
	cudaMalloc(&zD, N * sizeof(int));
	cudaMalloc(&sumD, N * sizeof(int));
	cudaMemset(xD, 10, N * sizeof(int));
	cudaMemset(yD, 10, N * sizeof(int));
	cudaMemset(zD, 10, N * sizeof(int));
	
	for (i = 0; i < 10000; i++)
	{
		AtomKernel<<<N / BLK_SIZE + 1, BLK_SIZE>>>(atomsD, sumD);
		cudaDeviceSynchronize();
	}
	
	for (i = 0; i < 10000; i++)
	{
		CoalescedKernel<<<N / BLK_SIZE + 1, BLK_SIZE>>>(xD, yD, zD, sumD);
		cudaDeviceSynchronize();
	}
	
	return 0;
}

__global__ void AtomKernel(Atom *atoms, int *sum)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	
	sum[idx] = 0;
	sum[idx] += atoms[idx].x * atoms[idx].x;
	sum[idx] += atoms[idx].y * atoms[idx].y;
	sum[idx] += atoms[idx].z * atoms[idx].z;
}

__global__ void CoalescedKernel(int *x, int *y, int *z, int *sum)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	
	sum[idx] = 0;
	sum[idx] += x[idx] * x[idx];
	sum[idx] += y[idx] * y[idx];
	sum[idx] += z[idx] * z[idx];
}
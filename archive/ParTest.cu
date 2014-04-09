#include <math.h>
#include <time.h>
#include <iostream>
#define BLOCK_SIZE 4
#define CHANGED_MOL 5

using namespace std;

struct Atom
{
    double x;
    double y;
    double z;
    double sigma;
    double epsilon;
    double charge;
};

struct Molecule
{
    Atom *atoms;
	int numOfAtoms;
};

struct Environment
{
    int primaryAtomIndex;
	double x;
	double y;
	double z;
	double cutoff;
};

__global__ void calcInterMolecularEnergy(Molecule *molecules, int currentMol, int numM, Environment *enviro, double *energies, int segmentSize);
__global__ void calcInterAtomicEnergy(Molecule *molecules, int currentMol, int otherMol, Environment *enviro, double *energies, int segmentSize);
__global__ void calcIntraMolecularEnergy(Molecule *molecules, int currentMol, int numE, Environment *enviro, double *energies, int segmentSize);
__device__ double calc_lj(Atom atom1, Atom atom2, double r2);
__device__ double calcCharge(double charge1, double charge2, double r);
__device__ double makePeriodic(double x, double box);
__device__ double calcBlending(double d1, double d2);
__device__ int getXFromIndex(int idx);
__device__ int getYFromIndex(int x, int idx);

int main()
{
	int numA = 100, numM = 10, numEnergies = numA * numA, i, maxMolSize = 0;
	double *energiesH, *energiesD, totalEnergy;
	Atom *atomsH, *atomsD;
	Molecule *moleculesH, *molTrans, *moleculesD;
	Environment *envH, *envD;
	
	energiesH = (double*) malloc(numEnergies*sizeof(double));
	atomsH = (Atom*) malloc(numA*sizeof(Atom));
	moleculesH = (Molecule*) malloc(numM*sizeof(Molecule));
	molTrans = (Molecule*) malloc(numM*sizeof(Molecule));
	envH = (Environment*) malloc(sizeof(Environment));
	
	envH->primaryAtomIndex = 0;
	envH->x = 10;
	envH->y = 10;
	envH->z = 10;
	envH->cutoff = 1.1;
	
	for (i = 0; i < numA; i++)
	{
		atomsH[i].x = i % 10 + 0.5;
		atomsH[i].y = i / 10 + 0.5;
		atomsH[i].z = 0;
		atomsH[i].charge = 1;
		atomsH[i].sigma = 1;
		atomsH[i].epsilon = 1;
	}
	
	cudaMalloc(&atomsD, numA*sizeof(Atom));
	cudaMalloc(&moleculesD, numM*sizeof(Molecule));
	cudaMalloc(&envD, sizeof(Environment));
	cudaMalloc(&energiesD, numEnergies*sizeof(double));
	
	cudaMemcpy(atomsD, atomsH, numA*sizeof(Atom), cudaMemcpyHostToDevice);
	cudaMemcpy(envD, envH, sizeof(Environment), cudaMemcpyHostToDevice);
	
	for (i = 0; i < numEnergies; i++)
	{
		energiesH[i] = 0;
	}
	
	cudaMemcpy(energiesD, energiesH, numEnergies*sizeof(double), cudaMemcpyHostToDevice);
	
	//sets up host molecules (will be pre-done)
	for (i = 0; i < numM; i++)
	{
		moleculesH[i].atoms = atomsH + i * numA / numM;
		moleculesH[i].numOfAtoms = numA / numM;
	}
	
	//sets up device molecules for transfer
	//copies host molecules exactly except for *atoms, which is calculated here
	Atom *a = atomsD;
	for (i = 0; i < numM; i++)
	{
		molTrans[i].atoms = a;
		molTrans[i].numOfAtoms = moleculesH[i].numOfAtoms;
		a += moleculesH[i].numOfAtoms;
		
		if (moleculesH[i].numOfAtoms > maxMolSize)
		{
			maxMolSize = moleculesH[i].numOfAtoms;
		}
	}
	
	cudaMemcpy(moleculesD, molTrans, numM*sizeof(Molecule), cudaMemcpyHostToDevice);
	
	for (i = 0; i < numM; i++)
	{
		moleculesH[i].atoms = atomsH + i * numA / numM;
		moleculesH[i].numOfAtoms = numA / numM;
	}
	
	//calculate intermolecular energies (cutoff check for each molecule)
	calcInterMolecularEnergy<<<numM / BLOCK_SIZE + 1, BLOCK_SIZE>>>(moleculesD, CHANGED_MOL, numM, envD, energiesD, maxMolSize * maxMolSize);
	//calculate intramolecular energies for changed molecule
	int numIntraEnergies = (numA / numM) * (numA / numM - 1) / 2;
	calcIntraMolecularEnergy<<<numIntraEnergies / BLOCK_SIZE + 1, BLOCK_SIZE>>>(moleculesD, CHANGED_MOL, numIntraEnergies, envD, energiesD, maxMolSize * maxMolSize);
						
	cudaMemcpy(energiesH, energiesD, numEnergies*sizeof(double), cudaMemcpyDeviceToHost);
	
	totalEnergy = 0;
	for (i = 0; i < numEnergies; i++)
	{
		totalEnergy += energiesH[i];
	}
	
	cout << "Total energy from molecule " << CHANGED_MOL << ": " << totalEnergy << endl;
}

__global__ void calcInterMolecularEnergy(Molecule *molecules, int currentMol, int numM, Environment *enviro, double *energies, int segmentSize)
{
	int otherMol = blockIdx.x * blockDim.x + threadIdx.x;
	
	if (otherMol < numM && otherMol != currentMol)
	{
		Atom atom1 = molecules[currentMol].atoms[enviro->primaryAtomIndex];
		Atom atom2 = molecules[otherMol].atoms[enviro->primaryAtomIndex];
			
		//calculate distance between atoms
		double deltaX = makePeriodic(atom1.x - atom2.x, enviro->x);
		double deltaY = makePeriodic(atom1.y - atom2.y, enviro->y);
		double deltaZ = makePeriodic(atom1.z - atom2.z, enviro->z);
		//double deltaX = box->makePeriodic(atom1.x - atom2.x, enviro->x);
		//double deltaY = box->makePeriodic(atom1.y - atom2.y, enviro->y);
		//double deltaZ = box->makePeriodic(atom1.z - atom2.z, enviro->z);
		
		double r2 = (deltaX * deltaX) +
					(deltaY * deltaY) + 
					(deltaZ * deltaZ);

		if (r2 < enviro->cutoff * enviro->cutoff)
		{
			//calculate intermolecular energies
			calcInterAtomicEnergy
			<<<molecules[currentMol].numOfAtoms, molecules[otherMol].numOfAtoms>>>
			(molecules, currentMol, otherMol, enviro, energies, segmentSize);
		}
	}
}

__global__ void calcInterAtomicEnergy(Molecule *molecules, int currentMol, int otherMol, Environment *enviro, double *energies, int segmentSize)
{
	Atom atom1 = molecules[currentMol].atoms[blockIdx.x],
		 atom2 = molecules[otherMol].atoms[threadIdx.x];
	//energyIdx will need to take into account different mol sizes (in atoms)
	int energyIdx = otherMol * segmentSize + blockIdx.x * blockDim.x + threadIdx.x;
	
	if (!(currentMol == otherMol && threadIdx.x == blockIdx.x) && atom1.sigma >= 0 && atom1.epsilon >= 0 && atom2.sigma >= 0 && atom2.epsilon >= 0)
	{
		double totalEnergy = 0;
		
		//calculate difference in coordinates
		double deltaX = atom1.x - atom2.x;
		double deltaY = atom1.y - atom2.y;
		double deltaZ = atom1.z - atom2.z;
	  
		//calculate distance between atoms
		deltaX = makePeriodic(deltaX, enviro->x);
		deltaY = makePeriodic(deltaY, enviro->y);
		deltaZ = makePeriodic(deltaZ, enviro->z);
		//deltaX = box->makePeriodic(deltaX, enviro->x);
		//deltaY = box->makePeriodic(deltaY, enviro->y);
		//deltaZ = box->makePeriodic(deltaZ, enviro->z);
		
		double r2 = (deltaX * deltaX) +
			 (deltaY * deltaY) + 
			 (deltaZ * deltaZ);
			
		//gets the fValue if in the same molecule
		double fvalue = 1.0;
		/*if(mol1 == mol2)
		{
			int ** hopTab1 = box->tables[mol1 % box->molecTypenum].hopTable;
			fValue = box->getFValue(i, j, hopTab1);
		}*/
		
		totalEnergy += calc_lj(atom1, atom2, r2) * fvalue;
		totalEnergy += calcCharge(atom1.charge, atom2.charge, sqrt(r2)) * fvalue;
		
		energies[energyIdx] = totalEnergy;
	}
}

__global__ void calcIntraMolecularEnergy(Molecule *molecules, int currentMol, int numE, Environment *enviro, double *energies, int segmentSize)
{
	Molecule cMol = molecules[currentMol];
	int energyIdx = blockIdx.x * blockDim.x + threadIdx.x,
		x = getXFromIndex(energyIdx);
	Atom atom1 = cMol.atoms[x], atom2 = cMol.atoms[getYFromIndex(x, energyIdx)];
	
	if (energyIdx < numE)
	{
		energyIdx += currentMol * segmentSize;
		
		double totalEnergy = 0;
			
		//calculate difference in coordinates
		double deltaX = atom1.x - atom2.x;
		double deltaY = atom1.y - atom2.y;
		double deltaZ = atom1.z - atom2.z;
	  
		//calculate distance between atoms
		deltaX = makePeriodic(deltaX, enviro->x);
		deltaY = makePeriodic(deltaY, enviro->y);
		deltaZ = makePeriodic(deltaZ, enviro->z);
		//deltaX = box->makePeriodic(deltaX, enviro->x);
		//deltaY = box->makePeriodic(deltaY, enviro->y);
		//deltaZ = box->makePeriodic(deltaZ, enviro->z);
		
		double r2 = (deltaX * deltaX) +
			 (deltaY * deltaY) + 
			 (deltaZ * deltaZ);
			
		//gets the fValue if in the same molecule
		double fvalue = 1.0;
		/*if(mol1 == mol2)
		{
			int ** hopTab1 = box->tables[mol1 % box->molecTypenum].hopTable;
			fValue = box->getFValue(i, j, hopTab1);
		}*/
		
		totalEnergy += calc_lj(atom1, atom2, r2) * fvalue;
		totalEnergy += calcCharge(atom1.charge, atom2.charge, sqrt(r2)) * fvalue;
		
		energies[energyIdx] = totalEnergy;
	}
}

__device__ double calc_lj(Atom atom1, Atom atom2, double r2)
{
    //store LJ constants locally
    double sigma = calcBlending(atom1.sigma, atom2.sigma);
    double epsilon = calcBlending(atom1.epsilon, atom2.epsilon);
    
    if (r2 == 0.0)
    {
        return 0.0;
    }
    else
    {
    	//calculate terms
    	const double sig2OverR2 = pow(sigma, 2) / r2;
		const double sig6OverR6 = pow(sig2OverR2, 3);
    	const double sig12OverR12 = pow(sig6OverR6, 2);
    	const double energy = 4.0 * epsilon * (sig12OverR12 - sig6OverR6);
        return energy;
    }

}

__device__ double calcCharge(double charge1, double charge2, double r)
{  
    if (r == 0.0)
    {
        return 0.0;
    }
    else
    {
    	// conversion factor below for units in kcal/mol
    	const double e = 332.06;
        return (charge1 * charge2 * e) / r;
    }
}

__device__ double makePeriodic(double x, double box)
{
    
    while(x < -0.5 * box)
    {
        x += box;
    }

    while(x > 0.5 * box)
    {
        x -= box;
    }

    return x;

}

__device__ double calcBlending(double d1, double d2)
{
    return sqrt(d1 * d2);
}

__device__ int getXFromIndex(int idx)
{
    int c = -2 * idx;
    int discriminant = 1 - 4 * c;
    int qv = (-1 + sqrtf(discriminant)) / 2;
    return qv + 1;
}

__device__ int getYFromIndex(int x, int idx)
{
    return idx - (x * x - x) / 2;
}
#include "baseTests.h"

double calculate_energy(double **coords,  int n_atoms,  double *box_size, double sigma,  double epsilon)
{

    // Loop over all pairs of atoms and calculate
    // the LJ energy
    double total_energy = 0;

    for (int i = 0; i < n_atoms-1; i = i + 1)
    {
        for (int j = i+1; j < n_atoms; j = j + 1)
        {
            double delta_x = coords[j][0] - coords[i][0];
            double delta_y = coords[j][1] - coords[i][1];
            double delta_z = coords[j][2] - coords[i][2];

            // Apply periodic boundaries
            delta_x = make_periodic(delta_x, box_size[0]);
            delta_y = make_periodic(delta_y, box_size[1]);
            delta_z = make_periodic(delta_z, box_size[2]);

             double r2 = (delta_x*delta_x) + (delta_y*delta_y) +
                              (delta_z*delta_z);

            // E_LJ = 4*epsilon[ (sigma/r)^12 - (sigma/r)^6 ]
             double sig2_over_r2 = (sigma*sigma) / r2;
             double sig6_over_r6 = sig2_over_r2*sig2_over_r2*sig2_over_r2;
             double sig12_over_r12 = sig6_over_r6 * sig6_over_r6;

             double e_lj = 4.0 * epsilon * ( sig12_over_r12 - sig6_over_r6 );

            total_energy = total_energy + e_lj;
        }
    }

    // return the total energy of the atoms
    return total_energy;
}

double calculate_energy(Atom *atoms, Environment *enviro, Molecule *molecules)
{
    int atomNumber = enviro->numOfAtoms;

    double totalEnergy = 0;
    
    int i;
    for(i = 0; i < atomNumber - 1; i++)
    {
        int j;
        for(j = i + 1; j < atomNumber; j++)
        {
            double sigma = sqrt(atoms[i].sigma * atoms[j].sigma);
            double epsilon = sqrt(atoms[i].epsilon * atoms[j].epsilon);
    
            double deltaX = atoms[j].x - atoms[i].x;
            double deltaY = atoms[j].y - atoms[i].y;
            double deltaZ = atoms[j].z - atoms[i].z;
          
            deltaX = make_periodic(deltaX, enviro->x);
            deltaY = make_periodic(deltaY, enviro->y);
            deltaZ = make_periodic(deltaZ, enviro->z);

            double r2 = (deltaX * deltaX) +
                              (deltaY * deltaY) + 
                              (deltaZ * deltaZ);

            double r = sqrt(r2);

            double sig2OverR2 = pow(sigma, 2) / r2;
            double sig6OverR6 = pow(sig2OverR2, 3);
            double sig12OverR12 = pow(sig6OverR6, 2);
            double lj_energy = 4.0 * epsilon * (sig12OverR12 - sig6OverR6);

            double charge_energy = calc_charge(atoms[i], atoms[j], *enviro);
            double fValue = 1.0;
            if (molecules != NULL)
            {
                fValue = getFValueLinear(atoms[i], atoms[j], molecules, enviro);
            }
            
            if (r2 == 0.0)
            {
                lj_energy = 0.0;
                charge_energy = 0.0;
            }

            totalEnergy += fValue * (lj_energy + charge_energy);
        }
    }

    return totalEnergy;
}

double make_periodic(double x,  double box)
{
    while (x < -0.5*box)
    {
        x = x + box;
    }

    while (x > 0.5*box)
    {
        x = x - box;
    }

    return x;
}

double wrap_into_box(double x, double box)
{
    while (x > box)
    {
        x = x - box;
    }

    while (x < 0)
    {
        x = x + box;
    }

    return x;
}

long timevaldiff(struct timeval *starttime, struct timeval *finishtime)
{
    long msec;
    msec = (finishtime->tv_sec - starttime->tv_sec)*1000;
    msec += (finishtime->tv_usec - starttime->tv_usec)/1000;
    return msec;
}

double calc_r_value(Atom a1, Atom a2, Environment enviro)
{
    double dx = make_periodic(a1.x - a2.x, enviro.x);
    double dy = make_periodic(a1.y - a2.y, enviro.y);
    double dz = make_periodic(a1.z - a2.z, enviro.z);

    return sqrt(dx * dx + dy * dy + dz * dz);
}

double calc_charge(Atom a1, Atom a2, Environment enviro)
{
    double e = 332.06;  

    return (a1.charge * a2.charge * e) / calc_r_value(a1, a2, enviro);
}

int getMoleculeFromIDLinear(Atom a1, Molecule *molecules, Environment enviro)
{
    int atomId = a1.id;
    int currentIndex = enviro.numOfMolecules - 1;
    int molecId = molecules[currentIndex].id;
    while(atomId < molecId && currentIndex > 0)
    {
        currentIndex -= 1;
        molecId = molecules[currentIndex].id;
    }
    return molecId;

}

double getFValueLinear(Atom atom1, Atom atom2, Molecule *molecules, Environment *enviro)
{
    int m1 = getMoleculeFromIDLinear(atom1, molecules, *enviro);
    int m2 = getMoleculeFromIDLinear(atom2, molecules, *enviro);
    Molecule molec = molecules[0];

    if(m1 != m2)
    {
        return 1.0;
    }
	else
    {
        int hopChain = hopGE3Linear(atom1.id, atom2.id, molecules[m1]);
        if (hopChain == 3)
            return 0.5;
        else if (hopChain > 3)
            return 1.0;
        else
            return 0.0;
     }
}

int hopGE3Linear(int atom1, int atom2, Molecule molecule)
{
    for(int x=0; x< molecule.numOfHops; x++)
    {
		Hop myHop = molecule.hops[x];
		if((myHop.atom1==atom1 && myHop.atom2==atom2) || (myHop.atom1==atom2 && myHop.atom2==atom1))
        {
			return myHop.hop;
        }
	 }
	 return 0;
}

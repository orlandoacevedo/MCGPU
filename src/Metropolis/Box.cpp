/*
	SuperClass to the SerialBox and ParallelBox classes.
*/

#include "Box.h"

Box::Box(IOUtilities ioUtil)
{
	atoms = ioUtil.atompool;
	molecules = ioUtil.molecules;
	environment = ioUtil.currentEnvironment;
	
	bonds = ioUtil.bondpool;
	angles = ioUtil.anglepool;
	dihedrals = ioUtil.dihedralpool;
	hops = ioUtil.hoppool;

	atomCount = environment->numOfAtoms;
	moleculeCount = environment->numOfMolecules;

	std::cout << "# of atoms: " << atomCount << std::endl;
	std::cout << "# of molecules: " << moleculeCount << std::endl;
}

Box::~Box()
{
	FREE(atoms);
	FREE(environment);
	FREE(molecules);
	FREE(energies);
}

int Box::chooseMolecule()
{
	return (int) randomReal(0, environment->numOfMolecules);
}

int Box::changeMolecule(int molIdx)
{
	double maxTranslation = environment->maxTranslation;
	double maxRotation = environment->maxRotation;
		
	saveChangedMole(molIdx);
		
	//Pick an atom in the molecule about which to rotate
	int atomIndex = randomReal(0, molecules[molIdx].numOfAtoms);
	Atom vertex = molecules[molIdx].atoms[atomIndex];

	const double deltaX = randomReal(-maxTranslation, maxTranslation);
	const double deltaY = randomReal(-maxTranslation, maxTranslation);
	const double deltaZ = randomReal(-maxTranslation, maxTranslation);

	const double degreesX = randomReal(-maxRotation, maxRotation);
	const double degreesY = randomReal(-maxRotation, maxRotation);
	const double degreesZ = randomReal(-maxRotation, maxRotation); 

	moveMolecule(molecules[molIdx], vertex, deltaX, deltaY, deltaZ,
		degreesX, degreesY, degreesZ);

	keepMoleculeInBox(&molecules[molIdx], environment);

	return molIdx;
}

void Box::keepMoleculeInBox(Molecule *molecule, Environment *environment)
{		
		for (int j = 0; j < molecule->numOfAtoms; j++)
        {
		    //X axis
			wrapBox(molecule->atoms[j].x, environment->x);
            //Y axis
			wrapBox(molecule->atoms[j].y, environment->y);
            //Z axis
			wrapBox(molecule->atoms[j].z, environment->z);
		}
}

int Box::rollback(int moleno)
{
	return copyMolecule(&molecules[moleno],&changedMol);
}

int Box::saveChangedMole(int moleno)
{
	Molecule *mole_src = &molecules[moleno];

	//free memory of changedMol before allocate memory
	FREE(changedMol.atoms);
	FREE(changedMol.bonds);
	FREE(changedMol.angles);
	FREE(changedMol.dihedrals);
	FREE(changedMol.hops);

	memcpy(&changedMol,mole_src,sizeof(changedMol));

	changedMol.atoms = (Atom *)malloc(sizeof(Atom) * mole_src->numOfAtoms);
	changedMol.bonds = (Bond *)malloc(sizeof(Bond) * mole_src->numOfBonds);
	changedMol.angles = (Angle *)malloc(sizeof(Angle) * mole_src->numOfAngles);
	changedMol.dihedrals = (Dihedral *)malloc(sizeof(Dihedral) * mole_src->numOfDihedrals);
	changedMol.hops = (Hop *)malloc(sizeof(Hop) * mole_src->numOfHops);

	copyMolecule(&changedMol,mole_src);

	return 0;
}

int Box::copyMolecule(Molecule *mole_dst, Molecule *mole_src)
{
    mole_dst->numOfAtoms = mole_src->numOfAtoms;
    mole_dst->numOfBonds = mole_src->numOfBonds;
    mole_dst->numOfAngles = mole_src->numOfAngles;
    mole_dst->numOfDihedrals = mole_src->numOfDihedrals;
    mole_dst->numOfHops = mole_src->numOfHops;
    mole_dst->moleculeIdentificationNumber = mole_src->moleculeIdentificationNumber;
    
    for(int i = 0; i < mole_src->numOfAtoms; i++)
    {
        mole_dst->atoms[i] = mole_src->atoms[i];
    }

    for(int i = 0; i < mole_src->numOfBonds; i++)
    {
        mole_dst->bonds[i] = mole_src->bonds[i];
    }

    for(int i = 0; i < mole_src->numOfAngles; i++)
    {
        mole_dst->angles[i] = mole_src->angles[i];
    }

    for(int i = 0; i < mole_src->numOfDihedrals; i++)
    {
        mole_dst->dihedrals[i] = mole_src->dihedrals[i];
    }
	
    for(int i = 0; i < mole_src->numOfHops; i++)
    {
        mole_dst->hops[i] = mole_src->hops[i];
    }
  
    return 0;  	
}

double Box::wrapBox(double x, double box)
{

    while(x > box)
    {
        x -= box;
    }
    while(x < 0)
    {
        x += box;
    }

    return x;
}
/// @file Box.cpp
///
/// Represents a simulation box, holding environment and molecule data.

#include "Box.h"
#include "Metropolis/Utilities/MathLibrary.h"

Box::Box()
{
	changedMol = Molecule();

	environment = NULL;
	atoms = NULL;
	molecules = NULL;
	energies = NULL;
	
	bonds = NULL;
	angles = NULL;
	dihedrals = NULL;
	hops = NULL;

	atomCount = 0;
	moleculeCount = 0;
}

Box::~Box()
{
	FREE(environment);
	FREE(atoms);
	FREE(molecules);
	FREE(energies);

	FREE(bonds);
	FREE(angles);
	FREE(dihedrals);
	FREE(hops);
}

int Box::chooseMolecule()
{
	return (int) randomReal(0, environment->numOfMolecules);
}

int Box::changeMolecule(int molIdx)
{
	Real maxTranslation = environment->maxTranslation;
	Real maxRotation = environment->maxRotation;
		
	saveChangedMol(molIdx);
		
	//Pick an atom in the molecule about which to rotate
	int atomIndex = randomReal(0, molecules[molIdx].numOfAtoms);
	Atom vertex = molecules[molIdx].atoms[atomIndex];

	const Real deltaX = randomReal(-maxTranslation, maxTranslation);
	const Real deltaY = randomReal(-maxTranslation, maxTranslation);
	const Real deltaZ = randomReal(-maxTranslation, maxTranslation);

	const Real degreesX = randomReal(-maxRotation, maxRotation);
	const Real degreesY = randomReal(-maxRotation, maxRotation);
	const Real degreesZ = randomReal(-maxRotation, maxRotation); 

	moveMolecule(molecules[molIdx], vertex, deltaX, deltaY, deltaZ,
		degreesX, degreesY, degreesZ);

	keepMoleculeInBox(molIdx);

	return molIdx;
}

void Box::keepMoleculeInBox(int molIdx)
{		
		for (int j = 0; j < molecules[molIdx]numOfAtoms; j++)
        {
		    //X axis
			wrapBox(molecules[molIdx]atoms[j].x, environment->x);
            //Y axis
			wrapBox(molecules[molIdx]atoms[j].y, environment->y);
            //Z axis
			wrapBox(molecules[molIdx]atoms[j].z, environment->z);
		}
}

int Box::rollback(int molIdx)
{
	return copyMolecule(&molecules[molIdx],&changedMol);
}

void Box::saveChangedMol(int molIdx)
{
	Molecule *mole_src = &molecules[molIdx];

	//free memory of changedMol before allocate memory
	delete[] changedMol.atoms;
	delete[] changedMol.bonds;
	delete[] changedMol.angles;
	delete[] changedMol.dihedrals;
	delete[] changedMol.hops;

	memcpy(&changedMol,mole_src,sizeof(changedMol));

	changedMol.atoms = new Atom[mole_src->numOfAtoms];
	changedMol.bonds = new Bond[mole_src->numOfBonds];
	changedMol.angles = new Angle[mole_src->numOfAngles];
	changedMol.dihedrals = new Dihedral[mole_src->numOfDihedrals];
	changedMol.hops = new Hop[mole_src->numOfHops];]

	copyMolecule(&changedMol,mole_src);
}

void Box::copyMolecule(Molecule *mol_dst, Molecule *mol_src)
{
    mole_dst->numOfAtoms = mole_src->numOfAtoms;
    mole_dst->numOfBonds = mole_src->numOfBonds;
    mole_dst->numOfAngles = mole_src->numOfAngles;
    mole_dst->numOfDihedrals = mole_src->numOfDihedrals;
    mole_dst->numOfHops = mole_src->numOfHops;
    mole_dst->id = mole_src->id;
    
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
}

Real Box::wrapBox(Real x, Real boxDim)
{

    while(x > boxDim)
    {
        x -= boxDim;
    }
    while(x < 0)
    {
        x += boxDim;
    }

    return x;
}
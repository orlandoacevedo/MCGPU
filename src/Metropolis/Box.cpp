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
		for (int j = 0; j < molecules[molIdx].numOfAtoms; j++)
        {
		    //X axis
			wrapBox(molecules[molIdx].atoms[j].x, environment->x);
            //Y axis
			wrapBox(molecules[molIdx].atoms[j].y, environment->y);
            //Z axis
			wrapBox(molecules[molIdx].atoms[j].z, environment->z);
		}
}

int Box::rollback(int molIdx)
{
	copyMolecule(&molecules[molIdx],&changedMol);

	return molIdx;
}

void Box::saveChangedMol(int molIdx)
{
	Molecule *mol_src = &molecules[molIdx];

	//free memory of changedMol before allocate memory
	delete[] changedMol.atoms;
	delete[] changedMol.bonds;
	delete[] changedMol.angles;
	delete[] changedMol.dihedrals;
	delete[] changedMol.hops;

	memcpy(&changedMol,mol_src,sizeof(changedMol));

	changedMol.atoms = new Atom[mol_src->numOfAtoms];
	changedMol.bonds = new Bond[mol_src->numOfBonds];
	changedMol.angles = new Angle[mol_src->numOfAngles];
	changedMol.dihedrals = new Dihedral[mol_src->numOfDihedrals];
	changedMol.hops = new Hop[mol_src->numOfHops];

	copyMolecule(&changedMol,mol_src);
}

void Box::copyMolecule(Molecule *mol_dst, Molecule *mol_src)
{
    mol_dst->numOfAtoms = mol_src->numOfAtoms;
    mol_dst->numOfBonds = mol_src->numOfBonds;
    mol_dst->numOfAngles = mol_src->numOfAngles;
    mol_dst->numOfDihedrals = mol_src->numOfDihedrals;
    mol_dst->numOfHops = mol_src->numOfHops;
    mol_dst->id = mol_src->id;
    
    for(int i = 0; i < mol_src->numOfAtoms; i++)
    {
        mol_dst->atoms[i] = mol_src->atoms[i];
    }

    for(int i = 0; i < mol_src->numOfBonds; i++)
    {
        mol_dst->bonds[i] = mol_src->bonds[i];
    }

    for(int i = 0; i < mol_src->numOfAngles; i++)
    {
        mol_dst->angles[i] = mol_src->angles[i];
    }

    for(int i = 0; i < mol_src->numOfDihedrals; i++)
    {
        mol_dst->dihedrals[i] = mol_src->dihedrals[i];
    }
	
    for(int i = 0; i < mol_src->numOfHops; i++)
    {
        mol_dst->hops[i] = mol_src->hops[i];
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
/**
 * Represents a simulation box, holding environment and molecule data.
 *
 * Superclass to SerialBox. This class serves as a temporary intermediate class
 * between constants files, z-matrices, config files, and state files, and the
 * SimBox class, which now performs energy calculations.
 */

#include "Box.h"
#include "Metropolis/Utilities/MathLibrary.h"

Box::Box() {
	changedMol = Molecule();

	environment = NULL;
	atoms = NULL;
	molecules = NULL;
	neighborList = NULL;

	bonds = NULL;
	angles = NULL;
	dihedrals = NULL;
	hops = NULL;

	atomCount = 0;
	moleculeCount = 0;
}

Box::~Box() {
	FREE(environment);
	FREE(atoms);
	FREE(molecules);
	FREE(neighborList);

	FREE(bonds);
	FREE(angles);
	FREE(dihedrals);
	FREE(hops);
}

void Box::createNeighborList() {
	neighborList = new NeighborList(molecules, environment);
}

int Box::chooseMolecule() {
	return (int) randomReal(0, environment->numOfMolecules);
}

int Box::changeMolecule(int molIdx) {
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

int Box::changeMolecule(int molIdx, int vIdx, Real dX, Real dY, Real dZ, Real rX, Real rY, Real rZ) {
	Atom vertex = molecules[molIdx].atoms[vIdx];

	saveChangedMol(molIdx);
	moveMolecule(molecules[molIdx], vertex, dX, dY, dZ, rX, rY, rZ);
	keepMoleculeInBox(molIdx);
	return molIdx;
}

void Box::keepMoleculeInBox(int molIdx) {
    int primaryIndex = (*(*(environment->primaryAtomIndexArray))[molecules[molIdx].type])[0];
    Atom primaryAtom = molecules[molIdx].atoms[primaryIndex];

    int positionX = isOutOfBounds(primaryAtom.x, environment->x);
    int positionY = isOutOfBounds(primaryAtom.y, environment->y);
    int positionZ = isOutOfBounds(primaryAtom.z, environment->z);

    for (int i = 0; i < molecules[molIdx].numOfAtoms; i++) {
		//X axis
		molecules[molIdx].atoms[i].x = wrapBox(molecules[molIdx].atoms[i].x, environment->x, positionX);
		//Y axis
		molecules[molIdx].atoms[i].y = wrapBox(molecules[molIdx].atoms[i].y, environment->y, positionY);
		//Z axis
		molecules[molIdx].atoms[i].z = wrapBox(molecules[molIdx].atoms[i].z, environment->z, positionZ);
    }
}

int Box::isOutOfBounds(Real coor, Real boxDim) {
    if (coor < 0)
	return BELOW_ZERO;
    else if (coor > boxDim)
	return ABOVE_BOX_DIM;

    return IN_BOX;
}

int Box::rollback(int molIdx) {
	copyMolecule(&molecules[molIdx],&changedMol);

	return molIdx;
}

void Box::saveChangedMol(int molIdx) {
	Molecule *mol_src = &molecules[molIdx];

	//free memory of changedMol before allocate memory
	if (changedMol.numOfAtoms != 0) {
		delete[] changedMol.atoms;
		delete[] changedMol.bonds;
		delete[] changedMol.angles;
		delete[] changedMol.dihedrals;
		delete[] changedMol.hops;
	}

	memcpy(&changedMol,mol_src,sizeof(changedMol));

	changedMol.atoms = new Atom[mol_src->numOfAtoms];
	changedMol.bonds = new Bond[mol_src->numOfBonds];
	changedMol.angles = new Angle[mol_src->numOfAngles];
	changedMol.dihedrals = new Dihedral[mol_src->numOfDihedrals];
	changedMol.hops = new Hop[mol_src->numOfHops];

	copyMolecule(&changedMol,mol_src);
}

void Box::copyMolecule(Molecule *mol_dst, Molecule *mol_src) {
  mol_dst->numOfAtoms = mol_src->numOfAtoms;
  mol_dst->numOfBonds = mol_src->numOfBonds;
  mol_dst->numOfAngles = mol_src->numOfAngles;
  mol_dst->numOfDihedrals = mol_src->numOfDihedrals;
  mol_dst->numOfHops = mol_src->numOfHops;
  mol_dst->id = mol_src->id;
  mol_dst->type = mol_src->type;

  for(int i = 0; i < mol_src->numOfAtoms; i++) {
    mol_dst->atoms[i] = mol_src->atoms[i];
  }

  for(int i = 0; i < mol_src->numOfBonds; i++) {
    mol_dst->bonds[i] = mol_src->bonds[i];
  }

  for(int i = 0; i < mol_src->numOfAngles; i++) {
    mol_dst->angles[i] = mol_src->angles[i];
  }

  for(int i = 0; i < mol_src->numOfDihedrals; i++) {
    mol_dst->dihedrals[i] = mol_src->dihedrals[i];
  }

  for(int i = 0; i < mol_src->numOfHops; i++) {
    mol_dst->hops[i] = mol_src->hops[i];
  }
}

Real Box::wrapBox(Real x, Real boxDim, int position) {
	if (position == IN_BOX)
		return x;
	else if(position == ABOVE_BOX_DIM)
		x -= boxDim;
	else if (position == BELOW_ZERO)
		x += boxDim;

  return x;
}

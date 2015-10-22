//using namespace std;

/*!\file
  \brief Structures and functions used to do geometric transformations.
  \author Alexander Luchs, Riley Spahn, Seth Wooten

 */

#include "MathLibrary.h"

//_________________________________________________________________________________________________________________
//  Specific namespace/using requirements
//_________________________________________________________________________________________________________________
		//to avoid "using namespace std;"
using std::vector;
using std::string;
using std::map;
using std::stringstream;
using std::cout;


/**
  Global stringstream used to write data to the output Log file
*/
stringstream output;

void seed(int seed) {
	srand(seed);
}

Real randomReal(const Real start, const Real end) {
	return (end-start) * ((Real) rand() / RAND_MAX) + start;
}

Point createPoint(double X, double Y, double Z) {
  Point p;
  p.x = X;
  p.y = Y;
  p.z = Z;
  return p;
}

Plane createPlane(Atom a1, Atom a2, Atom a3) {
  Plane p;
  p.atom1 = a1;
  p.atom2 = a2;
  p.atom3 = a3;
  return p;
}

void printPoint(Point p) {
  printf("%f, %f, %f\n", p.x, p.y, p.z);
}

Atom getAtom(vector<Atom> atoms, unsigned long atomID) {
  for (int i = 0; i < atoms.size(); i++) {
    if (atoms[i].id == atomID) {
      return atoms[i];
    }
  }
	unsigned long inval = (unsigned long) -1;
  return createAtom(inval,-1,-1,-1);
}

Bond getBond(const vector<Bond>& bonds, unsigned long a1, unsigned long a2) {
  for (int i = 0; i < bonds.size(); i++) {
    if (bonds[i].atom1 == a1 && bonds[i].atom2 == a2) {
      return bonds[i];
    }
  }
  return Bond(-1, -1, -1.0, false);
}

vector<unsigned long> getAllBonds(vector<Bond> bonds, unsigned long atomID) {
  vector<unsigned long> toReturn;

  for (int i = 0; i < bonds.size(); i++) {
    unsigned long oppositeAtom = getOppositeAtom(bonds[i], atomID);

		unsigned long inval = (unsigned long) -1;
    if (oppositeAtom != inval) {
      toReturn.push_back(oppositeAtom);
    }
  }
  return toReturn;
}

vector<unsigned long> getIntersection(vector<unsigned long> v1, vector<unsigned long> v2) {

  vector<unsigned long> intersection;

  //not effecient but I will be working with small data sets.
  for (int i = 0; i < v1.size(); i++) {
    for (int j = 0; j < v2.size(); j++) {
      if (v1[i] == v2[j]) {
        intersection.push_back(v1[i]);
      }
    }
  }
  return intersection;
}

bool isMember(vector<unsigned long> atoms, unsigned long toCheck) {
  for (int i = 0; i < atoms.size(); i++) {
    if(atoms[i] == toCheck) {
      return true;
    }
  }
  return false;
}

double degreesToRadians(double degrees) {
  return degrees * PI / 180.0;
}

double radiansToDegrees(double radians) {
  return radians * 180 / PI;
}

unsigned long getOppositeAtom(Bond bond, unsigned long atomID) {
  if (bond.atom1 == atomID) {
    return bond.atom2;
  } else if(bond.atom2 == atomID) {
    return bond.atom1;
  } else {
    return (unsigned long) -1;
  }
}

unsigned long getOppositeAtom(Angle angle, unsigned long atomID) {
  if (angle.atom1 == atomID) {
    return angle.atom2;
  } else if (angle.atom2 == atomID) {
    return angle.atom1;
  } else {
    return (unsigned long) -1;
  }
}

unsigned long getOppositeAtom(Dihedral dihedral, unsigned long atomID) {
  if (dihedral.atom1 == atomID) {
    return dihedral.atom2;
  } else if (dihedral.atom2 == atomID) {
    return dihedral.atom1;
  } else {
    return (unsigned long) -1;
  }
}

unsigned long getCommonAtom(vector<Bond> bonds, unsigned long atom1, unsigned long atom2) {
  vector<unsigned long> atom1Bonds; // atom ids bonded to atom1
  vector<unsigned long> atom2Bonds; // atom ids bonded to atom2

  for (int i = 0; i < bonds.size(); i++) {
    unsigned long opp1 = getOppositeAtom(bonds[i], atom1);
    unsigned long opp2 = getOppositeAtom(bonds[i], atom2);
		unsigned long inval = (unsigned long) -1;


    if (opp1 != inval) {
      atom1Bonds.push_back(opp1);
    }
    if (opp2 != inval) {
      atom2Bonds.push_back(opp2);
    }
  }

  for (int i = 0; i < atom1Bonds.size(); i++) {
    unsigned long currentAtom1 = atom1Bonds[i];

    for (int j = 0; j < atom2Bonds.size(); j++) {
      if (currentAtom1 == atom2Bonds[j]) {
        return currentAtom1;
      }
    }
  }

  return (unsigned long) -1;
}

bool inXZPlane(Atom atom) {
  if (fabs(atom.y) < .00001) {
    return true;
  } else {
    return false;
  }
}

double getDistance(Atom atom1, Atom atom2) {
  return sqrt( pow(atom1.x - atom2.x, 2) +
               pow(atom1.y - atom2.y, 2) +
               pow(atom1.z - atom2.z, 2));
}

double getAngle(Atom atom1, Atom atom2, Atom atom3) {
  // uses law of cosines, (d3_1)^2 = (d1_2)^2 + (d2_3)^2 - 2(d1_2)(d2_3)*cos(theta)
  // theta = arccos(((d3_1)^2 - (d1_2)^2)/(2(d1_2)(d2_3)))
  double d1_2 = getDistance(atom1, atom2); // distance atom1 to atom2
  double d2_3 = getDistance(atom2, atom3); // distance atom2 to atom3
  double d3_1 = getDistance(atom3, atom1); // distance atom3 to atom1

  double numerator = pow(d3_1, 2) - pow(d1_2, 2) - pow(d2_3, 2);
  double denominator = -2 * d1_2 * d2_3;

  double radians = acos(numerator / denominator);

  double degrees = radiansToDegrees(radians);

  if (degrees > 180) {
    degrees -= 180;
  }

  return degrees;
}

Point getNormal(Plane p) {
  //normal vector = (point2 - point1) CROSS (point3 - point1)
  double atomAX = p.atom2.x - p.atom1.x;
  double atomAY = p.atom2.y - p.atom1.y;
  double atomAZ = p.atom2.z - p.atom1.z;
  Point atomA = createPoint(atomAX, atomAY, atomAZ);

  double atomBX = p.atom3.x - p.atom1.x;
  double atomBY = p.atom3.y - p.atom1.y;
  double atomBZ = p.atom3.z - p.atom1.z;
  Point atomB = createPoint(atomBX, atomBY, atomBZ);

  //calculate the cross product of atomA and atomB
  double normX = atomA.y * atomB.z - atomA.z * atomB.y;
  double normY = atomA.z * atomB.x - atomA.x * atomB.z;
  double normZ = atomA.x * atomB.y - atomA.y * atomB.x;
  return createPoint(normX, normY, normZ);
}


double getDistance(Atom atom0, Atom atom1, Atom atom2, Atom atom3) {

  if (fabs(getAngle(atom1, atom2, atom3)) < .00001) {
    //determine if line passes through origin
    // see : http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
		unsigned long inval = (unsigned long) -1;

    Atom a1_a0 = createAtom(inval, atom1.x - atom0.x, atom1.y - atom0.y,
        atom1.z - atom0.z);
    Atom a2_a1 = createAtom(inval, atom2.x - atom1.x, atom2.y - atom1.y,
        atom2.z - atom1.z);

    //calc dot product of a2_a0 a2_a1
    double numerator = (a1_a0.x * a2_a1.x) + (a1_a0.y * a2_a1.y) +
        (a1_a0.z * a2_a1.z);

    //find the square of the madnitude of a2 - a1 (distanc between atom1 and atom2)
    double denom = pow(getDistance(atom2, atom1), 2);

    double t = -1 * numerator / denom;

    //find the distance^2 from the origin to the line generated by the atoms.
    double distance_2_a0 = pow((atom1.x - atom0.x) + (atom2.x - atom1.x) * t, 2) +
        pow((atom1.y - atom0.y) + (atom2.y - atom1.y) * t, 2) +
        pow((atom1.z - atom0.z) + (atom2.z - atom1.z) * t, 2);

    return sqrt(distance_2_a0);
  } else {
    // atom1 atom2 atom3 do not define a line
    return -1;
  }
}

Point getNormal(Atom atom1, Atom atom2, Atom atom3) {
  //atoms are co-linear
  if (fabs(getAngle(atom1, atom2, atom3)) < .00001) {
    //will try to use 0,0,0 then 1,0,0 then 0,1,0 as normals
		unsigned long inval = (unsigned long) -1;
    //determine if line passes through origin
    Atom atom0 = createAtom(inval, 0, 0, 0);
    double distance_2_origin = getDistance(atom0, atom1, atom2, atom3);

    atom0 = createAtom(inval, 1,0,0);
    double distance_2_x1 = getDistance(atom0, atom1, atom2, atom3);

    atom0 = createAtom(inval, 0, 1, 0);
    double distance_2_y1 = getDistance(atom0, atom1, atom2, atom3);

    if (distance_2_origin > .00001) {
      // use origin
      atom0 = createAtom(inval, 0, 0, 0);
    } else if (distance_2_x1 > .00001) {
      // use 1,0,0
      atom0 = createAtom(inval, 1, 0, 0);
    } else {
      // use 0,1,0
      atom0 = createAtom(inval, 0, 1, 0);
    }

    Plane p = createPlane(atom0, atom2, atom1);
    return getNormal(p);
  }

  //Needs a return of some type
  return Point(); //This is just a temp filler
}

double getAngle(Plane p1, Plane p2) {
  //the normal vectors for each plane defined by a point
  Point normal1 = getNormal(p1);
  Point normal2 = getNormal(p2);

	unsigned long inval = (unsigned long) -1;

  Atom a1 = createAtom(inval, normal1.x, normal1.y, normal1.z);
  Atom a2 = createAtom(inval, 0, 0, 0);
  Atom a3 = createAtom(inval, normal2.x, normal2.y, normal2.z);

  return getAngle(a1, a2, a3);
}

Atom translateAtom(Atom atom, double x, double y, double z) {
  atom.x += x;
  atom.y += y;
  atom.z += z;

  return atom;
}

Atom rotateAboutX(Atom atom, double theta) {
  double thetaRadians = degreesToRadians(theta);
  Atom returnAtom = atom;
  returnAtom.y = atom.y * cos(thetaRadians) + atom.z * sin(thetaRadians);
  returnAtom.z = atom.z * cos(thetaRadians) - atom.y * sin(thetaRadians);
  return returnAtom;
}

Atom rotateAboutY(Atom atom, double theta) {
  double thetaRadians = degreesToRadians(theta);
  Atom returnAtom = atom;
  returnAtom.z = atom.z * cos(thetaRadians) + atom.x * sin(thetaRadians);
  returnAtom.x = atom.x * cos(thetaRadians) - atom.z * sin(thetaRadians);
  return returnAtom;
}

Atom rotateAboutZ(Atom atom, double theta) {
  double thetaRadians = degreesToRadians(theta);
  Atom returnAtom = atom;
  returnAtom.x = atom.x * cos(thetaRadians) + atom.y * sin(thetaRadians);
  returnAtom.y = atom.y * cos(thetaRadians) - atom.x * sin(thetaRadians);
  return returnAtom;
}

Atom rotateAtomInPlane(Atom atom1, Atom atom2, Atom atom3, double theta) {
  //Create a plane.
  Plane atomPlane = createPlane(atom1, atom2, atom3);
  //Find the normal vector.
  Point normal;
  //Arbitrarily assigned percent diff.
  if(fabs(getAngle(atom1, atom2, atom3)) < .0001) {
    //find a line that is perpendicular to the line of atoms as the normal.
    normal = getNormal(atom1, atom2, atom3);
  } else {
    normal = getNormal(atomPlane);
  }

	unsigned long inval = (unsigned long) -1;
  Atom vectorEnd = createAtom(inval, atom2.x + normal.x, atom2.y + normal.y,
      atom2.z + normal.z);

  //Rotate about that normal vector
  return rotateAtomAboutVector(atom1, atom2, vectorEnd, theta);
}

Atom rotateAtomAboutVector(Atom atom1, Atom atom2, Atom atom3, double theta) {
  theta *= -1;

  //Translate all atoms so that atom2 is at the origin.
  //The rotation axis needs to pass through the origin
  Atom originalAtom2 = atom2;

  atom1 = translateAtom(atom1, -originalAtom2.x, -originalAtom2.y, -originalAtom2.z);
  atom2 = translateAtom(atom2, -originalAtom2.x, -originalAtom2.y, -originalAtom2.z);
  atom3 = translateAtom(atom3, -originalAtom2.x, -originalAtom2.y, -originalAtom2.z);

  //find the angle between the vector and xz plane
  double xzAngle = 0.0;

  if (!inXZPlane(atom2) || !inXZPlane(atom3)) {
    // rotation axis not in xz plane
    //determin if the rotation axis is vertical
    Atom xzVector;

		unsigned long inval = (unsigned long) -1;

    if (compareDoubleDifference(atom2.x, atom3.x, DOUBPREC)
        && compareDoubleDifference(atom2.z, atom3.z, DOUBPREC)) {

      xzVector = createAtom(inval, atom3.x + 1, atom2.y, atom3.z);
    } else {
      xzVector = createAtom(inval, atom3.x, 0, atom3.z);
    }
    xzAngle = getAngle(atom3, atom2, xzVector);

    //rotate about z axis so that vector is parallel to xz plane
    atom1 = rotateAboutZ(atom1, xzAngle);
    //atom2 should not change because atom2 is at the origin
    atom2 = rotateAboutZ(atom2, xzAngle);
    atom3 = rotateAboutZ(atom3, xzAngle);
  }

	unsigned int inval = (unsigned int) -1;

  //find the angle between the vector and the z axis
  Atom zAxis = createAtom(inval, 0, 0, 1);
  double zAngle = getAngle(atom3, atom2, zAxis);

  //rotate about y axis so that the vector is parallel to z axis
  atom1 = rotateAboutY(atom1, zAngle);
  atom2 = rotateAboutY(atom2, zAngle);
  atom3 = rotateAboutY(atom3, zAngle);

  //rotate atom1 theta about the z axis.
  atom1 = rotateAboutZ(atom1, theta);

  //invert rotation about y axis
  atom1 = rotateAboutY(atom1, -zAngle);
  atom2 = rotateAboutY(atom2, -zAngle);
  atom3 = rotateAboutY(atom3, -zAngle);

  //invert rotation about z axis
  atom1 = rotateAboutZ(atom1, -xzAngle);
  atom2 = rotateAboutZ(atom2, -xzAngle);
  atom3 = rotateAboutZ(atom3, -xzAngle);

  //invert translation to origin
  atom1 = translateAtom(atom1, originalAtom2.x, originalAtom2.y, atom3.z);
  atom2 = translateAtom(atom2, originalAtom2.x, originalAtom2.y, atom3.z);
  atom3 = translateAtom(atom3, originalAtom2.x, originalAtom2.y, atom3.z);

  //the inversions for atoms 2 and 3 are not neccesary b/c of copy by value.
  return atom1;
}

bool compareDoubleDifference(double a, double b, double precision) {
  if (fabs(a - b) < precision) {
    return true;
  } else {
    return false;
  }
}

Molecule moveMolecule(Molecule molec, Atom pivot, double xTrans, double yTrans,
        double zTrans, double xRot, double yRot, double zRot) {

  for (int i = 0; i < molec.numOfAtoms; i++) {
    //translate molecule to the origin to rotate
    molec.atoms[i] = translateAtom(molec.atoms[i], -pivot.x, -pivot.y, -pivot.z);
    //rotateAboutX
    molec.atoms[i] = rotateAboutX(molec.atoms[i], xRot);
    //rotateAboutY
    molec.atoms[i] = rotateAboutY(molec.atoms[i], yRot);
    //rotateAboutZ
    molec.atoms[i] = rotateAboutZ(molec.atoms[i], zRot);
    //translate to original position
    molec.atoms[i] = translateAtom(molec.atoms[i], pivot.x, pivot.y, pivot.z);

    //translate atom to final position
    molec.atoms[i] = translateAtom(molec.atoms[i], xTrans, yTrans, zTrans);
  }

  return molec;
}

double randomNUM(const double start, const double end) {
  return (end-start) * (double(rand()) / RAND_MAX) + start;
}

Molecule translateMolecule(Molecule molec, double maxTranslation) {
  const double xTrans = randomNUM(-maxTranslation, maxTranslation);
	const double yTrans = randomNUM(-maxTranslation, maxTranslation);
	const double zTrans = randomNUM(-maxTranslation, maxTranslation);

  for (int i = 0; i < molec.numOfAtoms; i++) {
    molec.atoms[i] = translateAtom(molec.atoms[i], xTrans, yTrans, zTrans);
  }
  return molec;
}

Molecule rotateMolec(Molecule molec, Atom pivot, double maxRotation) {
	const double xRot = randomNUM(-maxRotation, maxRotation);
	const double yRot = randomNUM(-maxRotation, maxRotation);
	const double zRot = randomNUM(-maxRotation, maxRotation);

  for (int i = 0; i < molec.numOfAtoms; i++) {
    //translate molecule to the origin to rotate
    molec.atoms[i] = translateAtom(molec.atoms[i], -pivot.x, -pivot.y, -pivot.z);
    //rotateAboutX
    molec.atoms[i] = rotateAboutX(molec.atoms[i], xRot);
    //rotateAboutY
    molec.atoms[i] = rotateAboutY(molec.atoms[i], yRot);
    //rotateAboutZ
    molec.atoms[i] = rotateAboutZ(molec.atoms[i], zRot);
    //translate to original position
    molec.atoms[i] = translateAtom(molec.atoms[i], pivot.x, pivot.y, pivot.z);
  }
  return molec;
}

//returns the index of the first Bond found where Atom1 = atomID
//else returns -1
int hasBond(vector<Bond> Bonds, unsigned long atomId) {
  for (int x = 0; x<Bonds.size(); x++) {
		if (atomId == Bonds[x].atom1) {
      return x;
    }
  }
  return -1;
}

//returns the index of the first Angle found where Atom1 = atomID
//else returns -1
int hasAngle(vector<Angle> Angles, unsigned long atomId) {
  for (int x = 0; x<Angles.size(); x++) {
	  if ((int)atomId == Angles[x].atom1) {
      return x;
    }
  }
  return -1;
}

//returns the index of the first Dihedral found where Atom1 = atomID
//else returns -1
int hasDihedral(vector<Dihedral> Dihedrals, unsigned long atomId) {
  for (int x = 0; x<Dihedrals.size(); x++) {
	  if ((int)atomId == Dihedrals[x].atom1) {
      return x;
    }
  }
  return -1;
}

void setMoleculeVectors(Molecule *molec,  int numBonded, unsigned long lineAtomId, vector<Bond> &bondVector,
    vector<Angle> &angleVector, vector<Dihedral> &dihedralVector) {

	bondVector.clear();
  angleVector.clear();
  dihedralVector.clear();
	int atomId = (int) lineAtomId;

  for (int m = 0; m < numBonded; m++) {
    for (int x = 0; x < molec[m].numOfBonds; x++) {
      if (molec[m].bonds[x].atom1 <= atomId) {
		    bondVector.push_back(molec[m].bonds[x]);
      }
    }

    for (int x = 0; x < molec[m].numOfAngles; x++) {
	    if (molec[m].angles[x].atom1 <= atomId) {
        angleVector.push_back(molec[m].angles[x]);
      }
    }

    for (int x = 0; x < molec[m].numOfDihedrals; x++) {
	    if (molec[m].dihedrals[x].atom1 <= atomId) {
        dihedralVector.push_back(molec[m].dihedrals[x]);
      }
    }
	}
	 //cout << "Sizes Bond Angle Dihedral vectors\n" << bondVector.size()<<" | "<<angleVector.size()<<" | "<<dihedralVector.size()<<endl;
}

void buildMoleculeInSpace(Molecule *molec, int numBonded) {
  vector<Molecule> retVector;
  vector<Atom> atomVector;
  vector<Bond> bondVector;
  vector<Angle> angleVector;
  vector<Dihedral> dihedralVector;
  vector<unsigned long> dummies;
    //run a build on each bonded molecule
	for (int m = 0; m < numBonded; m++) {

    double centX, centY, centZ=0.0; //centerPoint to build molecule around
    Atom lineAtom;
    //run build on each atom in the molecule
    for (int x = 0; x < molec[m].numOfAtoms; x++) {
    	lineAtom = molec[m].atoms[x];
    	//set the vectors with appropiate contents as if in Zmatrix
    	setMoleculeVectors(molec, numBonded, lineAtom.id, bondVector, angleVector, dihedralVector);

    	//use the first atoms position as the centerpoint
    	if (x==0) {
        centX=lineAtom.x;
        centY=lineAtom.y;
        centZ=lineAtom.z;
    	}

    	lineAtom.x=centX;
    	lineAtom.y=centY;
    	lineAtom.z=centZ;

    	int bondFound = hasBond(bondVector,lineAtom.id);

      if (bondFound >= 0) {
    		Bond lineBond =  bondVector[bondFound];
    		// Get other atom in bond
        unsigned long otherID = getOppositeAtom(lineBond, lineAtom.id);

        Atom otherAtom = getAtom(atomVector, otherID);

        // Move newAtom bond distance away from other atom in y direction.
        lineAtom.x = otherAtom.x;
        lineAtom.y = otherAtom.y + lineBond.distance;
        lineAtom.z = otherAtom.z;
    	}

    	int angleFound = hasAngle(angleVector,lineAtom.id);

      if (angleFound >=0 ) {
    		Angle lineAngle = angleVector[angleFound];
    		// Get other atom listed in angle
				unsigned long inval = (unsigned long) -1;
        Atom otherAtom = createAtom(inval, -1, -1, -1);
        unsigned long otherID = getOppositeAtom(lineAngle, lineAtom.id);
        otherAtom = getAtom(atomVector, otherID);

        // Get common atom that lineAtom and otherAtom are bonded to
        //it will be the vertex of the angle.
        unsigned long commonID = getCommonAtom(bondVector, lineAtom.id, otherID);
        Atom commonAtom = getAtom(atomVector, commonID);

        double currentAngle = getAngle(lineAtom, commonAtom, otherAtom);
        double angleChange = lineAngle.value - currentAngle;

        lineAtom = rotateAtomInPlane(lineAtom, commonAtom, otherAtom, angleChange);
    	}

    	int dihedralFound = hasDihedral(dihedralVector, lineAtom.id);

      if (dihedralFound >= 0) {
    		Dihedral lineDihedral = dihedralVector[dihedralFound];
    		//get other atom in the dihedral
        unsigned long otherID = getOppositeAtom(lineDihedral, lineAtom.id);
        Atom otherAtom = getAtom(atomVector, otherID);

        //There are guranteed to be 4 atoms involved in the dihedral
        //because it takes at least 4 atoms to define two non equal
        //planes.

        //get all of the atoms bonded to lineAtom
        vector<unsigned long> bondedToLineAtom = getAllBonds(bondVector, lineAtom.id);
        //get all of the atoms bonded to  otherAtom
        vector<unsigned long> bondedToOtherAtom = getAllBonds(bondVector, otherAtom.id);
        //find bond that bonds together two of the atoms in the intersection
        Bond linkingBond = Bond(-1, -1, -1, false);
        bool foundBond = false;

        // this could possibly be abstracted into its own function and may made not to be and n^3 algorithm. ugh
        for (int i = 0; i < bondedToLineAtom.size(); i++) {
          unsigned long currentToLine = bondedToLineAtom[i];

          for (int j = 0; j < bondedToOtherAtom.size(); j++) {
            unsigned long currentToOther = bondedToOtherAtom[j];

            // cout << "Needs bond between atom " << currentToLine << " and " << currentToOther << endl;
            for (int k = 0; k < bondVector.size(); k++) {
              Bond currentBond = bondVector[k];
              //printf("Current bond: atom1=%d, atom2=%d\n", currentBond.atom1, currentBond.atom2);
              if (getOppositeAtom(currentBond, currentToOther) == currentToLine) {
                linkingBond = currentBond;
                foundBond = true;
              }
              if (foundBond) {
                break;
              }
            }

            if (foundBond) {
              break;
            }
          }

          if (foundBond) {
            break;
          }
        }

        //find atom bonded to the common atom that is not line atom or otherAtom
        if (linkingBond.atom1 == -1 || linkingBond.atom2 == -1) {
          unsigned long commonAtom = getCommonAtom(bondVector, otherAtom.id, lineAtom.id);

          for (int i = 0; i < bondVector.size(); i++) {
            unsigned long opposite = getOppositeAtom(bondVector[i], commonAtom);

            if (opposite != (unsigned long) -1 && opposite != otherAtom.id && opposite != lineAtom.id) {
              linkingBond = bondVector[i];
              break;
            }
          }
        }

    		//plane 1 is lineAtom and atoms in linking bond and will be rotated
        //plane 2 is otherAtom and atoms in linking bond
        //the bond creates the vector about which lineAtom will be rotated.

        Plane rotatePlane = createPlane(lineAtom,
        getAtom(atomVector, linkingBond.atom1),
        getAtom(atomVector, linkingBond.atom2));

        Plane nonMovingPlane = createPlane(otherAtom,
        getAtom(atomVector, linkingBond.atom1),
        getAtom(atomVector, linkingBond.atom2));
        //find the angle between the planes.
        double initialAngle = getAngle(rotatePlane, nonMovingPlane);
        //find the angle needed to rotate.
        double toRotate = initialAngle - lineDihedral.value;
        //rotate lineAtom needed degrees about linkbond.
        //determine which atom in linkingBond is vector head and tail
        Atom vectorHead;
        Atom vectorTail;
        Bond temp = getBond(bondVector, lineAtom.id, linkingBond.atom1);

        if (temp.atom1 == -1 && temp.atom2 == -1) {
          // linkingBond.atom1 is not bonded to line atom and is the tail(start)
          vectorTail = getAtom(atomVector, linkingBond.atom1);
          vectorHead = getAtom(atomVector, linkingBond.atom2);
        } else {
          vectorTail = getAtom(atomVector, linkingBond.atom2);
          vectorHead = getAtom(atomVector, linkingBond.atom1);
        }

        lineAtom = rotateAtomAboutVector(lineAtom, vectorTail, vectorHead, toRotate);
    	}
      atomVector.push_back(lineAtom);
    }//for loop

    //modify/adjust the molecules x,y,z values to the newly calculated ones
    int i=0;
    int x=0;

    if (m>0) {
      x=molec[m-1].numOfAtoms;
    }

    while (i < molec[m].numOfAtoms) {
      molec[m].atoms[i].x = atomVector[x].x;
      molec[m].atoms[i].y = atomVector[x].y;
      molec[m].atoms[i].z = atomVector[x].z;
      i++;
      x++;
    }

  } //for loop for molecule
}

void buildMoleculeXYZ(Molecule *molec, int numBonded) {
  vector<Molecule> retVector;
  vector<Atom> atomVector;
  vector<Bond> bondVector;
  vector<Angle> angleVector;
  vector<Dihedral> dihedralVector;
  vector<unsigned long> dummies;

	Atom lineAtom;

	//run a build on each molecule in zmatrix
	for (int m = 0; m < numBonded; m++) {
	  //run build on atom 1 in the molecule

	  lineAtom = molec[m].atoms[0];
	  //set the vectors with appropiate contents as if in Zmatrix

	  setMoleculeVectors(molec, numBonded, lineAtom.id, bondVector, angleVector, dihedralVector);

	  // First atom at (0,0,0)
	  molec[m].atoms[0].x = 0.0;
	  molec[m].atoms[0].y = 0.0;
	  molec[m].atoms[0].z = 0.0;

	  //run build on atom 2 in the molecule
	  lineAtom = molec[m].atoms[1];
	  setMoleculeVectors(molec, numBonded, lineAtom.id, bondVector, angleVector, dihedralVector);

	  int bondFound = hasBond(bondVector,lineAtom.id);

    if (bondFound >= 0) {
		  Bond lineBond =  bondVector[bondFound];
      // Second atom on x-axis
		  molec[m].atoms[1].x = lineBond.distance;
		  molec[m].atoms[1].y = 0.0;
		  molec[m].atoms[1].z = 0.0;
	  } else {
	    break;
    }

		//run build on atom 3 in the molecule
		lineAtom = molec[m].atoms[2];
		setMoleculeVectors(molec, numBonded, lineAtom.id, bondVector, angleVector, dihedralVector);

		int angleFound = hasAngle(angleVector,lineAtom.id);

    if (angleFound >=0 ) {
			int bondFound = hasBond(bondVector,lineAtom.id);
			Bond lineBond =  bondVector[bondFound];
			Angle lineAngle = angleVector[angleFound];
			unsigned long BondedTo = getOppositeAtom(lineBond, lineAtom.id);
			double thetaRadians = degreesToRadians(lineAngle.value);

			//adjust BondedTo index
			for (int i = 0; i < m; i++)
				BondedTo -= molec[i].numOfAtoms;
			// Third atom in XY plane
			if (BondedTo == 1) {
				molec[m].atoms[2].x = lineBond.distance * cos(thetaRadians);
      } else if (BondedTo == 2) {
				molec[m].atoms[2].x = molec[m].atoms[1].x - lineBond.distance * cos(thetaRadians);
			}
			molec[m].atoms[2].y = lineBond.distance * sin(thetaRadians);
			molec[m].atoms[2].z = 0.0;

    } else {
			break;
    }

		//check if atom 4 exists
		if (molec[m].numOfAtoms < 4) {
		  break;
    }

		//run build on atoms 4 and above in the molecule
		double xbs, ybs, zbs, sinval;
		double ia[3], ib[3], ic[3];
		double vecangdih[3],vecangbnd[3];
		double vx[3], vy[3], vz[3];

		for (int N = 3; N < molec[m].numOfAtoms; N++) {
			lineAtom = molec[m].atoms[N];
			setMoleculeVectors(molec, numBonded, lineAtom.id, bondVector, angleVector, dihedralVector);
			int bondFound = hasBond(bondVector,lineAtom.id);
			Bond lineBond =  bondVector[bondFound];
			int angleFound = hasAngle(angleVector,lineAtom.id);
			Angle lineAngle = angleVector[angleFound];
			double thetaRadians = degreesToRadians(lineAngle.value);
			// Dihedral angle is defined counterclockwise
			int dihedralFound = hasDihedral(dihedralVector,lineAtom.id);
			Dihedral lineDihedral = dihedralVector[dihedralFound];
			double phiRadians = degreesToRadians(lineDihedral.value);

			// x|y|z in the local system
			// xbs = bndlgth * sin(ang) * cos(dih)
			// ybs = bndlgth * sin(ang) * sin(dih)
			// zbs = -bndlgth * cos(ang)

			sinval = sin(thetaRadians);
			xbs = lineBond.distance * sinval * cos(phiRadians);
			ybs = lineBond.distance * sinval * sin(phiRadians);
			zbs = -lineBond.distance * cos(thetaRadians);

	    // adjustAtomsIDs(molec, m);
			// Determine transformation (direction cosine) matrix
			unsigned long BondedTo = getOppositeAtom(lineBond, lineAtom.id);
			unsigned long AngleWith = getOppositeAtom(lineAngle, lineAtom.id);
			unsigned long DihedralWith = getOppositeAtom(lineDihedral, lineAtom.id);

			//adjust indexing for BondedTo, AngleWith, and DihedralWith
			for (int i = 0; i < m; i++) {
				BondedTo -= molec[i].numOfAtoms;
				AngleWith -= molec[i].numOfAtoms;
				DihedralWith -= molec[i].numOfAtoms;
			}

			ia[0] = molec[m].atoms[DihedralWith - 1].x;
			ia[1] = molec[m].atoms[DihedralWith - 1].y;
			ia[2] = molec[m].atoms[DihedralWith - 1].z;
			ib[0] = molec[m].atoms[AngleWith - 1].x;
			ib[1] = molec[m].atoms[AngleWith - 1].y;
			ib[2] = molec[m].atoms[AngleWith - 1].z;
			ic[0] = molec[m].atoms[BondedTo - 1].x;
			ic[1] = molec[m].atoms[BondedTo - 1].y;
			ic[2] = molec[m].atoms[BondedTo - 1].z;

			for (int i = 0; i < 3; i++) {
				vecangdih[i] = ia[i] - ic[i];
				vecangbnd[i] = ib[i] - ic[i];
    	}

			cross(vecangdih,vecangbnd,vy);
			cross(vecangbnd,vy,vx);
			cross(vx,vy,vz);

			//Map the coordinates in this basis set to our cartesian coordinates
			//     x = xbs*vx(0) + ybs*vy(0) + zbs*vz(0)
			//     y = xbs*vx(1) + ybs*vy(1) + zbs*vz(1)
			//     z = xbs*vx(2) + ybs*vy(2) + zbs*vz(2)
			//
			//These coordinates are based at the origin - they need to be based from
			//the coordinates of the bond atom e.g.
			// x += xbnd     y+= ybnd   z+= zbnd

			molec[m].atoms[N].x = vx[0]*xbs + vy[0]*ybs + vz[0]*zbs + ic[0];
			molec[m].atoms[N].y = vx[1]*xbs + vy[1]*ybs + vz[1]*zbs + ic[1];
			molec[m].atoms[N].z = vx[2]*xbs + vy[2]*ybs + vz[2]*zbs + ic[2];

		}
	  //revertAtomIDs(molec, m);
	} //End for molecules loop
}

void adjustAtomIDs(Molecule* molec, int m) {
  if (m > 0) {
		std::cout << "Adjusting ID..." << std::endl;

    for (int i = 0; i < molec[m].numOfAtoms; i++) {
			std::cout << "ID of atom " << i << ": " << molec[m].atoms[i].id << std::endl;
		}

		for (int i = 0; i < molec[m].numOfAtoms; i++) {
			for (int j = 0; j < m; j++) {
		    molec[m].atoms[i].id -= molec[j].numOfAtoms;
			}
		}

  	for (int i = 0; i < molec[m].numOfAtoms; i++) {
			std::cout << "ID of atom " << i << ": " << molec[m].atoms[i].id << std::endl;
    }

  }
}

void revertAtomIDs(Molecule* molec, int m) {
  if (m > 0) {
		std::cout << "Reverting ID bigger than 3 back to normal..." << std::endl;

    for (int i = 0; i < molec[m].numOfAtoms; i++) {
			std::cout << "ID of atom " << i << ": " << molec[m].atoms[i].id << std::endl;
    }

		for (int i = 0; i < molec[m].numOfAtoms; i++) {
			for (int j = 0; j < m; j++) {
			  molec[m].atoms[i].id += molec[j].numOfAtoms;
			}
		}

    for (int i = 0; i < molec[m].numOfAtoms; i++) {
		  std::cout << "ID of atom " << i << ": " << molec[m].atoms[i].id << std::endl;
    }

	} //End for dihedral atoms loop
}

void cross(double *A, double *B, double *C) {
	// C = A X B; C returned normalized
	const double tiny = 1.00e-36;
	double x, y, z;

	x = A[1]*B[2] - B[1]*A[2];
	y = B[0]*A[2] - A[0]*B[2];
	z = A[0]*B[1] - B[0]*A[1];

	//Avoid division by zero
	double length = std::max(tiny, ((x*x) + (y*y) + (z*z)));
	double lengthSQRT = sqrt(length);

	C[0] = x / lengthSQRT;
	C[1] = y / lengthSQRT;
	C[2] = z / lengthSQRT;
}

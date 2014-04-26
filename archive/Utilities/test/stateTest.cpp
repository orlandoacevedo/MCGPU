#include "stateTest.h"

void testGetHopFromLine()
{
    cout << "Testing getHopFromLine" << endl;
    string line = "3 4 2";
    string line2 = "6 3 6";

    Hop h1 = getHopFromLine(line);
    Hop h2 = getHopFromLine(line2);

    assert(h1.atom1 == 3);
    assert(h1.atom2 == 4);
    assert(h1.hop == 2);

    assert(h2.atom1 == 6);
    assert(h2.atom2 == 3);
    assert(h2.hop == 6);

    cout << "Testing getHopFromLine Complete\n" << endl;

}

void testGetDihedralFromLine()
{
    cout << "Testing getDihedralFromLine" << endl;
    string line = "3 5 3.456 0";
    string line2 = "89 21 45.21 1";

    Dihedral d = getDihedralFromLine(line);
    Dihedral d2 = getDihedralFromLine(line2);

    assert(d.atom1 == 3);
    assert(d.atom2 == 5);
    assert(d.value == 3.456);
    assert(d.variable == false);

    assert(d2.atom1 = 89);
    assert(d2.atom2 = 21);
    assert(d2.value = 3.456);
    assert(d2.variable = true);

    cout << "Testing getDihedralFromLine Complete\n" << endl;
}

void testGetAtomFromLine()
{
    cout << "Testing getAtomFromLine" << endl;
    
    //id x y z sigma epsilon charg
    string line1 = "2, 4.234, 56.34, 1.23, 4.12, 9.76 -2.4";
    string line2 = "4, 7.89, 2.12345, 34.5, 65.7, 7.96 5.46";

    Atom a1 = getAtomFromLine(line1);
    Atom a2 = getAtomFromLine(line2);

    assert(a1.id == 2);
    assert(a1.x == 4.234);
    assert(a1.y == 56.34);
    assert(a1.z == 1.23);
    assert(a1.sigma == 4.12);
    assert(a1.epsilon == 9.76);
    assert(a1.charge == -2.4);

    assert(a2.id == 4);
    assert(a2.x == 7.89);
    assert(a2.y == 2.12345);
    assert(a2.z == 34.5);
    assert(a2.sigma == 65.7);
    assert(a2.epsilon == 7.96);
    assert(a2.charge == 5.46);

    cout << "Testing getAtomFromLine Complete\n" << endl;
}

void testGetBondFromLine()
{
    cout << "Testing getBondFromLine" << endl;
    string line = "3 5 3.456 0";
    string line2 = "89 21 45.21 1";

    Bond d = getBondFromLine(line);
    Bond d2 = getBondFromLine(line2);

    assert(d.atom1 == 3);
    assert(d.atom2 == 5);
    assert(d.distance == 3.456);
    assert(d.variable == false);

    assert(d2.atom1 = 89);
    assert(d2.atom2 = 21);
    assert(d2.distance = 3.456);
    assert(d2.variable = true);

    cout << "Testing getBondFromLine Complete\n" << endl;
}

void testGetAngleFromLine()
{
    cout << "Testing getAngleFromLine" << endl;
    string line = "3 5 3.456 0";
    string line2 = "89 21 45.21 1";

    Angle d = getAngleFromLine(line);
    Angle d2 = getAngleFromLine(line2);

    assert(d.atom1 == 3);
    assert(d.atom2 == 5);
    assert(d.value == 3.456);
    assert(d.variable == false);

    assert(d2.atom1 = 89);
    assert(d2.atom2 = 21);
    assert(d2.value = 3.456);
    assert(d2.variable = true);

    cout << "Testing getAngleFromLine Complete\n" << endl;
}

void testGetEnvironmentFromLine()
{
    cout << "Testing getEnvironmentFromLine" << endl;
    string line1 = "10.0 11.2 32.34 90 150.4";
    string line2 = "11.4 34.1 12.54 45 160.9";

    Environment e1 = getEnvironmentFromLine(line1);
    Environment e2 = getEnvironmentFromLine(line2);

    assert(e1.x == 10.0);
    assert(e1.y == 11.2);
    assert(e1.z == 32.34);
    assert(e1.numOfAtoms == 90);
    assert(e1.temperature == 150.4);

    assert(e2.x == 11.4);
    assert(e2.y == 34.1);
    assert(e2.z == 12.54);
    assert(e2.numOfAtoms == 45);
    assert(e2.temperature == 160.9);

    cout << "Testing getEnvironmentFromLine Complete\n" << endl;
}

void testWriteOutReadInState()
{
    cout << "Testing readInEnvironment" << endl;
    Atom atom1 = createAtom(0, 1.23f, 2.3f, 4.3f, 5.6f, 4.34f);
    Atom atom2 = createAtom(1, 3.22f, 4.2f, 6.5f, 8.6f, 6.36f);
    Atom atom3 = createAtom(2, 2.21f, 7.4f, 2.8f, 4.2f, 2.22f);
    Atom atomArray[3] = {atom1, atom2, atom3};

    Bond b1 = createBond(2,1,32.4, true);
    Bond b2 = createBond(0,1,3.1, false);
    Bond bondArray[2] = {b1, b2};

    Angle a1 = createAngle(1,0, 56.0, true);
    Angle a2 = createAngle(2,1,23.2, false);
    Angle angleArray[2] = {a1, a2};

    Dihedral d1 = createDihedral(1,2,4.35,false);
    Dihedral d2 = createDihedral(0,2,5.34,true);
    Dihedral dihedralArray[2] = {d1, d2};
    
    Hop h1 = createHop(1,2,3);
    Hop h2 = createHop(4,2,5);
    Hop hopArray[2] = {h1, h2};

    Environment enviro = createEnvironment(10.f, 10.f, 10.f, .5, 273, 2, 9.0, 15.0);

    Molecule molec1 = createMolecule(1, atomArray, angleArray, bondArray, dihedralArray, hopArray,
            3,2,2,2,2);

    Molecule molec2 = createMolecule(2, atomArray, &a1, &b1, dihedralArray, hopArray,
            3,1,1,2, 2);

    Molecule molecArray[2] = {molec1, molec2};
    string filename = "STATE_TEST.state";
    printState(&enviro, molecArray, 2, filename);

    printf("State written out\n");

    /*==================
      STATE HAS BEEN WRITTEN OUT
    ==================*/

    Environment enviroTest = readInEnvironment(filename);
    assert(enviro.x == enviroTest.x);
    assert(enviro.y == enviroTest.y);
    assert(enviro.z == enviroTest.z);
    assert(enviro.numOfAtoms == enviroTest.numOfAtoms);

    vector<Molecule> molecVector = readInMolecules(filename); 

    
    for(int j = 0; j < molecVector.size(); j++)
    { 
        double floatThreshold = .00001;
        //print atom
        Molecule molec = molecVector[j];
        for(int i = 0; i < molec.numOfAtoms; i++)
        {
            Atom a = molec.atoms[i];
            assert( fabs(a.z - atomArray[i].z) < floatThreshold);
            assert( fabs(a.x - atomArray[i].x) < floatThreshold);
            assert( fabs(a.y - atomArray[i].y) < floatThreshold);
            assert( fabs(a.sigma - atomArray[i].sigma) < floatThreshold);
            assert( fabs(a.epsilon - atomArray[i].epsilon) < floatThreshold);
            assert(a.id == atomArray[i].id);
        }
        //print dihedrals
        for(int i = 0; i < molec.numOfDihedrals; i++)
        {
            Dihedral d = molec.dihedrals[i];
            assert(d.atom1 == dihedralArray[i].atom1);
            assert(d.atom2 == dihedralArray[i].atom2);
            assert(fabs(d.value - dihedralArray[i].value) < floatThreshold);
        }
        //print bonds
        for(int i = 0; i < molec.numOfBonds; i++)
        {
            
            Bond d = molec.bonds[i];
            assert(d.atom1 == bondArray[i].atom1);
            assert(d.atom2 == bondArray[i].atom2);
            assert(fabs(d.distance - bondArray[i].distance) < floatThreshold);

        }
        // test angles 
        for(int i = 0; i < molec.numOfAngles; i++)
        {
            
            Angle d = molec.angles[i];
            assert(d.atom1 == angleArray[i].atom1);
            assert(d.atom2 == angleArray[i].atom2);
            assert(fabs(d.value - angleArray[i].value) < floatThreshold);
        }
        //test hops
        for(int i = 0; i < molec.numOfHops; i++)
        {
            Hop h = molec.hops[i];
            assert(h.atom1 == hopArray[i].atom1);
            assert(h.atom2 == hopArray[i].atom2);
            assert(h.hop == hopArray[i].hop);
        }
    }

    cout << "Testing readInEnvironment Complete\n" << endl;
}



void runStateTests()
{
    testGetHopFromLine();
    testGetDihedralFromLine();
    testGetAtomFromLine();
    testGetBondFromLine();
    testGetAngleFromLine();
    testGetEnvironmentFromLine();
    testWriteOutReadInState();
}


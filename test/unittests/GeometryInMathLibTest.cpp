
#include "Metropolis/Utilities/MathLibrary.h"
#include "gtest/gtest.h"
#include <vector>

#define PRECISION .0001

// Descr: evident
// Implementation details: See gtest/samples for GTest syntax and usage
TEST(GeometryTest, GetNormal)
{

	Atom atom1, atom2, atom3;
	atom1.x = 3.2;
	atom1.y = 4.6;
	atom1.z = 0.2;
	atom2.x = 7.1;
	atom2.y = 1.4;
	atom2.z = 10.0;
	atom3.x = 0.0;
	atom3.y = 0.0;
	atom3.z = 0.0;

    Plane testPlane = createPlane(atom1, atom2, atom3);

    Point expected = createPoint(45.72, -30.58, -28.18);

    Point test = getNormal(testPlane);

    test.x = ((double) ((int) (test.x * 100))) / 100.0;
    test.y = ((double) ((int) (test.y * 100))) / 100.0;
    test.z = ((double) ((int) (test.z * 100))) / 100.0;

    EXPECT_NEAR(expected.x, test.x, .01);
    EXPECT_NEAR(expected.y, test.y, .01);
    EXPECT_NEAR(expected.z, test.z, .01);

}

// Descr: evident
TEST(GeometryTest, GetAngleBetweenPlanes)
{
    srand(time(NULL));

    int numberOfTests = 100;

    for (int i = 0; i < numberOfTests; i++)
    {
        Plane a, b;

        a.atom1.x = ((double) rand() / RAND_MAX) * 10;
        a.atom1.y = ((double) rand() / RAND_MAX) * 10;
        a.atom1.z = ((double) rand() / RAND_MAX) * 10;

        a.atom2.x = ((double) rand() / RAND_MAX) * 10;
        a.atom2.y = ((double) rand() / RAND_MAX) * 10;
        a.atom2.z = ((double) rand() / RAND_MAX) * 10;

        a.atom3.x = ((double) rand() / RAND_MAX) * 10;
        a.atom3.y = ((double) rand() / RAND_MAX) * 10;
        a.atom3.z = ((double) rand() / RAND_MAX) * 10;

        b.atom1.x = ((double) rand() / RAND_MAX) * 10;
        b.atom1.y = ((double) rand() / RAND_MAX) * 10;
        b.atom1.z = ((double) rand() / RAND_MAX) * 10;

        b.atom2.x = ((double) rand() / RAND_MAX) * 10;
        b.atom2.y = ((double) rand() / RAND_MAX) * 10;
        b.atom2.z = ((double) rand() / RAND_MAX) * 10;

        b.atom3.x = ((double) rand() / RAND_MAX) * 10;
        b.atom3.y = ((double) rand() / RAND_MAX) * 10;
        b.atom3.z = ((double) rand() / RAND_MAX) * 10;

        Point aNormal = getNormal(a);
        Point bNormal = getNormal(b);

        double numerator = aNormal.x * bNormal.x + aNormal.y * bNormal.y + aNormal.z * bNormal.z;
        double aMag = sqrt(aNormal.x * aNormal.x + aNormal.y * aNormal.y + aNormal.z * aNormal.z);
        double bMag = sqrt(bNormal.x * bNormal.x + bNormal.y * bNormal.y + bNormal.z * bNormal.z);
        double denominator = aMag * bMag;

        double thetaR = acos(numerator / denominator);
        double expectedTheta = radiansToDegrees(thetaR);

        double testTheta = getAngle(a, b);

		ASSERT_LT( (testTheta - expectedTheta) / expectedTheta, PRECISION);
    }
}


// Descr: evident
TEST(GeometryTest, GetBond)
{
    vector<Bond> bonds;


    Bond bond1 = Bond(1, 2, 3.0, false);
    Bond bond2 = Bond(3, 4, 3.0, false);
    Bond bond3 = Bond(2, 3, 3.0, false);

    bonds.push_back(bond1);
    EXPECT_EQ( 1, bonds[0].atom1);
    bonds.push_back(bond2);
    bonds.push_back(bond3);

    Bond testBond1 = getBond(bonds, 1, 2);
    Bond testBond2 = getBond(bonds, 3, 4);
    Bond testBond3 = getBond(bonds, 2, 3);
    Bond testBond4 = getBond(bonds, 1, 4);

    EXPECT_EQ( bond1.atom1, testBond1.atom1);
    EXPECT_EQ( -1 , testBond4.atom1);

	EXPECT_TRUE( bond1.atom1 == testBond1.atom1 && bond1.atom2 == testBond1.atom2 && bond1.distance == testBond1.distance );
	EXPECT_TRUE( bond2.atom1 == testBond2.atom1 && bond2.atom2 == testBond2.atom2 && bond2.distance == testBond2.distance );
	EXPECT_TRUE( bond3.atom1 == testBond3.atom1 && bond3.atom2 == testBond3.atom2 && bond3.distance == testBond3.distance );
	EXPECT_TRUE( testBond4.atom1 == -1 && testBond4.atom2 == -1 && testBond4.distance == -1.0 );

}

// Descr: evident
TEST(GeometryTest, GetAllBonds)
{
	vector<Bond> bonds;

    Bond bond1 = Bond(1, 2, 3.0, false);
    Bond bond2 = Bond(3, 4, 3.0, false);
    Bond bond3 = Bond(2, 3, 3.0, false);

    bonds.push_back(bond1);
    bonds.push_back(bond2);
    bonds.push_back(bond3);

    vector<unsigned long> testBonds1 = getAllBonds(bonds, 1);
    vector<unsigned long> testBonds2 = getAllBonds(bonds, 2);
    vector<unsigned long> testBonds3 = getAllBonds(bonds, 3);
    vector<unsigned long> testBonds4 = getAllBonds(bonds, 4);
    vector<unsigned long> testBonds5 = getAllBonds(bonds, 5);

	EXPECT_TRUE(testBonds1[0] == 2 && testBonds1.size() == 1);
	EXPECT_TRUE(testBonds2[0] == 1 && testBonds2[1] == 3 && testBonds2.size() == 2);
	EXPECT_TRUE(testBonds3[0] == 4 && testBonds3[1] == 2 && testBonds3.size() == 2);
	EXPECT_TRUE(testBonds4[0] == 3 && testBonds4.size() == 1);
	EXPECT_TRUE(testBonds5.size() == 0);
}

// Descr: evident
TEST(GeometryTest, GetIntersection)
{
    vector<unsigned long> section1, section2;

    section1.push_back(1);
    section1.push_back(2);
    section1.push_back(3);
    section1.push_back(4);

    section2.push_back(2);
    section2.push_back(4);
    section2.push_back(6);
    section2.push_back(8);

    vector<unsigned long> intersection = getIntersection(section1, section2);

	EXPECT_TRUE(intersection[0] == 2 && intersection[1] == 4 && intersection.size() == 2);
}

// Descr: tests the isMember function
TEST(GeometryTest, IsMember)
{

    vector<unsigned long> section;

    section.push_back(1);
    section.push_back(2);
    section.push_back(3);
    section.push_back(4);

	EXPECT_TRUE(isMember(section, 1) );
	EXPECT_TRUE(isMember(section, 2) );
	EXPECT_TRUE(isMember(section, 3) );
	EXPECT_TRUE(isMember(section, 4) );
	EXPECT_TRUE(!isMember(section, 5));

}

// Descr: test degrees to radians and radians to degrees
// Implementation details: See gtest/samples for GTest syntax and usage
TEST(GeometryTest, D2RandR2D)
{
    //test converting Degrees to Radians
    // I know this 'should' be EXPECT_DOUBLE_EQ, but even though the values were good to the 5th or 6th decimal place, it was failing.
    // So, there you go.
    EXPECT_FLOAT_EQ(degreesToRadians(1),0.0174532925 );
    EXPECT_FLOAT_EQ(degreesToRadians(254),4.4331363 );
    EXPECT_FLOAT_EQ(degreesToRadians(360),6.283185307 );
    EXPECT_FLOAT_EQ(degreesToRadians(15),0.261799388 );
    EXPECT_FLOAT_EQ(radiansToDegrees(3.14),179.908747671 );
    EXPECT_FLOAT_EQ(radiansToDegrees(0.174532925),10 );
    EXPECT_FLOAT_EQ(radiansToDegrees(1.745329252),100 );
    EXPECT_FLOAT_EQ(radiansToDegrees(.1234567),7.073547863 );
    EXPECT_FLOAT_EQ(radiansToDegrees(6.195918845),355 );
}


// Descr: evident
TEST(GeometryTest, GetOppositeAtom)
{
    Bond testBond = Bond(1, 3, 5.5, true);
    Angle testAngle = Angle(2, 6, 30, false);
    Dihedral testDihed = Dihedral(3, 7, 180, true);

    //test with bonds
    EXPECT_EQ(3, getOppositeAtom(testBond, 1));
    EXPECT_EQ(1, getOppositeAtom(testBond, 3));
    EXPECT_EQ(-1, getOppositeAtom(testBond, 4));

    EXPECT_EQ(6, getOppositeAtom(testAngle, 2));
    EXPECT_EQ(2, getOppositeAtom(testAngle, 6));
    EXPECT_EQ(-1, getOppositeAtom(testAngle, 5));

    EXPECT_EQ(7, getOppositeAtom(testDihed, 3));
    EXPECT_EQ(3, getOppositeAtom(testDihed, 7));
    EXPECT_EQ(-1, getOppositeAtom(testDihed, 6));

    // Second suite
    testBond = Bond(11, 12, 3.5, true);
    testAngle = Angle(1, 22, 30, false);
    testDihed = Dihedral(5, 9, 180, true);


    EXPECT_EQ(12, getOppositeAtom(testBond, 11));
    EXPECT_EQ(11, getOppositeAtom(testBond, 12));
    EXPECT_EQ(-1, getOppositeAtom(testBond, 13));

    EXPECT_EQ(22, getOppositeAtom(testAngle, 1));
    EXPECT_EQ(1, getOppositeAtom(testAngle, 22));
    EXPECT_EQ(-1, getOppositeAtom(testAngle, 15));

    EXPECT_EQ(9, getOppositeAtom(testDihed, 5));
    EXPECT_EQ(5, getOppositeAtom(testDihed, 9));
    unsigned long inval = (unsigned long) -1;
    EXPECT_EQ(-1, getOppositeAtom(testDihed, inval));
}


// Descr: evident
TEST(GeometryTest, GetCommonAtom)
{
	vector<Bond> bondVect(7);
    bondVect[0] =Bond(2,1,0.5,false);
    bondVect[1] =Bond(3,2,0.5,false);
    bondVect[2] =Bond(4,1,1.336532,true);
    bondVect[3] =Bond(5,1,1.8119,true);
    bondVect[4] =Bond(6,5,1.090187,true);
    bondVect[5] =Bond(7,5,1.090135,true);
    bondVect[6] =Bond(8,5,1.090135,true);

	EXPECT_EQ(1, getCommonAtom(bondVect,4,5) );
	EXPECT_EQ(2, getCommonAtom(bondVect,1,3) );
	EXPECT_EQ(5, getCommonAtom(bondVect,6,8) );
	EXPECT_EQ(1, getCommonAtom(bondVect,2,5) );
	EXPECT_EQ(-1, getCommonAtom(bondVect,3,5) );
	EXPECT_EQ(-1, getCommonAtom(bondVect,8,2) );

}

// Descr: evident
TEST(GeometryTest, GetDistance)
{
    Atom atom1= Atom(1,5,6,7);
    Atom atom2= Atom(2,10,11,12);
    EXPECT_NEAR(8.66025, getDistance(atom1,atom2) , .001);

    atom1= Atom(1,8,12,21);
    atom2= Atom(2,4,5,10);
    EXPECT_FLOAT_EQ( 13.638181, getDistance(atom1,atom2) );

    atom1= Atom(1,45,2,22);
    atom2= Atom(2,37,22,18);
    EXPECT_FLOAT_EQ( 21.9089023002, getDistance(atom1,atom2) );
}


// Descr: evident
TEST(GeometryTest, GetAngle)
{
    Atom atom1= Atom(1,5,6,7);
    Atom atom2= Atom(2,10,11,12);
    Atom atom3= Atom(3,14,22,9);
    EXPECT_NEAR( 124.986, getAngle(atom1,atom2,atom3) , .001);

    atom1= Atom(1,15,23,8);
    atom2= Atom(2,5,3,12);
    atom3= Atom(3,9,18,7);
    EXPECT_NEAR( 13.6609, getAngle(atom1,atom2,atom3) , .001 );
}

// Descr: evident
TEST(GeometryTest, TranslateAtom)
{
	srand(time(NULL));
    Atom testAtom;
    for (int i = 0; i < 2; i++)
    {
        double random = ((double) rand() / RAND_MAX) * 15;
        switch(i)
        {
            case 0:
                testAtom.x = random;
                break;
            case 1:
                testAtom.y = random;
                break;
            case 2:
                testAtom.z = random;
                break;
            default:
                break;
        }
    }

    double oldX = testAtom.x;
    double oldY = testAtom.y;
    double oldZ = testAtom.z;

    double translateX = ((double) rand() / RAND_MAX) * 15;
    double translateY = ((double) rand() / RAND_MAX) * 15;
    double translateZ = ((double) rand() / RAND_MAX) * 15;

    testAtom = translateAtom(testAtom, translateX, translateY, translateZ);

	EXPECT_NEAR( oldX + translateX , testAtom.x, .1);
	EXPECT_NEAR( oldY + translateY , testAtom.y, .1);
	EXPECT_NEAR( oldZ + translateZ , testAtom.z, .1);
}

// Descr: evident
TEST(GeometryTest, RotateAboutX)
{
    srand(time(NULL));
    Atom testAtom;
    for (int i = 0; i < 2; i++)
    {
        double random = ((double) rand() / RAND_MAX) * 15;
        switch(i)
        {
            case 0:
                testAtom.x = random;
                break;
            case 1:
                testAtom.y = random;
                break;
            case 2:
                testAtom.z = random;
                break;
            default:
                break;
        }
    }

    double dtr = PI / 180.0;

    double oldY = testAtom.y;
    double oldZ = testAtom.z;

    double theta = ((double) rand() / RAND_MAX) * 180;

    testAtom = rotateAboutX(testAtom, theta);

    double expectedY = cos(theta * dtr) * oldY + sin(theta * dtr) * oldZ;
    double expectedZ = cos(theta * dtr) * oldZ - sin(theta * dtr) * oldY;

	EXPECT_NEAR( expectedY, testAtom.y, .001  );
	EXPECT_NEAR( expectedZ, testAtom.z, .001  );
}

// Descr: evident
TEST(GeometryTest, RotateAboutY)
{
   srand(time(NULL));
    Atom testAtom;
    for (int i = 0; i < 2; i++)
    {
        double random = ((double) rand() / RAND_MAX) * 15;
        switch(i)
        {
            case 0:
                testAtom.x = random;
                break;
            case 1:
                testAtom.y = random;
                break;
            case 2:
                testAtom.z = random;
                break;
            default:
                break;
        }
    }

    double dtr = PI / 180.0;

    double oldX = testAtom.x;
    double oldZ = testAtom.z;

    double theta = ((double) rand() / RAND_MAX) * 180;

    testAtom = rotateAboutY(testAtom, theta);

    double expectedX = cos(theta * dtr) * oldX - sin(theta * dtr) * oldZ;
    double expectedZ = cos(theta * dtr) * oldZ + sin(theta * dtr) * oldX;


	EXPECT_NEAR( expectedX, testAtom.x, .001  );
	EXPECT_NEAR( expectedZ, testAtom.z, .001  );
}

// Descr: evident
TEST(GeometryTest, RotateAboutZ)
{
    srand(time(NULL));
    Atom testAtom;
    for (int i = 0; i < 2; i++)
    {
        double random = ((double) rand() / RAND_MAX) * 15;
        switch(i)
        {
            case 0:
                testAtom.x = random;
                break;
            case 1:
                testAtom.y = random;
                break;
            case 2:
                testAtom.z = random;
                break;
            default:
                break;
        }
    }

    double dtr = PI / 180.0;

    double oldX = testAtom.x;
    double oldY = testAtom.y;

    double theta = ((double) rand() / RAND_MAX) * 180;

    testAtom = rotateAboutZ(testAtom, theta);

    double expectedX = cos(theta * dtr) * oldX + sin(theta * dtr) * oldY;
    double expectedY = cos(theta * dtr) * oldY - sin(theta * dtr) * oldX;


	EXPECT_NEAR( expectedX, testAtom.x, .001  );
	EXPECT_NEAR( expectedY, testAtom.y, .001  );
}

// Descr: evident
TEST(GeometryTest, RotateAboutVector)
{
    double rotation = 90;
    unsigned long inval = (unsigned long) -1;
    Atom vertex = createAtom(inval, 0, 0, 0);
    Atom head = createAtom(inval, 1, 0, 0);
    Atom toRotate = createAtom(inval, 0, 1, 0);

    Atom rotated = rotateAtomAboutVector(toRotate, vertex, head, rotation);
    EXPECT_LT( fabs(rotated.y - 0.0), .01);
    EXPECT_LT( fabs(rotated.x - 0.0), .01);
    EXPECT_LT( fabs(rotated.z - 1.0), .01);

    rotation = 45;

    rotated = rotateAtomAboutVector(toRotate, vertex, head, rotation);
    EXPECT_LT( fabs(rotated.y - sqrt(.5)), .01);
    EXPECT_LT( fabs(rotated.x - 0.0), .01);
    EXPECT_LT( fabs(rotated.z - sqrt(.5)), .01);


    rotation = 90;
    //test rotating about atom with 0 intial angle
    toRotate = createAtom(inval, 2, 0, 0);
    rotated = rotateAtomAboutVector(toRotate, vertex, head, rotation);
    EXPECT_LT( fabs(rotated.y - 0), .01);
    EXPECT_LT( fabs(rotated.z - 0), .01);
    EXPECT_LT( fabs(rotated.x - 2), .01);

}

// Descr: evident
TEST(GeometryTest, RotateInPlane)
{
    double rotation = 90;
    unsigned long inval = (unsigned long) -1;
    Atom vertex = Atom(inval, 0, .5, 0);
    Atom head = Atom(inval, 0, 0, 0);
    Atom toRotate = Atom(inval, 0, 1, 0);

    Atom rotated = rotateAtomInPlane(toRotate, vertex, head, rotation);

    EXPECT_LT( fabs(rotated.y - .5), .01);

    rotation = 45;

    vertex = Atom(inval, 0, 0, 0);
    head = Atom(inval, 0, 0, -1);
    toRotate = Atom(inval, 0, 1, 0);

    rotated = rotateAtomInPlane(toRotate, vertex, head, rotation);
    EXPECT_LT( fabs(rotated.y - sqrt(.5)), .01);
    EXPECT_LT( fabs(rotated.z - sqrt(.5)), .01);
    EXPECT_LT( fabs(rotated.x - 0.0), .01);
}


#include "Metropolis/Utilities/MathLibrary.h"
#include "gtest/gtest.h"


// Tests factorial of negative numbers.

// Tests
// Descr: evident
// Implementation details: See gtest/samples for GTest syntax and usage
TEST(GeometryTest, GetNormal) {
  
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
    
    EXPECT_EQ(expected.x, test.x);
    EXPECT_EQ(expected.y, test.y);
    EXPECT_EQ(expected.z, test.z);
    
}

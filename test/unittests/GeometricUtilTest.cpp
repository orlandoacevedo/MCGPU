#include "Utilities/geometricUtil.h"
#include "gtest/gtest.h"

TEST(GeometricUtilTest, IsMemberInvalidInputs)
{
	vector<unsigned long> list;
	EXPECT_FALSE(isMember(list,1));

	list.push_back(1);
	list.push_back(3);
	EXPECT_FALSE(isMember(list,2));
	EXPECT_FALSE(isMember(list,0));
}

TEST(GeometricUtilTest, IsMemberValidInputs)
{
	vector<unsigned long> list;
	list.push_back(1);
	list.push_back(3);
	list.push_back(1);

	EXPECT_TRUE(isMember(list,3));
	EXPECT_TRUE(isMember(list,1));
}

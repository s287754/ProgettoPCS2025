/*#pragma once 

#include <gtest/gtest.h>

#include "Geometry.hpp"

namespace GeometryLibrary {
	
	TEST(TestGeometry,TestCell0Ds) 
	{
		Cell0Ds v{0,1.5,-2.3,4.0};
		EXPECT_EQ(v.id,0.0);
		EXPECT_DOUBLE_EQ(v.x,1.5);
		EXPECT_DOUBLE_EQ(v.y,-2.3);
		EXPECT_DOUBLE_EQ(v.z,4.0);
	}
	
	TEST(TestGeometry,TestCell1Ds)
	{
		Cell1Ds e{1,0,3};
		EXPECT_EQ(e.id,1);
		EXPECT_EQ(e.origin,0);
		EXPECT_EQ(e.end,3);
	}
	
}*/
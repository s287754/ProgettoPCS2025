#include "Geometry.hpp"
#include "Utils.hpp"

#include <gtest/gtest.h>

using namespace PolyhedralLibrary;

TEST(GeometryTest,CoordinatesBefore)
{
	Polyhedral p;
	ASSERT_EQ(p.NumCell0Ds,0);
	ASSERT_EQ(p.NumCell1Ds,0);
	ASSERT_EQ(p.NumCell2Ds,0);
	
	ASSERT_EQ(p.Cell0DsCoordinates.size(),0);
	ASSERT_EQ(p.Cell1DsVertices.size(),0);
	ASSERT_TRUE(p.Cell2DsVertices.empty());
	ASSERT_TRUE(p.Cell2DsEdges.empty());
}
TEST(GeometryTest,CoordinatesAfter)
{
	Polyhedral p1;
	Polyhedral p2;
	Polyhedral p3;
	
	createTetrahedron(p1);
	ASSERT_EQ(p1.NumCell0Ds,4);
	ASSERT_EQ(p1.NumCell1Ds,6);
	ASSERT_EQ(p1.NumCell2Ds,4);
	
	ASSERT_EQ(p1.Cell0DsCoordinates.rows(),p1.NumCell0Ds);
	ASSERT_EQ(p1.Cell1DsVertices.rows(),p1.NumCell1Ds);
	ASSERT_EQ(p1.Cell2DsVertices.size(),p1.NumCell2Ds);
	
	createCube(p2);
	ASSERT_EQ(p2.NumCell0Ds,8);
	ASSERT_EQ(p2.NumCell1Ds,12);
	ASSERT_EQ(p2.NumCell2Ds,6);
	
	ASSERT_EQ(p2.Cell0DsCoordinates.rows(),p2.NumCell0Ds);
	ASSERT_EQ(p2.Cell1DsVertices.rows(),p2.NumCell1Ds);
	ASSERT_EQ(p2.Cell2DsVertices.size(),p2.NumCell2Ds);
	
	createIcosahedron(p3);
	ASSERT_EQ(p3.NumCell0Ds,12);
	ASSERT_EQ(p3.NumCell1Ds,30);
	ASSERT_EQ(p3.NumCell2Ds,20);
	
	ASSERT_EQ(p3.Cell0DsCoordinates.rows(),p3.NumCell0Ds);
	ASSERT_EQ(p3.Cell1DsVertices.rows(),p3.NumCell1Ds);
	ASSERT_EQ(p3.Cell2DsVertices.size(),p3.NumCell2Ds);
}
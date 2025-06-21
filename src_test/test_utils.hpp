#include "Utils.hpp"

#include <gtest/gtest.h>

using namespace PolyhedralLibrary;

TEST(UtilsTest,DualT)
{
	Polyhedral polygon;
	createTetrahedron(polygon);
	
	Polyhedral dual = DualPolygon(polygon);
	
	EXPECT_EQ(dual.NumCell0Ds, polygon.NumCell2Ds);
	EXPECT_EQ(polygon.NumCell0Ds, dual.NumCell2Ds);
	EXPECT_EQ(dual.NumCell3Ds,1);
	
}
TEST(UtilsTest,DualCoordinatesT)
{
	Polyhedral polygon;
	createTetrahedron(polygon);
	Polyhedral dual = DualPolygon(polygon);
	
	const auto& face = polygon.Cell2DsVertices[0];
	
	Vector3d barycenter = Vector3d :: Zero();
	
	for (auto v : face)
		barycenter += polygon.Cell0DsCoordinates.row(v).transpose();
	
	barycenter /= face.size();
	
	Vector3d dualVertices = dual.Cell0DsCoordinates.row(0).transpose();
	
	EXPECT_EQ((dualVertices - barycenter).norm(),0.0);
		
}
TEST(UtilsTest,DualC)
{
	Polyhedral polygon;
	createCube(polygon);
	
	Polyhedral dual = DualPolygon(polygon);
	
	EXPECT_EQ(dual.NumCell0Ds, polygon.NumCell2Ds);
	EXPECT_EQ(polygon.NumCell0Ds, dual.NumCell2Ds);
	EXPECT_EQ(dual.NumCell3Ds,1);
	
}
TEST(UtilsTest,DualCoordinatesC)
{
	Polyhedral polygon;
	createCube(polygon);
	Polyhedral dual = DualPolygon(polygon);
	
	const auto& face = polygon.Cell2DsVertices[0];
	
	Vector3d barycenter = Vector3d :: Zero();
	
	for (auto v : face)
		barycenter += polygon.Cell0DsCoordinates.row(v).transpose();
	
	barycenter /= face.size();
	
	Vector3d dualVertices = dual.Cell0DsCoordinates.row(0).transpose();
	
	EXPECT_EQ((dualVertices - barycenter).norm(),0.0);
		
}
TEST(UtilsTest,DualOfDual)
{
	Polyhedral polygon;
	createCube(polygon);
	Polyhedral dual_1 = DualPolygon(polygon);
	Polyhedral dual_2 = DualPolygon(dual_1);
	
	EXPECT_EQ(dual_2.NumCell0Ds, polygon.NumCell0Ds);
	EXPECT_EQ(dual_2.NumCell1Ds, polygon.NumCell1Ds);
	EXPECT_EQ(dual_2.NumCell2Ds, polygon.NumCell2Ds);
}
TEST(UtilsTest,DualI)
{
	Polyhedral polygon;
	createIcosahedron(polygon);
	
	Polyhedral dual = DualPolygon(polygon);
	
	EXPECT_EQ(dual.NumCell0Ds, polygon.NumCell2Ds);
	EXPECT_EQ(polygon.NumCell0Ds, dual.NumCell2Ds);
	EXPECT_EQ(dual.NumCell3Ds,1);
	
}
TEST(UtilsTest,DualCoordinatesI)
{
	Polyhedral polygon;
	createIcosahedron(polygon);
	Polyhedral dual = DualPolygon(polygon);
	
	const auto& face = polygon.Cell2DsVertices[0];
	
	Vector3d barycenter = Vector3d :: Zero();
	
	for (auto v : face)
		barycenter += polygon.Cell0DsCoordinates.row(v).transpose();
	
	barycenter /= face.size();
	
	Vector3d dualVertices = dual.Cell0DsCoordinates.row(0).transpose();
	
	EXPECT_EQ((dualVertices - barycenter).norm(),0.0);
		
}

#include "Geometry.hpp"

#include <gtest/gtest.h>

using namespace std;
using namespace PolyhedralLibrary;

TEST(PolyhedralTest,Inizialization)
{
	Polyhedral polygon;
	
	EXPECT_EQ(polygon.NumCell0Ds,0);
	EXPECT_TRUE(polygon.Cell0DsId.empty());
	EXPECT_TRUE(polygon.MarkerCell0Ds.empty());
	EXPECT_EQ(polygon.Cell0DsCoordinates.size(),0);
	
	EXPECT_EQ(polygon.NumCell1Ds,0);
	EXPECT_TRUE(polygon.Cell1DsId.empty());
	EXPECT_TRUE(polygon.MarkerCell1Ds.empty());
	EXPECT_EQ(polygon.Cell1DsVertices.size(),0);
	
	EXPECT_EQ(polygon.NumCell2Ds,0);
	EXPECT_EQ(polygon.NumVertices,0);
	EXPECT_EQ(polygon.NumEdges,0);
	EXPECT_TRUE(polygon.Cell2DsId.empty());
	EXPECT_EQ(polygon.Cell2DsVertices.size(),0);
	EXPECT_EQ(polygon.Cell2DsEdges.size(),0);
	
	EXPECT_EQ(polygon.NumCell3Ds,0);
	EXPECT_TRUE(polygon.Cell3DsId.empty());
	EXPECT_EQ(polygon.Cell3DsVertices.size(),0);
	EXPECT_EQ(polygon.Cell3DsEdges.size(),0);
	EXPECT_EQ(polygon.Cell3DsFaces.size(),0);
}
TEST(PolyhedraTest,AddData)
{
	Polyhedral polygon;
	
	polygon.NumCell0Ds = 2;
	polygon.Cell0DsId = {0,1};
	
	polygon.Cell0DsCoordinates = MatrixXd (2,4);
	polygon.Cell0DsCoordinates << 0, 0.0, 0.0, 0.0,
								  1, 1.0, 1.0, 1.0;
  
    EXPECT_EQ(polygon.Cell0DsCoordinates.rows(),2);
	EXPECT_EQ(polygon.Cell0DsCoordinates(1,1),1.0);
	
	polygon.NumCell1Ds = 2;
	polygon.Cell1DsId = {0,1};
	
	polygon.Cell1DsVertices = MatrixXi (2,3);
	polygon.Cell1DsVertices << 0, 0,1,
							   1, 0,2;
							   
	EXPECT_EQ(polygon.Cell1DsVertices.rows(),2);
	EXPECT_EQ(polygon.Cell1DsVertices(1,1),0);
	
	polygon.NumCell2Ds = 2;
	polygon.NumVertices = 2;
	polygon.NumEdges = 2;
	polygon.Cell2DsId = {0,1};
	
	polygon.Cell2DsVertices = vector<vector<unsigned int>>{};
	polygon.Cell2DsEdges = vector<vector<unsigned int>>{};
	polygon.Cell2DsVertices = {{0,1,2},
							   {0,1,3}};
	polygon.Cell2DsEdges = {{0,3,2},
						    {0,3,4}};
	
	EXPECT_EQ(polygon.Cell2DsVertices[1][1],1);
	EXPECT_EQ(polygon.Cell2DsEdges[1][1],3);
	
	polygon.NumCell3Ds = 2;
	polygon.Cell3DsId = {0,1};
	
	polygon.Cell3DsVertices = vector<vector<unsigned int>>{};
	polygon.Cell3DsEdges = vector<vector<unsigned int>>{};
	polygon.Cell3DsFaces = vector<vector<unsigned int>>{};
	polygon.Cell3DsVertices = {{0,1,2},
							   {0,1,3}};
	polygon.Cell3DsEdges = {{0,3,2},
							{0,3,4}};
	polygon.Cell3DsFaces = {{0,1,2,3},
							{0,2,3,4}};
	
	EXPECT_EQ(polygon.Cell3DsVertices[1][1],1);
	EXPECT_EQ(polygon.Cell3DsEdges[1][1],3);
	EXPECT_EQ(polygon.Cell3DsFaces[1][1],2);
	
}
	
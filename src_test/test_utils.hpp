#include "Utils.hpp"

#include <gtest/gtest.h>

using namespace PolyhedralLibrary;

TEST(orientTest, RightOrientation)
{
	Eigen::MatrixXi M(6,2) ;
	M << 0,1,0,2,0,3,1,2,1,3,2,3;
	vector<vector<unsigned int>> v ={{0,3,1}} ;
	unsigned int n = v.size();
	ASSERT_TRUE(checkOrient(n,M,v));
}
TEST(orientTest, WrongOrientation)
{
	Eigen :: MatrixXi M(6,2);
	M << 0,1,0,2,0,3,1,2,1,3,2,3;
	vector<vector<unsigned int>> v ={{0,5,1}};
	unsigned int n = v.size();
	ASSERT_FALSE(checkOrient(n,M,v));
}
TEST(DualTest,RightCoordinates)
{
	Eigen::MatrixXd C(3,3);
	Eigen::Vector3d barycenter;
	Eigen::Vector3d v;
	C << +1,+1,+1,-1,-1,+1,-1,+1,-1;
	barycenter = (C.colwise().sum())/3;
	v = {-1.0/3,+1.0/3,+1.0/3};
	ASSERT_EQ(barycenter,v);
}
TEST(DualTest,WrongCoordinates)
{
	Eigen::MatrixXd C(3,3);
	Eigen::Vector3d barycenter;
	Eigen::Vector3d v;
	Eigen::Vector3d b;
	Eigen::Vector3d z;
	C << +1,+1,+1,-1,-1,+1,-1,+1,-1;
	barycenter = (C.colwise().sum())/3;
	v = {-1.0/3,+1.0/3,+1.0/3};
	z = {+1.0,+1.0,+1.0};
	b = v+z;
	ASSERT_NE(barycenter,b);
}
TEST(GeodesicTest,CorrectGeneration)
{
	Polyhedral polygon;
	polygon.NumCell0Ds = 3;
	polygon.NumCell1Ds = 1;
	polygon.NumCell2Ds = 1;
	polygon.Cell2DsVertices.push_back({0,1,2});
	
	polygon.Cell0DsCoordinates.resize(3,3);
	polygon.Cell0DsCoordinates << 1,0,0,0,1,0,0,0,1;
	polygon.Cell1DsVertices.resize(3,2);
	polygon.Cell1DsVertices << 0,1,0,2,1,2;
	
	GeodesicPolyhedron(3,3,0,2,polygon);
	
	ASSERT_EQ(polygon.Cell0DsCoordinates.rows(),6);
	ASSERT_EQ(polygon.Cell1DsVertices.rows(),9);
	ASSERT_EQ(polygon.Cell2DsVertices.size(),4);
}
TEST(VectorCompareTest,RightOrdering)
{
	Eigen::Vector3d v1 = {1.0,0.0,0.0};
	Eigen::Vector3d v2 = {1.0,1.0,0.0};
	Eigen::Vector3d v3 = {1.0,0.0,0.0};
	
	Vector3dCompare comp;
	
	ASSERT_TRUE(comp(v1,v2));
	ASSERT_FALSE(comp(v2,v1));
	ASSERT_FALSE(comp(v1,v3));
}
TEST(shortestPathTest,RightPath)
{
	Polyhedral polygon;
	polygon.Cell0DsCoordinates.resize(3,3);
	polygon.Cell0DsCoordinates << 0,0,0,+1,+1,+1,-2,-2,-2;
	polygon.Cell1DsVertices.resize(3,2);
	polygon.Cell1DsVertices << 0,1,0,2,1,2;
	unsigned int start = 0;
	unsigned int end = 2;
	vector<unsigned int> path;
	vector<unsigned int> crossedEdges;
	shortestPath(polygon,start,end,path,crossedEdges);
	
	vector<unsigned int> expected_path = {0,1,2} ;
	vector<unsigned int> expected_crossedEdges = {0,1,2};
	ASSERT_EQ(path,expected_path);
	ASSERT_EQ(crossedEdges,expected_crossedEdges);
}
TEST(shortestPathTest,WrongPath)
{
	Polyhedral polygon;
	polygon.Cell0DsCoordinates.resize(3,3);
	polygon.Cell0DsCoordinates << 0,0,0,+1,+1,+1,-2,-2,-2;
	polygon.Cell1DsVertices.resize(3,2);
	polygon.Cell1DsVertices << 0,1,0,2,1,2;
	unsigned int start = 0;
	unsigned int end = 2;
	vector<unsigned int> path;
	vector<unsigned int> crossedEdges;
	shortestPath(polygon,start,end,path,crossedEdges);
	
	vector<unsigned int> expected_path = {0,2,1} ;
	vector<unsigned int> expected_crossedEdges = {1,0,2};
	ASSERT_NE(path,expected_path);
	ASSERT_NE(crossedEdges,expected_crossedEdges);
}
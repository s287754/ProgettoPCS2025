#include <iostream>
#include <vector>
#include <Eigen/Eigen>
#include <cmath>
#include "Utils.hpp"
#include "Geometry.hpp"

using namespace std;
using namespace Eigen;

namespace PolyhedralLibrary {
bool createTetrahedron (Polyhedral& polygon)
{
	polygon.NumCell0Ds = 4;
	
	for (unsigned int i = 0; i < polygon.NumCell0Ds; i++)
		polygon.Cell0DsId.push_back(i);
	
	const double x = +1.0/sqrt(3);
	polygon.Cell0DsCoordinates.resize(polygon.NumCell0Ds,3);
	polygon.Cell0DsCoordinates << +x,+x,+x, /*A: id 0*/ 
								  -x,-x,+x, /*B: id 1*/ 
								  -x,+x,-x, /*C: id 2*/ 
								  +x,-x,-x; /*D: id 3*/ 
	
	for(unsigned int i = 0; i < polygon.NumCell0Ds; ++i)
		polygon.MarkerCell0Ds[1].push_back(i);
    
    polygon.NumCell1Ds = 6;
	
	for (unsigned int i = 0; i < polygon.NumCell1Ds; i++)
		polygon.Cell1DsId.push_back(i);
	
	polygon.Cell1DsVertices.resize(polygon.NumCell1Ds,2);
	polygon.Cell1DsVertices << 0,1, /*AB: id 0*/
							   0,2, /*AC: id 1*/
							   0,3, /*AD: id 2*/
							   1,2, /*BC: id 3*/
							   1,3, /*BD: id 4*/
							   2,3; /*CD: id 5*/
							   
  
   
   for(unsigned int i = 0; i < polygon.NumCell1Ds; ++i)
		polygon.MarkerCell1Ds[1].push_back(i);

	polygon.NumCell2Ds = 4;
	
	for (unsigned int i = 0; i < polygon.NumCell2Ds; i++)
		polygon.Cell2DsId.push_back(i);
	
	polygon.NumVertices = 3;
	
	polygon.Cell2DsVertices = {{0,1,2}, /*ABC*/ 
							   {2,0,3}, /*CAD*/							 
							   {3,2,1},
							   {1,3,0}}; /*DCB*/
	
	polygon.NumEdges = 3;

	polygon.Cell2DsEdges = {{0,3,1}, /* ABC */
						   {0,4,2}, /* CAD */
						   {1,5,2}, /* DCB */
						   {3,5,4}}; /* BDA */
	
	checkOrient(polygon); 
	
	polygon.NumCell3Ds = 1; /* Un unico tetraedro*/
	polygon.Cell3DsId.push_back(0); /* Un solo ID */
	polygon.Cell3DsVertices.push_back(polygon.Cell0DsId); /* ID dei vertici che lo compongono */
	polygon.Cell3DsEdges.push_back(polygon.Cell1DsId); /* ID dei lati che lo compongono */
	polygon.Cell3DsFaces.push_back(polygon.Cell2DsId); /* ID delle facce che lo compongono */
	return true;
}

bool createCube (Polyhedral& polygon)  /*Creo una funzione che popoli la mia struct Polyhedral*/
{
	polygon.NumCell0Ds = 8;
	
	for (unsigned int i = 0; i < polygon.NumCell0Ds; i++)
		polygon.Cell0DsId.push_back(i);
	const double x = +1.0/sqrt(3);
	polygon.Cell0DsCoordinates.resize(polygon.NumCell0Ds,3);
	polygon.Cell0DsCoordinates << +x,+x,+x, /*A*/
							      -x,+x,+x, /*B*/
							      -x,-x,+x, /*C*/
							      +x,-x,+x, /*D*/
							      +x,+x,-x, /*E*/
							      +x,-x,-x, /*F*/
							      -x,-x,-x, /*G*/
							      -x,+x,-x; /*H*/
							   
   polygon.NumCell1Ds = 12;
   
   for (unsigned int i = 0; i < polygon.NumCell1Ds; i++)
	   polygon.Cell1DsId.push_back(i);
   
   polygon.Cell1DsVertices.resize(polygon.NumCell1Ds,2);/*Controlla bene il resize e paragonalo al reserve*/
   polygon.Cell1DsVertices << 0,1, /* AB */
							  1,2, /* BC */
							  2,3, /* CD */
							  3,0, /* DA */
							  4,5, /* EF */
							  5,6, /* FG */
							  6,7, /* GH */
							  7,4, /* HE */
						      0,4, /* AE */
							  1,5, /* BF */
							  2,6, /* CG */
							  3,7; /* DH */
						   
   polygon.NumCell2Ds = 6;
   
   for (unsigned int i = 0; i < polygon.NumCell2Ds; i++)
	   polygon.Cell2DsId.push_back(i);
   
   polygon.NumVertices = 4;
   polygon.NumEdges = 4;
   
   polygon.Cell2DsVertices = {{0,1,2,3}, /* ABCD */
							  {4,5,6,7}, /* EFGH */
							  {0,1,5,4}, /* ABFE */
							  {1,2,6,5}, /* BCGF */
							  {2,3,7,6}, /* CDEH */
							  {3,0,4,7}}; /* DAEH */
							  
   polygon.Cell2DsEdges = {{0,1,2,3}, /* AB-BC-CD-DA */ 
						   {4,5,6,7}, /* EF-FG-GH-HE */
						   {0,9,4,8}, /* AB-BF-EF-AE */
						   {1,10,5,9}, /* BC-CG-FG-BF */
						   {2,11,6,10}, /* CD-DH-GH-CG */
						   {3,8,7,11}}; /* DA-AE-HE-DH */
						   
   checkOrient(polygon);
   
   polygon.NumCell3Ds = 1; /* Un unico cubo*/
   polygon.Cell3DsId.push_back(0); /* Un solo ID */
   polygon.Cell3DsVertices.push_back(polygon.Cell0DsId); /* ID dei vertici che lo compongono */
   polygon.Cell3DsEdges.push_back(polygon.Cell1DsId); /* ID dei lati che lo compongono */
   polygon.Cell3DsFaces.push_back(polygon.Cell2DsId); /* ID delle facce che lo compongono */
							   
   return true;
}

bool createIcosahedron (Polyhedral& polygon)
{
	polygon.NumCell0Ds = 12 ;
	for (unsigned int i = 0; i < polygon.NumCell0Ds; i++)
		polygon.Cell0DsId.push_back(i);
	
	polygon.Cell0DsCoordinates.resize(polygon.NumCell0Ds,3);
	const double phi = (1.0 + sqrt(5.0)) / 2.0;
	const double k = 1.0/(sqrt(phi*phi + 1));
	polygon.Cell0DsCoordinates << -1.0 * k, +phi * k,  0,
								  +1.0 * k, +phi * k,  0,
								  -1.0 * k, -phi * k,  0,
								  +1.0 * k, -phi * k,  0,
								   0,  -1.0 * k, +phi * k,
								   0,  +1.0 * k, +phi * k,
								   0,  -1.0 * k, -phi * k,
								   0,  +1.0 * k, -phi * k,
								  +phi * k,  0, -1.0 * k,
								  +phi * k,  0, +1.0 * k,
								  -phi * k,  0, -1.0 * k,
								  -phi * k,  0, +1.0 * k;
								  
								  
								  
	polygon.NumCell1Ds = 30;
	
	for(unsigned int i = 0; i < polygon.NumCell1Ds; i++)
		polygon.Cell1DsId.push_back(i);
	
	polygon.Cell1DsVertices.resize(polygon.NumCell1Ds,2);
	polygon.Cell1DsVertices <<   0, 1, //0
								 0, 5, //1
								 0, 7, //2
								 0, 10, //3
								 0, 11, //4
								 1, 5, //5
								 1, 7, //6
								 1, 8, //7
								 1, 9, //8
								 2, 3, //9
								 2, 4,//10
								 2, 6, //11
								 2, 10, //12
								 2, 11, //13
								 3, 4, //14
								 3, 6, //15
								 3, 8, //16
								 3, 9, //17
								 4, 5, //18
								 4, 9, //19
								 4, 11, //20
								 5, 9, //21
								 5, 11, //22
								 6, 7, //23
								 6, 8, //24
								 6, 10, //25
								 7, 8, //26
								 7, 10, //27
								 8, 9, //28
								 10, 11; //29 
						
	polygon.NumCell2Ds = 20;
	
	for (unsigned int i = 0; i < polygon.NumCell2Ds; i++)
		polygon.Cell2DsId.push_back(i);
	
	polygon.NumVertices = 3;
	
	polygon.Cell2DsVertices = {{0, 11, 5},
							   {0, 5, 1},
							   {0, 1, 7},
							   {0, 7, 10},
							   {0, 10, 11},
							   {1, 5, 9},
							   {5, 11, 4},
							   {11, 10, 2},
							   {10, 7, 6},
							   {7, 1, 8},
							   {3, 9, 4},
							   {3, 4, 2},
							   {3, 2, 6},
						       {3, 6, 8},
							   {3, 8, 9},
							   {4, 9, 5},
							   {2, 4, 11},
							   {6, 2, 10},
							   {8, 6, 7},
							   {9, 8, 1}};				  

	polygon.NumEdges = 3;
	
	polygon.Cell2DsEdges = {{4, 22, 1},
							{1, 5, 0},
							{0, 6, 2},
							{2, 27, 3},
							{3, 29, 4},
							{5, 21, 8},
							{22, 20, 18},
							{29, 12, 13},
							{27, 23, 25},
							{6, 7, 26},
							{17, 19, 14},
							{14, 10, 9},
							{9, 11, 15},
							{15, 24, 16},
							{16, 28, 17},
							{19, 21, 18},
							{10, 20, 13},
							{11, 12, 25},
							{24, 23, 26},
							{28, 7, 8}};
							
	checkOrient (polygon);
	
	polygon.NumCell3Ds = 1; /* Un unico icosaedro*/
	polygon.Cell3DsId.push_back(0); /* Un solo ID */
	polygon.Cell3DsVertices.push_back(polygon.Cell0DsId); /* ID dei vertici che lo compongono */
	polygon.Cell3DsEdges.push_back(polygon.Cell1DsId); /* ID dei lati che lo compongono */
	polygon.Cell3DsFaces.push_back(polygon.Cell2DsId); /* ID delle facce che lo compongono */
	
	return true;
}
}
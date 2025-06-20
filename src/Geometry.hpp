#pragma once

#include <iostream>
#include <vector>
#include <array>
#include <Eigen/Eigen>

using namespace std;
using namespace Eigen;

namespace PolyhedralLibrary {
	
	struct Polyhedral
	{
		unsigned int NumCell0Ds = 0;
		vector<unsigned int> Cell0DsId = {}; /*ID vertici*/
		map<unsigned int,vector<unsigned int>> MarkerCell0Ds = {};
		MatrixXd Cell0DsCoordinates = {}; /* Coordinate vertici n:=numero di righe=numero di vertici + 1 di intestazione; m:=numero di colonne=3(x,y,z) + 1(ID) */
		
		unsigned int NumCell1Ds = 0;
		vector<unsigned int> Cell1DsId = {};
		map<unsigned int,vector<unsigned int>> MarkerCell1Ds = {};
		MatrixXi Cell1DsVertices = {}; /*ID vertici di origine e di fine n=numero di lati+1; m=3(2 IDOrigin/IDEnd + 1)*/

		unsigned int NumCell2Ds = 0;
		unsigned int NumVertices = 0;
		unsigned int NumEdges = 0;
		vector<unsigned int> Cell2DsId = {};
		vector<vector<unsigned int>> Cell2DsVertices = {}; /* ID Vertici n:=numeri dei "primi" vector = numero di facce che compongono il poliedro; m:=numero del "secondo" vector=numero di vertici della faccia */
		vector<vector<unsigned int>> Cell2DsEdges = {}; /*ID Lati n=numero di facce che compongono il poliedro; m=numero di lati della faccia; */
		
		unsigned int NumCell3Ds = 0;
		vector<unsigned int> Cell3DsId = {};
		vector<vector<unsigned int>> Cell3DsVertices = {}; /*ID Vertici*/
		vector<vector<unsigned int>> Cell3DsEdges = {}; /*ID Lati */
		vector<vector<unsigned int>> Cell3DsFaces = {}; /*ID Facce */

	};
}
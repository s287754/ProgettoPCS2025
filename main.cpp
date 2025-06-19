#include <iostream>
#include <array>
#include "Utils.hpp"
#include "Geometry.hpp"
#include "UCDUtilities.hpp"

using namespace std;
using namespace PolyhedralLibrary;

int main () {
	
	cout << "Inserisci una quadrupla di valori:\n" << endl;
	
	array <unsigned int,4> input;
	cin >> input[0] >> input[1] >> input[2] >> input[3];
	
	 
	Polyhedral tetrahedron;
	createTetrahedron(tetrahedron);
	cout << "Num facce prima di GeodesicPolyhedron: " << tetrahedron.NumCell2Ds << endl;
	DualPolygon(tetrahedron);
	GeodesicPolyhedron(input[0],input[1],input[2],input[3],tetrahedron);
	
	Polyhedral cube;
	createCube(cube);
	DualPolygon(cube);
	
	Polyhedral icosahedron;
	createIcosahedron(icosahedron);
	DualPolygon(icosahedron);
	
	exportToUCD(tetrahedron, "tetrahedron");

	return 0;
}
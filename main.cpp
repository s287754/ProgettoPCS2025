#include <iostream>
#include <array>
#include "Utils.hpp"
#include "Geometry.hpp"
#include "UCDUtilities.hpp"

using namespace std;
using namespace PolyhedralLibrary;

int main () {
	
	cout << "Inserisci al minimo una quadrupla o al massimo una sestupla di valori\nGli ultimi due saranno gli id dei vertici tra cui" 
				"calcolare il cammino minimo\nPremi invio per terminare:\n" << endl;
	
	vector <unsigned int> input;
	string riga;
	
	getline (cin, riga);
	istringstream ss(riga);
	
	unsigned int num;
	
	while (ss >> num)
	{
		input.push_back(num);
		if (input.size() == 6)
			break;
	}
	
	Polyhedral tetrahedron;
	createTetrahedron(tetrahedron);
	
	Polyhedral d_tetrahedron = DualPolygon(tetrahedron); 
	
	Polyhedral cube;
	createCube(cube);
	
	Polyhedral octahedron = DualPolygon(cube);
	
	Polyhedral icosahedron;
	createIcosahedron(icosahedron);
	
	Polyhedral dodecahedron = DualPolygon(icosahedron);
	
	unsigned int p = input[0];
	unsigned int q = input[1];
	unsigned int b = input[2];
	unsigned int c = input[3];
	unsigned int start = input[4];
	unsigned int end = input[5];
	
	if (p < 3 || q < 3)
	{
		cout << "Errore, i valori di p e q devono essere entrambi maggiori di 3!" << endl;
		return false;
	}
	
	if (b != 0 && c != 0)
	{
		cout << "Per ottendere il poliedro geodetico di ordine I è necessario che o b o c sia uguale a zero" << endl;
		return false;
	}
	
	Polyhedral poly;

	string name;
	
	if (q == 3)
	{
		poly = tetrahedron;
		name = " Tetraedro";
	}
	else if (q == 4)
	{
		poly = octahedron;
		name = "l' Ottaedro";
	}
	else if (q == 5)
	{
		poly = icosahedron;
		name = "l'Icosaedro";
	}
	else 
	{
		cout << "Valore di q non valido" << endl;
		return 1;
	}

	if (!GeodesicPolyhedron(p, q, b, c, poly))
	{
		cout << "Errore nell'input, p è diverso da tre ma uguale a " << p << endl;
		return false;
	}
	
	cout << "Ho creato il poliedro geodetico del" << name << endl ;
	
	Polyhedral dualPoly = DualPolygon(poly);
	
	for (unsigned int i = 0; i < dualPoly.Cell0DsCoordinates.rows(); ++i) {
    dualPoly.Cell0DsCoordinates.row(i).normalize();
	}
	cout << "Ho creato anche il duale del poliedro geodetico del" << name << endl;
	vector<unsigned int> path;
	vector<unsigned int> crossedEdges;

	if (shortestPath(poly,start,end,path,crossedEdges))
		exportToUCD (poly, name,path, crossedEdges);
	
	vector<unsigned int> emptyPath;
	vector<unsigned int> emptyEdges;
	exportToUCD(dualPoly,"duale_" + name,emptyPath,emptyEdges);
	
	return 0;
}
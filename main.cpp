#include <iostream>
#include <array>
#include "Utils.hpp"
#include "Geometry.hpp"
#include "UCDUtilities.hpp"

using namespace std;
using namespace PolyhedralLibrary;

int main () {
	
	cout << "Inserisci al minimo una quadrupla o al massimo una sestupla di valori\n" << endl;
	cout << "Gli ultimi due saranno gli id dei vertici tra cui calcolare il cammino minimo\nPremi invio per terminare:\n" << endl;
	
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
	
	DualPolygon(tetrahedron);
	
	Polyhedral cube;
	createCube(cube);
	
	DualPolygon(cube);
	
	Polyhedral icosahedron;
	createIcosahedron(icosahedron);
	
	DualPolygon(icosahedron);
	
	unsigned int p = input[0];
	unsigned int q = input[1];
	unsigned int b = input[2];
	unsigned int c = input[3];
	unsigned int start = input[4];
	unsigned int end = input[5];

	Polyhedral poly;

	string name;
	
	if (q == 3)
	{
		poly = tetrahedron;
		name = "Tetraedro";
	}
	else if (q == 4)
	{
		Polyhedral octahedron = DualPolygon(cube);
		poly = octahedron;
		name = "Ottaedro";
	}
	else if (q == 5)
	{
		poly = icosahedron;
		name = "Icosaedro";
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
	
	cout << "Ho creato il poliedro geodetico del " << name << endl ;
	
	vector<unsigned int> path;
	if (shortestPath(poly, start, end, path)) {
		cout << "Cammino minimo trovato:\n";
		for (int idx : path)
			cout << idx << " ";
		cout << std::endl;
	} else {
		cout << "Nessun cammino trovato tra " << start << " e " << end << std::endl;
	}
	
	exportToUCD (poly, name);

	return 0;
}
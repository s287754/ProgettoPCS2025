#include <iostream>
#include <fstream>
#include <iomanip>
#include <random>
#include <string>
#include <vector>
#include "DataGenerator.hpp"
#include "Geometry.hpp"

using namespace std;
using namespace GeometryLibrary;

void dataGeneratorCell0Ds (unsigned int n, const string& filename)
{
	
	ofstream file(filename);
	if (file.fail())
	
	{
		cerr <<"Error"<< endl;
	}

	file << "ID" << ";" << "x" << ";" << "y" << ";" << "z" << endl;
	file << scientific << setprecision(6);

	for (unsigned int i = 0; i < n; i++)
	{
		double x = static_cast<double>(rand())/RAND_MAX;
		double y = static_cast<double>(rand())/RAND_MAX;
		double z = static_cast<double>(rand())/RAND_MAX;
	
		file << i << ";" << x << ";" << y << ";" << z << endl;
	}
	
	file.close();
}

void dataGeneratorCell1Ds (unsigned int n, const string& filename, const string& filenameInput)
{
	ifstream fileIn(filenameInput);
	ofstream file(filename);
	if (file.fail() && fileIn.fail())
	{
		cerr << "Error" << endl;
		return;
	}
	
	string line;
	getline (fileIn,line); //Elimino la prima riga
	
	int* id = new int[n];
	while (!fileIn.eof())
	{
		for (unsigned int i = 0; i < n; i++)
		{
			getline (fileIn,line);
			istringstream convert(line);
			convert >> id[i];
		}
	}
	fileIn.close();
	
	file << "ID" << ";" << "Origin" << ";" << "End" << endl;
	for (unsigned int i = 0; i < n; i++)
	{
		unsigned int idOrigin = rand () % 20;
		unsigned int idEnd = rand () % 20;
		
		file << i << "\t" << idOrigin << "\t" << idEnd << endl;
	}
	
	delete [] id;
	file.close();
}
		

	

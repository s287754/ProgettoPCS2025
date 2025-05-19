#include <iostream>
#include <fstream>
#include <string>
#include "DataGenerator.hpp"
#include "Geometry.hpp"

using namespace std;
using namespace GeometryLibrary;

int main ()
{
	cout << "RandMax: " << RAND_MAX <<endl;
	
	int r = rand();
    double t = static_cast<double>(r) / RAND_MAX;
    cout << "rand: " << r << ", t: " << t << std::endl;
	dataGeneratorCell0Ds(20, "Cell0Ds.txt");
	dataGeneratorCell1Ds(20, "Cell1Ds.txt", "Cell0Ds.txt");
	return 0;
}
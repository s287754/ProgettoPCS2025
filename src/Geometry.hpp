#pragma once

#include <string>

using namespace std;

namespace GeometryLibrary 
{
	struct Cell0Ds {
		unsigned int id;
		double x,y,z;
	};
	
	struct Cell1Ds {
		unsigned int id;
		unsigned int origin;
		unsigned int end;
	};
};
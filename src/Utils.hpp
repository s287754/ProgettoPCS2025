#pragma once

#include <iostream>
#include "Geometry.hpp"

using namespace std;

namespace PolyhedralLibrary
{
	bool createCube (Polyhedral& polygon);
	bool createTetrahedron (Polyhedral& polygon);
	bool createIcosahedron (Polyhedral& polygon);
	
	bool checkOrient (const Polyhedral& polygon);
	Polyhedral DualPolygon (const Polyhedral& polygon);
	
	bool GeodesicPolyhedron(const unsigned int& p,const unsigned int& q,const unsigned int& b,const unsigned int& c, Polyhedral& poly);
	void exportToUCD(const Polyhedral &polygon, const string& basename);
	
	struct Vector3dCompare {
		bool operator()(const Vector3d& v1, const Vector3d& v2) const 
		{
			if (v1.x() != v2.x()) return v1.x() < v2.x();
			if (v1.y() != v2.y()) return v1.y() < v2.y();
			return v1.z() < v2.z();
		}
	};
}
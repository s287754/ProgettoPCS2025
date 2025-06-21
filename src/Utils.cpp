#include <iostream>
#include <fstream>
#include <set>
#include <algorithm>
#include <unordered_set>
#include <vector>
#include <queue>

#include "Utils.hpp"
#include "Geometry.hpp"
#include "UCDUtilities.hpp"

using namespace std;
using namespace Eigen;
namespace PolyhedralLibrary {

bool checkOrient(const Polyhedral& polygon)
{
	
	unsigned int F = polygon.NumCell2Ds;
	
	for (size_t i = 0; i < F; i ++)
	{
		
		const auto& faceEdges = polygon.Cell2DsEdges[i];
		size_t E = faceEdges.size();
		
		for (size_t e = 0; e < E; e++)
		{
			unsigned int currentEdgeId = faceEdges[e];
			unsigned int nextEdgeId = faceEdges [(e+1)%E];
			unsigned int currentEdgeStart = polygon.Cell1DsVertices(currentEdgeId,0);
			unsigned int currentEdgeEnd = polygon.Cell1DsVertices(currentEdgeId,1);
			unsigned int nextEdgeStart = polygon.Cell1DsVertices(nextEdgeId,0);
			unsigned int nextEdgeEnd = polygon.Cell1DsVertices(nextEdgeId,1);
			
			bool valid = false;
			
			if (currentEdgeEnd == nextEdgeStart)
			{
				valid = true;
			}
			else if (currentEdgeStart == nextEdgeEnd)
			{
				valid = true;
			}
			else if (currentEdgeEnd == nextEdgeEnd)
			{
				valid = true;
			}
			else if (currentEdgeStart == nextEdgeStart)
			{
				valid = true;
			}
			if (!valid)
			{
				cout << "Errore nella faccia " << i << ": lato " << currentEdgeId
				<< " termina in " << currentEdgeEnd << ", ma lato " << nextEdgeId << " inizia in " << nextEdgeStart << endl;
				return false;
			}
		}
		
		
	}
	return true;
}

Polyhedral DualPolygon (const Polyhedral& polygon)
{
	Polyhedral dual;
	
	dual.NumCell0Ds = polygon.NumCell2Ds;
	dual.Cell0DsCoordinates.resize(dual.NumCell0Ds,3);
	
	for (unsigned int i = 0; i < polygon.NumCell2Ds; i++)
	{
		//Per ogni faccia, devo considerare i vertici, e poi trovare il baricentro
		dual.Cell0DsId.push_back(i);
		
		const vector<unsigned int>& faceVertexId = polygon.Cell2DsVertices[i];
		
		MatrixXd faceVertices(faceVertexId.size(),3); //4 righe e 3 colonne!!
		
		
		for (unsigned int j = 0; j < faceVertexId.size(); j++)
			faceVertices.row(j) = polygon.Cell0DsCoordinates.row(faceVertexId[j]);
		
		Vector3d sum = faceVertices.colwise().sum();
		
		Vector3d barycenter = sum / faceVertexId.size();
		dual.Cell0DsCoordinates.row(i) = barycenter.transpose();
	}
	
	set<pair<unsigned int,unsigned int>> addedEdges;
	
	for (unsigned int i = 0; i < polygon.NumCell2Ds; i ++)
	{
		const vector<unsigned int>& edgesF = polygon.Cell2DsEdges[i];
		
		for (unsigned int j = i+1; j < polygon.NumCell2Ds; j++)
		{
			const vector<unsigned int>& edgesG = polygon.Cell2DsEdges[j];
			
			for (unsigned int e1 : edgesF)
			{
				for (unsigned int e2 : edgesG)
				{
					if (e1 == e2)
					{
						if (addedEdges.find({i,j}) == addedEdges.end())
							{
								addedEdges.insert({i,j});
								dual.Cell1DsId.push_back(dual.NumCell1Ds++);
								dual.Cell1DsVertices.conservativeResize(dual.NumCell1Ds,2);
								dual.Cell1DsVertices.row(dual.NumCell1Ds - 1) << i,j;
							}
					}
				}
			}
		}
	}
	
	dual.NumCell2Ds = polygon.NumCell0Ds;
	
	for (unsigned int v = 0; v < polygon.NumCell0Ds; v++)
	{
		vector<unsigned int> adjacentFaces;
		
		for (unsigned int f = 0; f < polygon.NumCell2Ds; f++ )
		{
			const auto& faceVertices = polygon.Cell2DsVertices[f];
			if (find(faceVertices.begin(),faceVertices.end(),v) != faceVertices.end())
				adjacentFaces.push_back(f);
		}
		
		vector<unsigned int> ordered;
		unordered_set<unsigned int> visited;
		ordered.push_back(adjacentFaces[0]);
		visited.insert(adjacentFaces[0]);
		
		while(ordered.size() < adjacentFaces.size())
		{
			unsigned int last = ordered.back();
			const auto& edgesF = polygon.Cell2DsEdges[last];
			bool found = false;
			for(unsigned int f : adjacentFaces)
			{
				if (visited.count(f)) 
					continue;
				const auto& edgesG = polygon.Cell2DsEdges[f];
				for (unsigned int eF : edgesF)
				{
					if (find(edgesG.begin(),edgesG.end(), eF) != edgesG.end())
					{
						ordered.push_back(f);
						visited.insert(f);
						found = true;
						break;
					}
				}
				if (found)
					break;
			}
			if (!found)
				break;
		}
		
		dual.Cell2DsVertices.push_back(ordered);
		dual.Cell2DsId.push_back(v);
	}
	
	for (unsigned int f = 0; f < dual.Cell2DsVertices.size() ; f++)
	{
		
		vector<unsigned int>& face_vertices = dual.Cell2DsVertices[f];
		vector<unsigned int> face_edges;
		for (unsigned int i = 0; i < face_vertices.size(); i++)
		{
			unsigned int vert_1 = face_vertices[i];
			unsigned int vert_2 = face_vertices[(i+1)%face_vertices.size()];
			
			for (unsigned int k = 0; k < dual.Cell1DsVertices.rows(); k++)
			{
				unsigned int coord_origin = dual.Cell1DsVertices(k,0);
				unsigned int coord_end = dual.Cell1DsVertices(k,1);
				if ((coord_origin == vert_1 && coord_end == vert_2) || (coord_origin == vert_2 && coord_end == vert_1))
				{
					face_edges.push_back(k);
					break;
				}
			}
			
		}
		dual.Cell2DsEdges.push_back(face_edges);
	}
	
    dual.NumCell3Ds = 1;
    dual.Cell3DsId.push_back(0);
	dual.Cell3DsVertices.push_back(dual.Cell0DsId);
	dual.Cell3DsEdges.push_back(dual.Cell1DsId);
	dual.Cell3DsFaces.push_back(dual.Cell2DsId);
	
	return dual;
}
			
bool GeodesicPolyhedron(const unsigned int& p, const unsigned int& q, const unsigned int& b, const unsigned int& c, Polyhedral& poly)
{
	if (p != 3)
		return false;

	ofstream file_c0("Cell0Ds.txt");
	
	if (file_c0.fail())
	{
		cout << "Errore nell'apertura del file" << endl;
		return false;
	}
	
	file_c0 << "ID\tx\t\ty\t\tz\n" ;
	
	poly.Cell0DsId.clear();

	vector<vector<Vector3d>> SubdividedVertices(poly.NumCell2Ds); //Ne avremo 4 per il tetraedro, uno per ogni faccia
	MatrixXd original_coords = poly.Cell0DsCoordinates; //Creo una copia dei vertici originali
	
	vector<vector<vector<unsigned int>>> grids(poly.NumCell2Ds);
	
	unsigned int ID = 0;

	vector<Vector3d> ordered_vertices;
	map<Vector3d,unsigned int, Vector3dCompare> vertices_id;
	
	for (unsigned int f = 0; f < poly.NumCell2Ds; f++)
	{
		vector<unsigned int> face = poly.Cell2DsVertices[f]; //Es: {0,1,2}
		vector<Vector3d> face_vertices; 
		vector<vector<unsigned int>> grid(max(b,c)+1); //griglia in cui salvo gli id dei miei nuovi vertici
		
		for (unsigned int i = 0; i < max(b,c) + 1; i++)
		{
			grid[i].resize(max(b,c) - i + 1);
		}
		
		Vector3d A = original_coords.row(face[0]).transpose();
		Vector3d B = original_coords.row(face[1]).transpose();
		Vector3d C = original_coords.row(face[2]).transpose();
		
		for (unsigned int i = 0; i < max(b,c) + 1; i ++)
		{
			for (unsigned int j = 0; j < max(b,c) - i + 1; j++)
			{
				unsigned int k = max(b,c) - i - j;
				Vector3d P = (i * A + j * B + k * C) / double(max(b,c)); //Coordinate baricentriche
				P.normalize();
				face_vertices.push_back(P);
				
				auto it = vertices_id.find(P);
				if (it == vertices_id.end())
				{
					vertices_id[P] = ID;
					ordered_vertices.push_back(P);
					ID ++;
				}
				unsigned int idx = vertices_id[P];

				grid[i][j] = idx;
			}
		}
		
		grids[f] = grid;
		SubdividedVertices[f] = face_vertices;

	}
	
	poly.Cell0DsCoordinates.resize(ordered_vertices.size(),3);
		
	for (size_t i = 0; i < ordered_vertices.size() ; i++) 
	{	
		poly.Cell0DsCoordinates.row(i) = ordered_vertices[i];
		
		const Vector3d& v = ordered_vertices[i];
		file_c0 << i << "\t" << v(0) << "\t" << v(1) << "\t" << v(2) << "\n";
		
	}
	
	poly.NumCell0Ds = poly.Cell0DsCoordinates.rows();

	file_c0.close();

	ofstream file_c1("Cell1Ds.txt");
	
	if (file_c1.fail())
	{
		cout <<"Errore nell'apertura del file" << endl;
		return false;
	}
	
	file_c1 << "ID\tIDOrigin\tIDEnd\n";
	
	set<pair<unsigned int,unsigned int>> edges;
	
	for (unsigned int f = 0; f < poly.NumCell2Ds ; f++)
	{
		const auto& grid = grids[f];
		
		for(unsigned int i = 0; i < grid.size(); i++)
		{
			for (unsigned int j = 0; j < grid[i].size(); j++)
			{
				
				unsigned int v = grid[i][j];
				
				if (j+1 < grid[i].size())
					edges.insert(minmax(v,grid[i][j+1]));
				if (i+1 < grid.size() && j < grid[i+1].size())
					edges.insert(minmax(v,grid[i+1][j]));
				if (i+1 < grid.size() && j+1 < grid[i+1].size()-1)
					edges.insert(minmax(v,grid[i+1][j+1]));
				if (i+1 < grid.size() && j > 0 && j - 1 < grid [i+1].size())
					edges.insert(minmax(v,grid[i+1][j-1]));
			}
		}
	}
	
	unsigned int edgeID = 0;

	for (const auto& [a, b] : edges) {
		file_c1 << edgeID++ << "\t" << a << "\t\t" << b << "\n";
	}
	poly.Cell1DsVertices.resize(edges.size(),2);
	unsigned int v = 0;
	for (const auto& [a,b] : edges)
	{
			poly.Cell1DsVertices(v,0) = a;
			poly.Cell1DsVertices(v,1) = b;
			v++;
	}

	poly.Cell1DsId.clear();
	
	for (unsigned int i = 0; i < edges.size(); i++)
		poly.Cell1DsId.push_back(i); 

	for (unsigned int i = 0; i < poly.NumCell1Ds; i++)
		poly.MarkerCell1Ds[1].push_back(i);
	
	file_c1.close();

	ofstream file_c2("Cell2Ds.txt");
	
	if (file_c2.fail())
	{
		cout <<"Errore nell'apertura del file" << endl;
		return false;
	}
	
	file_c2 << "ID\tNumVertices\tNumEdges\tVertices\tEdges\n";
	
	poly.NumVertices = 3;
	poly.NumEdges = 3;
	
	poly.Cell2DsVertices.clear();
	
	for (unsigned int f = 0; f < poly.NumCell2Ds; f++)
	{
		const auto& grid = grids[f];
		
		for (unsigned int i = 0; i + 1 < grid.size(); i++)
		{
			for (unsigned int j = 0; j < grid[i].size(); j++)
			{
				if (j < grid[i+1].size())
				{
					unsigned int A = grid[i][j];
					unsigned int B = grid[i+1][j];
					unsigned int C = grid[i][j+1];
					poly.Cell2DsVertices.push_back({A,B,C});
				}
				
				if (j+1 < grid[i+1].size())
				{
					unsigned int A = grid[i][j+1];
					unsigned int B = grid[i+1][j];
					unsigned int C = grid[i+1][j+1];
					poly.Cell2DsVertices.push_back({A,B,C});
					
				}
			}
		}
	}
	
	poly.Cell2DsEdges.clear();
	
	for (unsigned int f = 0; f < poly.Cell2DsVertices.size() ; f++)
	{
		
		vector<unsigned int>& face_vertices = poly.Cell2DsVertices[f];
		vector<unsigned int> face_edges;
		for (unsigned int i = 0; i < face_vertices.size(); i++)
		{
			unsigned int vert_1 = face_vertices[i];
			unsigned int vert_2 = face_vertices[(i+1)%face_vertices.size()];
			
			for (unsigned int k = 0; k < poly.Cell1DsVertices.rows(); k++)
			{
				unsigned int coord_origin = poly.Cell1DsVertices(k,0);
				unsigned int coord_end = poly.Cell1DsVertices(k,1);
				if ((coord_origin == vert_1 && coord_end == vert_2) || (coord_origin == vert_2 && coord_end == vert_1))
				{
					face_edges.push_back(k);
					break;
				}
			}
			
		}
		poly.Cell2DsEdges.push_back(face_edges);
	}

	for (unsigned int faceID = 0; faceID < poly.Cell2DsVertices.size(); faceID ++ )
	{
		const auto& vertices = poly.Cell2DsVertices[faceID];
		const auto& edges = poly.Cell2DsEdges[faceID];
		
		file_c2 << faceID << "\t" << poly.NumVertices << "\t\t" << poly.NumEdges << "\t\t" ;
		for (const auto& v : vertices)
			file_c2 << v << " ";
		file_c2 << "\t\t";
		for (const auto& e : edges)
			file_c2 << e << " ";
		file_c2 << "\t" << "\n";
	}
	
	poly.NumCell2Ds = poly.Cell2DsVertices.size();
	
	file_c2.close();
	
	ofstream file_c3("Cell3Ds.txt");
	
	if (file_c3.fail())
	{
		cout << "Errore nell'apertura del file" << endl;
		return false;
	}
	
	

	return true;	
}

bool shortestPath(Polyhedral& polygon,unsigned int start, unsigned int end, vector<unsigned int>& path)
{
	int n = polygon.NumCell0Ds;
	vector<bool> visited(n,false);
	vector<unsigned int> predecessor(n,-1);
	
	vector<vector<unsigned int>> adj(n);
	
	for (unsigned int i = 0; i < polygon.Cell1DsVertices.rows(); i++)
	{
		unsigned int u = polygon.Cell1DsVertices(i,0);
		unsigned int v = polygon.Cell1DsVertices(i,1);
		adj[u].push_back(v);
		adj[v].push_back(u);
	}
	
	queue<unsigned int> Q;
	Q.push(start);
	visited[start] = true;

	Q.push(start);
	while(!Q.empty())
	{
		unsigned int u = Q.front();
		Q.pop();
		
		for (unsigned int w : adj[u])
		{
			if (!visited[w])
			{
				visited[w] = true;
				predecessor[w] = u;
				Q.push(w);
				if (w == end)
					break;
			}
		}
	}		
	
	if (!visited[end])
		return false;
	
	for (unsigned int at = end; at != -1; at = predecessor[at])
	{
		path.push_back(at);
	}
	
	reverse(path.end(),path.begin());
	
	std::cout << "Shortest path: ";
	for (int idx : path) {
    std::cout << idx << " ";
	}
	std::cout << std::endl;
	
	return true;
}
void exportToUCD(const Polyhedral &polygon, const string& basename)
{
	
	/// Per visualizzare online le mesh:
    /// 1. Convertire i file .inp in file .vtu con https://meshconverter.it/it
    /// 2. Caricare il file .vtu su https://kitware.github.io/glance/app/
	
	
	Gedim::UCDUtilities utilities;
    {
        vector<Gedim::UCDProperty<double>> pt_props(1);

        pt_props[0].Label = "Marker";
        pt_props[0].UnitLabel = "-";
        pt_props[0].NumComponents = 1;

        vector<double> markers(polygon.NumCell0Ds, 0.0);
        for(const auto &m : polygon.MarkerCell0Ds)
            for(const unsigned int id: m.second)
                markers[id] = m.first;

        pt_props[0].Data = markers.data();

        utilities.ExportPoints(basename + "_points.inp",
                               polygon.Cell0DsCoordinates.transpose(),
                               pt_props);
							   
	   
    }

    {

        vector<Gedim::UCDProperty<double>> e_props(1);

        e_props[0].Label = "Marker";
        e_props[0].UnitLabel = "-";
        e_props[0].NumComponents = 1;

        vector<double> m2(polygon.NumCell1Ds, 0.0);
        for(const auto &m : polygon.MarkerCell1Ds)
            for(const unsigned int id: m.second)
                m2[id] = m.first;

        e_props[0].Data = m2.data();

        utilities.ExportSegments(basename + "_edges.inp",
                                 polygon.Cell0DsCoordinates.transpose(),
                                 polygon.Cell1DsVertices.transpose(),
                                 {},
                                 e_props);
    }
}	
}	
	
	
	
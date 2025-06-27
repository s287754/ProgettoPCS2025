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

bool checkOrient(unsigned int& n, const MatrixXi& Vertices, const vector<vector<unsigned int>>& Edges) 
{
	for (size_t i = 0; i < n; i ++)
	{
		const auto& faceEdges = Edges[i]; //Prendo un vettore alla volta. Es: {0,3,1}
		size_t E = faceEdges.size();
		
		for (size_t e = 0; e < E; e++)
		{
			unsigned int currentEdgeId = faceEdges[e]; //considero un indice. Es 0
			unsigned int nextEdgeId = faceEdges [(e+1)%E]; //e considero il successivo Es 3
			unsigned int currentEdgeStart = Vertices(currentEdgeId,0); //All'interno della matrice dei vertici considero Origin
			unsigned int currentEdgeEnd = Vertices(currentEdgeId,1); //E End del vertice in cEI
			unsigned int nextEdgeStart = Vertices(nextEdgeId,0); //E Origin 
			unsigned int nextEdgeEnd = Vertices(nextEdgeId,1); //E End in nEI
			
			bool valid = false;
			//Controllo se la fine di un lato coincide con l'inizio del successivo a meno dell'orientamento
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
//Funzione che costruisce il duale
Polyhedral DualPolygon (const Polyhedral& polygon)
{
	Polyhedral dual;
	
	dual.NumCell0Ds = polygon.NumCell2Ds; 
	
	dual.Cell0DsCoordinates.resize(dual.NumCell0Ds,3);
	
	for (unsigned int i = 0; i < polygon.NumCell2Ds; i++)
	{
		//Per ogni faccia, devo considerare i vertici, e poi trovare il baricentro
		dual.Cell0DsId.push_back(i);
		
		const vector<unsigned int>& faceVertexId = polygon.Cell2DsVertices[i]; //Considero un vettore di Cell2DsVertices alla volta. Es : {0,1,2,3} per il cubo
		
		MatrixXd faceVertices(faceVertexId.size(),3); //4 righe e 3 colonne!! (per il cubo)
		
		
		for (unsigned int j = 0; j < faceVertexId.size(); j++)
			faceVertices.row(j) = polygon.Cell0DsCoordinates.row(faceVertexId[j]); //Vado a inserire i dati nella mia nuova matrice
		
		Vector3d sum = faceVertices.colwise().sum(); //Creo un vettore colonna in cui sommo tutte le colonne (sommo le coordinate di tutti i punti)
		
		Vector3d barycenter = sum / faceVertexId.size(); //Divido per il numero di vertici di ogni faccia e trovo il baricentro di ogni faccia
		dual.Cell0DsCoordinates.row(i) = barycenter.transpose(); //Inserisco nella matrice coordinate
	}
	
	set<pair<unsigned int,unsigned int>> addedEdges; //Evito duplicati per i lati nel duale
	
	for (unsigned int i = 0; i < polygon.NumCell2Ds; i ++)
	{
		const vector<unsigned int>& edgesF = polygon.Cell2DsEdges[i]; //Considero una faccia alla volta. Es: {0,1,2,3}
		
		for (unsigned int j = i+1; j < polygon.NumCell2Ds; j++) 
		{
			const vector<unsigned int>& edgesG = polygon.Cell2DsEdges[j]; //Considero le facce successive a quella considerata. Es: {4,5,6,7}
			
			for (unsigned int e1 : edgesF) //Prima faccia : controllo ogni elemento
			{
				for (unsigned int e2 : edgesG) //Seconda faccia : controllo ogni elemento
				{
					if (e1 == e2) //Controllo che condividano almeno un bordo e siano quindi adiacenti
					{
						if (addedEdges.find({i,j}) == addedEdges.end()) //Controllo che non sia stato già trattato
							{
								addedEdges.insert({i,j}); //Se non presente, aggiungilo
								dual.Cell1DsId.push_back(dual.NumCell1Ds++); //Aggiungo l'ID
								dual.Cell1DsVertices.conservativeResize(dual.NumCell1Ds,2);
								dual.Cell1DsVertices.row(dual.NumCell1Ds - 1) << i,j; //Aggiungo i nodi del duale alla matrice dei Vertici 
							}
					}
				}
			}
		}
	}
	
	dual.NumCell2Ds = polygon.NumCell0Ds;
	
	for (unsigned int v = 0; v < polygon.NumCell0Ds; v++) //Itero su tutti i vertici del poligono originale
	{
		vector<unsigned int> adjacentFaces; //Creo un vettore che memorizzerà le facce che contengono il vertice v (saranno adiacenti)
		
		for (unsigned int f = 0; f < polygon.NumCell2Ds; f++ )
		{
			const auto& faceVertices = polygon.Cell2DsVertices[f]; //Considero una faccia alla volta. Es: {0,1,2,3}
			if (find(faceVertices.begin(),faceVertices.end(),v) != faceVertices.end()) //Cerco tutte le facce che contengono il vertice v
				adjacentFaces.push_back(f); //Memorizzo tutte le facce che contengono v nel nuovo vettore
		}
		
		vector<unsigned int> ordered; //Mi serve questo vettore per mettere in ordine i vertici della faccia
		unordered_set<unsigned int> visited; //Mi serve per evitare ripetizioni
		ordered.push_back(adjacentFaces[0]); //Parto dalla prima faccia
		visited.insert(adjacentFaces[0]);
		//Nel ciclo ordino tutte le facce
		while(ordered.size() < adjacentFaces.size())
		{
			unsigned int last = ordered.back(); //Considero l'ultima faccia aggiunta 
			const auto& edgesF = polygon.Cell2DsEdges[last]; //Considero i lati di questa ultima faccia aggiunta
			bool found = false;
			for (unsigned int f : adjacentFaces) 
			{
				if (visited.count(f)) //Se la faccia f è già stata visitata, salto
					continue;
				const auto& edgesG = polygon.Cell2DsEdges[f]; //Prendo i lati della faccia f
				for (unsigned int eF : edgesF) //Per ogni elemento dell'ultima faccia aggiunta
				{
					if (find(edgesG.begin(),edgesG.end(), eF) != edgesG.end()) //Cerco un bordo comune : se lo trovo, le facce sonoa diacenti
					{
						ordered.push_back(f); //Inserisco in ordered
						visited.insert(f); //Inserisco in visited
						found = true;
						break;
					}
				}
				if (found) //Ho trovato la faccia da aggiungere quindi interrompo il ciclo
					break;
			}
			if (!found)
				break;
		}
		
		dual.Cell2DsVertices.push_back(ordered);
		dual.Cell2DsId.push_back(v);
	}
	//Creo adesso Cell2DsEdges
	for (unsigned int f = 0; f < dual.Cell2DsVertices.size() ; f++)
	{
		
		vector<unsigned int>& face_vertices = dual.Cell2DsVertices[f]; //Considero una faccia alla volta
		vector<unsigned int> face_edges; //Vettore locale
		for (unsigned int i = 0; i < face_vertices.size(); i++)
		{
			unsigned int vert_1 = face_vertices[i]; //Prendo un vertice*
			unsigned int vert_2 = face_vertices[(i+1)%face_vertices.size()]; //Prendo il successivo*
			
			for (unsigned int k = 0; k < dual.Cell1DsVertices.rows(); k++)
			{
				unsigned int coord_origin = dual.Cell1DsVertices(k,0); //Prendo Origin di un vertice
				unsigned int coord_end = dual.Cell1DsVertices(k,1); //Prendo End dello stesso vertice
				if ((coord_origin == vert_1 && coord_end == vert_2) || (coord_origin == vert_2 && coord_end == vert_1)) //Se l'origin o l'end del vertice coincide con vertice*
				{
					face_edges.push_back(k); //Inserisco nel vettore locale
					break;
				}
			}
			
		}
		dual.Cell2DsEdges.push_back(face_edges); //Salvo ogni vettore locale in Cell2DsEdges
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
	unsigned int real_p = poly.Cell2DsVertices[0].size();
	if (p != real_p)
		return false;
	
	(void) q; //Mi serve nel main, per decidere quale poliedro suddividere, non qui
	
	ofstream file_c0("Cell0Ds.txt");
	
	if (file_c0.fail())
	{
		cout << "Errore nell'apertura del file" << endl;
		return false;
	}
	
	file_c0 << "ID\tx\t\ty\t\tz\n" ;
	
	MatrixXd original_coords = poly.Cell0DsCoordinates; //Creo una copia dei vertici originali
	
	vector<vector<vector<unsigned int>>> grids(poly.NumCell2Ds);
	
	unsigned int ID = 0;

	vector<Vector3d> ordered_vertices; //Memorizzo i vettori ordinati
	map<Vector3d,unsigned int, Vector3dCompare> vertices_id; //Associo a ogni vertice un ID, se il vertice non è presente.
	
	for (unsigned int f = 0; f < poly.NumCell2Ds; f++)
	{
		vector<unsigned int> face = poly.Cell2DsVertices[f]; //Es: {0,1,2}
		vector<vector<unsigned int>> grid(max(b,c)+1); //griglia in cui salvo gli id dei miei nuovi vertici
		
		for (unsigned int i = 0; i < max(b,c) + 1; i++)
		{
			grid[i].resize(max(b,c) - i + 1); //Ogni riga della griglia ha una dimensione diversa(in base al numero di nuovi vertici)
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
				P.normalize(); //proiezione su sfera
				//face_vertices.push_back(P); 
				
				auto it = vertices_id.find(P); //Controllo che il vertice non sia già stato generato
				if (it == vertices_id.end())
				{
					vertices_id[P] = ID; //vertice non ancora generato : assegna un nuovo ID
					ordered_vertices.push_back(P); //inserisco il nuovo vertice nel vettore di vertici ordinati
					ID ++;
				}
				unsigned int idx = vertices_id[P]; //Recupera l'ID del vertice (nuovo o già generato in precedenza)

				grid[i][j] = idx; //inserisco l'ID nella griglia triangolare
			}
		}
		
		grids[f] = grid;
		
	}
	
	poly.Cell0DsCoordinates.resize(ordered_vertices.size(),3);
	poly.Cell0DsId.clear();
	
	for (size_t i = 0; i < ordered_vertices.size() ; i++) 
	{	
		poly.Cell0DsCoordinates.row(i) = ordered_vertices[i]; //Inserisco nella struttura tutti i vertici nuovi prodotti dalla geodetizzazione una e una sola volta
		poly.Cell0DsId.push_back(i);
		const Vector3d& v = ordered_vertices[i];
		file_c0 << poly.Cell0DsId[i] << "\t" << v(0) << "\t" << v(1) << "\t" << v(2) << "\n"; //Stampo nel file
	}

	poly.NumCell0Ds = poly.Cell0DsCoordinates.rows(); //Aggiorno la variabile

	file_c0.close();

	ofstream file_c1("Cell1Ds.txt");
	//Genero il file Cell1Ds.txt
	if (file_c1.fail())
	{
		cout <<"Errore nell'apertura del file" << endl;
		return false;
	}
	
	file_c1 << "ID\tIDOrigin\tIDEnd\n";
	
	set<pair<unsigned int,unsigned int>> edges;
	
	for (unsigned int f = 0; f < poly.NumCell2Ds ; f++)                                            //La mia griglia sarà così fatta (Es per il triangolo con c = 3)
	{
		const auto& grid = grids[f];                                                               // 0   1   2   3
		                                                                                           // 4   5   6
		for(unsigned int i = 0; i < grid.size(); i++)                                              // 7   8
		{                                                                                          // 9           
			for (unsigned int j = 0; j < grid[i].size(); j++)
			{
				
				unsigned int v = grid[i][j];
				
				if (i+1 < grid.size() && j < grid[i+1].size())
					edges.insert(minmax(v,grid[i+1][j])); //minmax serve a inserire gli edges in ordine (origin sarà sempre minore di end)
				//Collegamento in basso
				if (j+1 < grid[i].size())
					edges.insert(minmax(v,grid[i][j+1])); //Collegamento a dx
				
				if (j > 0 && i+1 < grid.size() && j-1 < grid[i+1].size())
					edges.insert(minmax(v,grid[i+1][j-1])); //Collegamento in basso a sx

			}
		}
	}
	
	unsigned int edgeID = 0;

	for (const auto& [a, b] : edges) {
		file_c1 << edgeID++ << "\t" << a << "\t\t" << b << "\n";
	}
	//Inserisco tutto all'interno di Cell1DsVertices
	poly.Cell1DsVertices.resize(edges.size(),2);
	unsigned int v = 0;
	for (const auto& [a,b] : edges)
	{
			poly.Cell1DsVertices(v,0) = a;
			poly.Cell1DsVertices(v,1) = b;
			v++;
	}

	poly.NumCell1Ds = poly.Cell1DsVertices.rows();
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
	
	poly.Cell2DsId.clear();
	poly.Cell2DsVertices.clear();
	//Creo la triangolazione
	for (unsigned int f = 0; f < poly.NumCell2Ds; f++)
	{
		const auto& grid = grids[f];
		
		for (unsigned int i = 0; i + 1 < grid.size(); i++)
		{
			for (unsigned int j = 0; j < grid[i].size(); j++)
			{
				if (j < grid[i+1].size())
				{
					unsigned int A = grid[i][j]; //0
					unsigned int B = grid[i+1][j]; //4
					unsigned int C = grid[i][j+1]; //1
					poly.Cell2DsVertices.push_back({A,B,C});
				}
				
				if (j+1 < grid[i+1].size())
				{
					unsigned int A = grid[i][j+1]; //1
					unsigned int B = grid[i+1][j]; //4
					unsigned int C = grid[i+1][j+1]; //5
					poly.Cell2DsVertices.push_back({A,B,C});
					
				}
			}
		}
	}
	
	poly.Cell2DsEdges.clear();
	//Come nel duale costruisco Cell2DsEdges
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
	
	poly.Cell2DsId.clear();
	for (unsigned int i = 0; i < poly.Cell2DsVertices.size(); i++)
		poly.Cell2DsId.push_back(i);
	
	poly.NumCell2Ds = poly.Cell2DsVertices.size();
	
	file_c2.close();
	
	ofstream file_c3("Cell3Ds.txt");
	
	if (file_c3.fail())
	{
		cout << "Errore nell'apertura del file" << endl;
		return false;
	}
	
	file_c3 << "ID\tNumVertices\tNumEdges\tNumFaces\tVertices\tEdges\t\tFaces\n";

	poly.NumCell3Ds = 1;
	
	poly.Cell3DsId.clear();
	poly.Cell3DsId.push_back(0);
	
	poly.Cell3DsVertices.clear();
	poly.Cell3DsVertices.push_back(poly.Cell0DsId);
	
	poly.Cell3DsEdges.clear();
	poly.Cell3DsEdges.push_back(poly.Cell1DsId);
	
	poly.Cell3DsFaces.clear();
	poly.Cell3DsFaces.push_back(poly.Cell2DsId);
	
	for (unsigned int i = 0; i < poly.NumCell3Ds; i++)
	{
		file_c3 << poly.Cell3DsId[i] << "\t" << poly.Cell3DsVertices[i].size() << "\t\t" 
		<< poly.Cell3DsEdges[i].size() << "\t\t" << poly.Cell3DsFaces[i].size() << "\t\t" ;
		
		for (const auto& v: poly.Cell3DsVertices[i])
			file_c3 << v << " ";
		file_c3 << "\t";
		for (const auto& e: poly.Cell3DsEdges[i])
			file_c3 << e <<" ";
		file_c3 << "\t";
		for (const auto& f: poly.Cell3DsFaces[i])
			file_c3 << f << " ";
		file_c3 << "\n";
	}

	return true;	
}
//Calcolo il cammino minimo con BFS
bool shortestPath(Polyhedral& polygon,unsigned int start, unsigned int end, vector<unsigned int>& path,vector<unsigned int>& crossedEdges)
{
	int n = polygon.NumCell0Ds; //Numero di vertici della mesh
	vector<bool> visited(n,false); //Nodi già visitati
	vector<unsigned int> predecessor(n,-1); //Salvo il nodo da cui arrivo per ricostruire il percorso
	
	vector<vector<unsigned int>> adj(n); //Memorizzo i nodi adiacenti per ogni nodo
	
	for (unsigned int i = 0; i < polygon.Cell1DsVertices.rows(); i++)
	{
		unsigned int u = polygon.Cell1DsVertices(i,0); //Origine del lato
		unsigned int v = polygon.Cell1DsVertices(i,1); //End del lato
		adj[u].push_back(v); //Grafo non orientato
		adj[v].push_back(u);
	}
	
	queue<unsigned int> Q; //creo la cosa
	Q.push(start);
	visited[start] = true; //Parto dal primo nodo (lo chiedo in input)

	while(!Q.empty())
	{
		unsigned int u = Q.front(); //Considera il primo nodo della coda
		Q.pop(); //Eliminalo dalla coda
		
		for (unsigned int w : adj[u]) //Considera tutti i nodi adiacenti a u
		{
			if (!visited[w]) //Se non è stato visitato
			{
				visited[w] = true; //Lo considero visitato
				predecessor[w] = u; //Salvo il predecessore 
				Q.push(w); //Lo inserisco nella coda
				if (w == end)
					break;
			}
		}
	}		
	
	if (!visited[end])
		return false;
	
	for (int at = end; at != -1; at = predecessor[at])
	{
		path.push_back(at);
	}
	
	reverse(path.end(),path.begin());
	
	std::cout << "Shortest path: ";
	for (int idx : path) {
    std::cout << idx << " ";
	}
	std::cout << std::endl;
	//Per ogni coppia di vertici cerco il lato che li connette effettivamente
	for (size_t i = 0; i < path.size()-1; i++)
	{
		unsigned int u = path[i];
		unsigned int v = path[i+1];
		bool found = false;
		
		for (unsigned int j = 0; j < polygon.Cell1DsVertices.rows(); j++)
		{
			unsigned int Origin = polygon.Cell1DsVertices(j,0);
			unsigned int End = polygon.Cell1DsVertices(j,1);
			
			if ((Origin == u && End == v) || (Origin == v && End == u))
			{
				crossedEdges.push_back(j); //Se il lato collega i due vertici della mesh, lo salvo
				found = true;
				break;
			}
		}
		
		if (!found)
		{
			cout << "Nessun percorso minimo trovato tra il nodo " << u << " e il nodo" << v << endl;
			return false;
		}
	}

	std::cout << "Lati attraversati (ID): ";
	for (unsigned int edgeID : crossedEdges)
	{
		std::cout << edgeID << " ";
	}
	std::cout << std::endl;
	
	return true;
}

void exportToUCD(const Polyhedral &polygon, const string& basename,vector<unsigned int>& path, vector<unsigned int>& crossedEdges)
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
	
	if (crossedEdges.size() >= 1)
	{
		unordered_map<unsigned int, unsigned int> vertexToIndex;
		
		vector<unsigned int> uniqueVertices;
		
		vector<Vector2i> edges;
		
		for(unsigned int edgeID : crossedEdges)
		{
			unsigned int u = polygon.Cell1DsVertices(edgeID,0);
			unsigned int v = polygon.Cell1DsVertices(edgeID,1);
			
			if (vertexToIndex.find(u) == vertexToIndex.end())
			{
				vertexToIndex[u] = uniqueVertices.size();
				uniqueVertices.push_back(u);
			}	
			if (vertexToIndex.find(v) == vertexToIndex.end())
			{
				vertexToIndex[v] = uniqueVertices.size();
				uniqueVertices.push_back(v);
			}	
				
			edges.push_back(Vector2i(vertexToIndex[u],vertexToIndex[v]));
		}
			
		MatrixXd pathCoordinates (path.size(),3);
		
		for (size_t i = 0; i < path.size(); i++)
			pathCoordinates.row(i) = polygon.Cell0DsCoordinates.row(uniqueVertices[i]);

		MatrixXi edgesMatrix(edges.size(),2);
		
		for (size_t i = 0; i < edges.size(); i++)
		{
			edgesMatrix.row(i) = edges[i];
		}
		
		vector<Gedim::UCDProperty<double>> pt_props_path;
		utilities.ExportPoints(basename + "_path_points.inp",
							   pathCoordinates.transpose(), 
							   pt_props_path);
							   
		vector<Gedim::UCDProperty<double>> e_props_path;
		utilities.ExportSegments(basename + "_path_edges.inp",
								 pathCoordinates.transpose(),
								 edgesMatrix.transpose(), 
								 {},
								 e_props_path);
	}
}	
}
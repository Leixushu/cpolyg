#include "PolyMesh.h"
#include <iostream>
#include "voro++.hh"

using namespace voro;
using namespace std;

// This function returns a random floating point number between 0 and 1
double rnd() {return double(rand())/RAND_MAX;}

PolyMesh::PolyMesh(vector<array<double, 2>> points)
{
    double x,y;
    
	container con(0, 1, 0, 1, 0, 1, 6, 6, 1, false, false, false, 8);
    
    cout << "Constructing Voronoi Diagram from " << points.size()
         << " generating points." << endl;
    
    for(int i = 0; i < points.size(); i++)
    {
        x = points[i][0];
        y = points[i][1];
        con.put(i,x,y,0);
    }
    
    
    con.draw_particles("plt/points.gnu");
    con.draw_cells_gnuplot("plt/edges.gnu");
}
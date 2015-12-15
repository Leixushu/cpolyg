#include "PolyMesh.h"
#include <iostream>
#include <voro++_2d.hh>
#include <set>
#include <vector>
#include <algorithm>
#include <iterator>

using namespace voro;
using namespace std;

#define EPS 1.e-10

template<typename T>
std::ostream & operator<<(std::ostream & os, std::vector<T> vec)
{
    os<<"{ ";
    std::copy(vec.begin(), vec.end(), std::ostream_iterator<T>(os, ", "));
    os<<"}";
    return os;
}

void PolyMesh::computep2p(void)
{
    int i, j, k, l;
    int a1, b1, a2, b2;
    int nv1, nv2;
    
    p2p.resize(np);
    
    // loop over all polygons
    for (i = 0; i < np; i++)
    {
        nv1 = p[i].size();
        // loop over all edges
        for (j = 0; j < nv1; j++)
        {
            a1 = p[i][j];
            b1 = p[i][(j+1)%nv1];
            
            for (k = i+1; k < np; k++)
            {
                nv2 = p[k].size();
                // loop over all edges
                for (l = 0; l < nv2; l++)
                {
                    a2 = p[k][l];
                    b2 = p[k][(l+1)%nv2];
                    
                    if ((a1 == a2 && b1 == b2) || (a1 == b2 && b1 == a2))
                    {
                        p2p[i].push_back(k);
                        p2p[k].push_back(i);
                        break;
                    }
                    
                }
            }
            
        }
    }
}

void PolyMesh::computebb()
{
    int i, j;
    int nv;
    array<double, 4> box;
    
    bb.resize(np);
    
    for (i = 0; i < np; i++)
    {
        nv = p[i].size();
        
        for (j = 0; j < nv; j++)
        {
            if (j == 0 || v[p[i][j]][0] < box[0])
                box[0] = v[p[i][j]][0];
            if (j == 0 || v[p[i][j]][1] < box[1])
                box[1] = v[p[i][j]][1];
            if (j == 0 || v[p[i][j]][0] > box[2])
                box[2] = v[p[i][j]][0];
            if (j == 0 || v[p[i][j]][1] > box[3])
                box[3] = v[p[i][j]][1];
        }
        
        bb[i][0] = box[0];
        bb[i][1] = box[1];
        bb[i][2] = box[2] - box[0];
        bb[i][3] = box[3] - box[1];
    }
}

void PolyMesh::computeTriangulation()
{
    int i, j;
    int nv;
    
    for (i = 0; i < np; i++)
    {
        nv = p[i].size();
        
        for (j = 0; j < nv; j++)
        {
        
        }
    }
}

int PolyMesh::addVertex(array<double, 2> vertex)
{
    int i;
    
    for (i = 0; i < v.size(); i++)
    {
        if (fabs(vertex[0] - v[i][0]) + fabs(vertex[1] - v[i][1]) < EPS)
        {
            return i;
        }
    }
    
    v.push_back(vertex);
    return v.size()-1;
}

PolyMesh::PolyMesh(vector<array<double, 2>> points)
{
    int i;
    double x,y;
    array<double, 2> vertex;
    vector<int> polygon;
    
	container_2d con(0, 1, 0, 1, 1, 1, false, false, 8);
	voronoicell_2d cell;
    
    cout << "Constructing Voronoi Diagram from " << points.size()
         << " generating points." << endl;
    
    cout << "ok..." << endl;
    
    for(i = 0; i < points.size(); i++)
    {
        x = points[i][0];
        y = points[i][1];
        con.put(i,x,y);
    }
    
    con.draw_particles("plt/points.gnu");
    con.draw_cells_gnuplot("plt/edges.gnu");
    
    for(i = 0; i < points.size(); i++)
    {
        if(con.compute_cell(cell, 0, i))
        {
            polygon.clear();
            
            int k=0;
            do
            {
                vertex[0] = points[i][0] + 0.5*cell.pts[2*k];
                vertex[1] = points[i][1] + 0.5*cell.pts[2*k + 1];
                
                int vi = addVertex(vertex);
                
                if(std::find(polygon.begin(), polygon.end(), vi) == polygon.end())
                {
                    polygon.push_back(vi);
                }
                
                k = cell.ed[2*k];
            } while (k!=0);
            
            p.push_back(polygon);
        }
    }
    
    cout << "Resulting in " << v.size() << " vertices, "
                            << p.size() << " polygons." << endl;
    
    np = p.size();
    computep2p();
    computebb();
    computeTriangulation();
}
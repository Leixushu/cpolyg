#include "PolyMesh.h"
#include <voro++_2d.hh>
#include <iostream>
#include <fstream>
#include <cassert>

#include "Quadrature.h"

using namespace voro;
using namespace std;

void PolyMesh::computep2p(void)
{
    int i, j, k, l;
    int a1, b1, a2, b2;
    int nv1, nv2;
    bool found;
    BoundaryInfo info;
    
    info.periodic = false;
    bc.emplace(-1, info);
    
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
            
            found = false;
            for (k = 0; k < np; k++)
            {
                if (i == k) continue;
                
                nv2 = p[k].size();
                // loop over all edges
                for (l = 0; l < nv2; l++)
                {
                    a2 = p[k][l];
                    b2 = p[k][(l+1)%nv2];
                    
                    if ((a1 == a2 && b1 == b2) || (a1 == b2 && b1 == a2))
                    {
                        p2p[i].push_back(k);
                        //p2p[k].push_back(i);
                        found = true;
                        break;
                    }
                }
                
                if (found) break;
            }
            
            // if this edge of the polygon has no neighboring polygons, it is on the 
            // exterior of the domain, so we mark it with a negative number
            if (!found)
            {
                p2p[i].push_back(-1);
            }
        }
    }
}

// compute and store the bounding boxes for all of polygons
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

// triangulate each polygon and store the triangulation
void PolyMesh::computeTriangulation()
{
    int i, j;
    int nv;
    vector<array<double, 2> > points;
    
    tri.resize(np);
    for (i = 0; i < np; i++)
    {
        nv = p[i].size();
        points.resize(nv);
        
        for (j = 0; j < nv; j++)
        {
            points[j][0] = v[p[i][j]][0];
            points[j][1] = v[p[i][j]][1];
        }
        
        tri[i] = Triangulation(points);
    }
}

int PolyMesh::findVertex(double x, double y)
{
    unsigned int i;
    
    for (i = 0; i < v.size(); i++)
    {
        if (fabs(x - v[i][0]) + fabs(y - v[i][1]) < kEPS)
        {
            return i;
        }
    }
    
    return -1;
}

// Add a vertex to the mesh and return its index. If the vertex already exists
// (to within a tolerance), don't add and return the index of the existing vertex
int PolyMesh::addVertex(array<double, 2> vertex)
{
    int i;
    
    i = findVertex(vertex[0], vertex[1]);
    
    if(i < 0)
    {
        v.push_back(vertex);
        return v.size()-1;
    } else
    {
        return i;
    }
}

PolyMesh::PolyMesh(vector<array<double, 2> > points, double width, double height)
{
    unsigned int i;
    double x,y;
    array<double, 2> vertex;
    vector<int> polygon;
    
	container_2d con(0, width+kEPS, 0, height+kEPS, 1, 1, false, false, 8);
	voronoicell_2d cell;
    
    cout << "Constructing Voronoi Diagram from " << points.size()
         << " generating points." << endl;
    
    for(i = 0; i < points.size(); i++)
    {
        x = points[i][0];
        y = points[i][1];
        con.put(i,x,y);
    }
    
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
    
    con.draw_particles("plt/vorpoints.gnu");
    con.draw_cells_gnuplot("plt/voredges.gnu");
    
    np = p.size();
    computep2p();
    computebb();
    computeTriangulation();
}

// Write out plot files for the mesh
void PolyMesh::gnuplot()
{
    ofstream edgeFile;
    int i, j, nv;
    
    edgeFile.open("plt/edges.gnu");
    
    for (i = 0; i < np; i++)
    {
        edgeFile << "# polygon number " << i << endl;
        nv = p[i].size();
        
        for (j = 0; j <= nv; j++)
        {
            edgeFile << v[p[i][j%nv]][0] << "\t" << v[p[i][j%nv]][1] << endl;
        }
        edgeFile << endl << endl;
    }
    
    edgeFile.close();
}

template <typename T> int sgn(T val)
{
    return (T(0) < val) - (val < T(0));
}

void PolyMesh::getOutwardNormal(int i, int a, int b, double &x, double &y)
{
    double dx, dy, cx, cy, length;
    unsigned int j;
    int s;
    
    cx = 0;
    cy = 0;
    for (j = 0; j < p[i].size(); j++)
    {
        cx += v[p[i][j]][0];
        cy += v[p[i][j]][1];
    }
    cx /= p[i].size();
    cy /= p[i].size();
    
    dx = v[a][0] - v[b][0];
    dy = v[a][1] - v[b][1];
    
    length = sqrt(dx*dx + dy*dy);
    
    s = sgn(-dy*(v[a][0] - cx) + dx*(v[a][1] - cy));
    
    x = -s*dy/length;
    y = s*dx/length;
}

PolyMesh PolyMesh::triangulate(std::vector<std::array<double, 2> > points)
{
    PolyMesh msh;
    int i, j;
    
    cout << "Constructing Delaunay Triangulation from " << points.size()
         << " generating points." << endl;
    
    Triangulation delaunay(points);
    
    msh.v = delaunay.points;
    msh.np = delaunay.triangles.size();
    msh.p.resize(msh.np, vector<int>(3));
    
    for (i = 0; i < msh.np; i++)
    {
        for (j = 0; j < 3; j++)
        {
            msh.p[i][j] = delaunay.triangles[i][j];
        }
    }
        
    cout << "Resulting in " << msh.v.size() << " vertices, "
                            << msh.p.size() << " polygons." << endl;
    
    msh.computep2p();
    msh.computebb();
    msh.computeTriangulation();
    
    return msh;
}

void PolyMesh::getPeriodicCoordinates(int p2, double xIn, double yIn,
                                      double &xOut, double &yOut) const
{
    double a1x, a1y, b1x, b1y, a2x, a2y, b2x, b2y, xNew, yNew;
    
    BoundaryInfo b = bc.at(p2);
    
    a1x = v[b.a1][0];
    a1y = v[b.a1][1];
        
    a2x = v[b.a2][0];
    a2y = v[b.a2][1];
        
    b1x = v[b.b1][0];
    b1y = v[b.b1][1];
        
    b2x = v[b.b2][0];
    b2y = v[b.b2][1];
    
    xNew = xIn + (a2x - a1x);
    yNew = yIn + (a2y - a1y);
            
    xOut = 2.0*(xNew - bb[b.p2][0])/bb[b.p2][2] - 1;
    yOut = 2.0*(yNew - bb[b.p2][1])/bb[b.p2][3] - 1;
}

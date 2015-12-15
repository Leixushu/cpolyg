#include "Triangulation.h"
#include <stdlib.h>
#include <iostream>

Triangulation::Triangulation(std::vector<std::array<double, 2> > pts)
{
    points = pts;
    doTriangulation();
}

void Triangulation::doTriangulation()
{
    struct triangulateio in, out;
    
    char switches[] = "zQ";
    int i;
    
    in.numberofpointattributes = 0;
    in.pointattributelist = NULL;
    in.pointmarkerlist = NULL;
    
    in.numberofpoints = points.size();
    in.pointlist = (double *)malloc(sizeof(double)*2*points.size());
    
    for (i = 0; i < points.size(); i++)
    {
        in.pointlist[2*i] = points[i][0];
        in.pointlist[2*i+1] = points[i][1];
    }
    
    in.trianglelist = NULL;
    in.triangleattributelist = NULL;
    in.trianglearealist = NULL;
    in.neighborlist = NULL;
    in.numberoftriangles = 0;
    in.numberofcorners = 0;
    in.numberoftriangleattributes = 0;
    
    in.segmentlist = NULL;
    in.segmentmarkerlist = NULL;
    in.numberofsegments = 0;
    
    in.holelist = NULL;
    in.numberofholes = 0;
    
    in.regionlist = NULL;
    in.numberofregions = 0;
    
    in.edgelist = 0;
    in.edgemarkerlist = NULL;
    in.normlist = NULL;
    
    out.pointlist = NULL;
    out.pointattributelist = NULL;
    out.pointmarkerlist = NULL;
    out.trianglelist = NULL;
    out.neighborlist = NULL;
    out.triangleattributelist = NULL;
    out.segmentlist = NULL;
    out.segmentmarkerlist = NULL;
    out.edgelist = NULL;
    out.edgemarkerlist = NULL;
    
    triangulate(switches, &in, &out, NULL);
    
    free(in.pointlist);
    
    points.resize(out.numberofpoints);
    for (i = 0; i < out.numberofpoints; i++)
    {
        points[i][0] = out.pointlist[2*i];
        points[i][1] = out.pointlist[2*i+1];
    }
    
    triangles.resize(out.numberoftriangles);
    p.resize(out.numberoftriangles);
    for (i = 0; i < out.numberoftriangles; i++)
    {
        triangles[i][0] = out.trianglelist[out.numberofcorners*i];
        triangles[i][1] = out.trianglelist[out.numberofcorners*i+1];
        triangles[i][2] = out.trianglelist[out.numberofcorners*i+2];
        
        p[i][0] = points[triangles[i][0]][0];
        p[i][1] = points[triangles[i][0]][1];
        p[i][2] = points[triangles[i][1]][0];
        p[i][3] = points[triangles[i][1]][1];
        p[i][4] = points[triangles[i][2]][0];
        p[i][5] = points[triangles[i][2]][1];
    }
    
    if(out.numberofpoints)
    {
        free(out.pointlist);
    }
    
    if(out.numberoftriangles)
    {
        free(out.trianglelist);
    }
}

#include <iostream>
#include "../debug/debug.h"
#include "unstructuredmesh.h"
using namespace std;

int main()
{

    UnstructuredMesh<2, size_t, double> mesh;

    mesh.GetVertices().resize(4);

    mesh.GetVertices().at(0).SetIndex(0);
    mesh.GetVertices().at(0) = {0,0};
    mesh.GetVertices().at(1).SetIndex(1);
    mesh.GetVertices().at(1) = {0,1};
    mesh.GetVertices().at(2).SetIndex(2);
    mesh.GetVertices().at(2) = {1,0};
    mesh.GetVertices().at(3).SetIndex(3);
    mesh.GetVertices().at(3) = {1,1};

    mesh.GetFaces().resize(5);
    mesh.GetFaces().at(0).VertexA = 0;
    mesh.GetFaces().at(0).VertexB = 1;
    mesh.GetFaces().at(1).VertexA = 0;
    mesh.GetFaces().at(1).VertexB = 2;
    mesh.GetFaces().at(2).VertexA = 1;
    mesh.GetFaces().at(2).VertexB = 2;
    mesh.GetFaces().at(3).VertexA = 1;
    mesh.GetFaces().at(3).VertexB = 3;
    mesh.GetFaces().at(4).VertexA = 2;
    mesh.GetFaces().at(4).VertexB = 3;
    for(size_t i = 0; i < 5; i++)
        mesh.GetFaces().at(i).SetIndex(i);

    mesh.GetFaces().at(0).SetNextBElem(1,0);
    mesh.GetFaces().at(1).SetNextBElem(2,0);
    mesh.GetFaces().at(2).SetNextBElem(0,0);

    mesh.GetFaces().at(2).SetNextBElem(3,1);
    mesh.GetFaces().at(3).SetNextBElem(4,1);
    mesh.GetFaces().at(4).SetNextBElem(2,1);


    mesh.GetCells().resize(2);
    mesh.GetCells().at(0).SetBoundaryElementIndex(0);
    mesh.GetCells().at(0).SetIndex(0);
    mesh.GetCells().at(1).SetBoundaryElementIndex(2);
    mesh.GetCells().at(1).SetIndex(1);

    mesh.InitializeCenters();

    for (size_t i = 0; i < 2; i++) {
        DBGMSG("Cell number" << i);
        size_t edge = mesh.GetCells().at(i).GetBoundaryElementIndex();
        size_t tmp_edge = edge;
        do {
            DBGVAR(tmp_edge, mesh.GetFaces().at(tmp_edge).VertexA,mesh.GetFaces().at(tmp_edge).VertexB,mesh.GetVertices().at(mesh.GetFaces().at(tmp_edge).VertexA), mesh.GetVertices().at(mesh.GetFaces().at(tmp_edge).VertexB), mesh.CalculateEdgeMeasure(tmp_edge));
            tmp_edge = mesh.GetFaces().at(tmp_edge).GetNextBElem(i);
        } while(edge != tmp_edge);

    }

    for(size_t i = 0; i < 2; i++) {
        DBGVAR(mesh.GetCells().at(i).GetCenter(), mesh.CalculateCellMeasure(i));
        DBGVAR(mesh.GetCells().at(i).GetFlag());
    }
    DBGVAR(mesh.CalculateEdgeMeasure(2)/mesh.CalculateCellDist(0,1));
    DBGVAR(mesh.CalculateFaceMeasureOverCellDist(2));

}

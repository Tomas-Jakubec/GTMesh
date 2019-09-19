#include <iostream>
#include "../debug/debug.h"
#include "unstructuredmesh.h"
#include "mesh_functions.h"
using namespace std;



int main()
{
    using sit = UnstructuredMesh<2, size_t, double>;
    sit mesh;

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
            DBGVAR(tmp_edge, mesh.GetFaces().at(tmp_edge).VertexA,mesh.GetFaces().at(tmp_edge).VertexB,mesh.GetVertices().at(mesh.GetFaces().at(tmp_edge).VertexA), mesh.GetVertices().at(mesh.GetFaces().at(tmp_edge).VertexB), mesh.CalculateEdgeMeasure(tmp_edge))
            tmp_edge = mesh.GetFaces().at(tmp_edge).GetNextBElem(i);
        } while(edge != tmp_edge);

    }

    for(size_t i = 0; i < 2; i++) {
        DBGVAR(mesh.GetCells().at(i).GetCenter(), mesh.CalculateCellMeasure(i))
        DBGVAR(mesh.GetCells().at(i).GetFlag())
    }
    DBGVAR(mesh.CalculateEdgeMeasure(2)/mesh.CalculateCellDist(0,1))
    DBGVAR(mesh.CalculateFaceMeasureOverCellDist(2))


    sit::MeshElementWrap<1> ele(&mesh, mesh.GetEdges().at(0));
    sit::MeshElementWrap<2> cell(&mesh, mesh.GetCells().at(0));

    DBGMSG("cell boundary iterator test");

    for(auto i : cell.GetSubelements()){
        DBGVAR(i)
    }


    for(auto i : mesh.GetElement<2>(1).GetSubelements()){
        DBGVAR(i)
    }

    DBGMSG("cell vertices using iterator");
    for(size_t i = 0; i < mesh.GetCells().size(); i++){
        for(auto j : mesh.GetElement<2>(i).GetSubelements()){
            auto e = mesh.GetElement<1>(j);
            DBGVAR(e.GetElement().GetVertexAIndex(), e.GetElement().GetVertexBIndex())
        }
    }

    CellsVertices<2,2,size_t, double>::run(mesh);

    DBGMSG("3D test");

    using sit3 = UnstructuredMesh<3, size_t, double, 6>;
    UnstructuredMesh<3, size_t, double, 6> mesh3;
    size_t index = 0;
DBGCHECK;
    mesh3.GetVertices().push_back({index, {0,0,0}});
    mesh3.GetVertices().push_back({1, {1,0,0}});
    mesh3.GetVertices().push_back({2, {0,1,0}});
    mesh3.GetVertices().push_back({3, {1,1,0}});
    mesh3.GetVertices().push_back({4, {0,0,1}});
    mesh3.GetVertices().push_back({5, {1,0,1}});
    mesh3.GetVertices().push_back({6, {0,1,1}});
    mesh3.GetVertices().push_back({7, {1,1,1}});
DBGCHECK;
    mesh3.GetEdges().push_back({0,0,1});
    mesh3.GetEdges().push_back({1,0,2});
    mesh3.GetEdges().push_back({2,1,3});
    mesh3.GetEdges().push_back({3,2,3});
    mesh3.GetEdges().push_back({4,0,4});
    mesh3.GetEdges().push_back({5,1,5});
    mesh3.GetEdges().push_back({6,3,4});
    mesh3.GetEdges().push_back({7,2,6});
    mesh3.GetEdges().push_back({8,4,5});
    mesh3.GetEdges().push_back({9,5,7});
    mesh3.GetEdges().push_back({10,6,7});
    mesh3.GetEdges().push_back({11,4,6});
DBGCHECK;
    mesh3.GetFaces().push_back(0);
    mesh3.GetFaces().at(0).GetSubelements().AddSubelement(0,true);
    mesh3.GetFaces().at(0).GetSubelements().AddSubelement(1,true);
    mesh3.GetFaces().at(0).GetSubelements().AddSubelement(2,true);
    mesh3.GetFaces().at(0).GetSubelements().AddSubelement(3,true);
    mesh3.GetFaces().push_back(1);
    mesh3.GetFaces().at(1).GetSubelements().AddSubelement(0,true);
    mesh3.GetFaces().at(1).GetSubelements().AddSubelement(3,true);
    mesh3.GetFaces().at(1).GetSubelements().AddSubelement(8,true);
    mesh3.GetFaces().at(1).GetSubelements().AddSubelement(4,true);
    mesh3.GetFaces().push_back(2);
    mesh3.GetFaces().at(2).GetSubelements().AddSubelement(1,true);
    mesh3.GetFaces().at(2).GetSubelements().AddSubelement(4,true);
    mesh3.GetFaces().at(2).GetSubelements().AddSubelement(11,true);
    mesh3.GetFaces().at(2).GetSubelements().AddSubelement(7,true);
    mesh3.GetFaces().push_back(3);
    mesh3.GetFaces().at(3).GetSubelements().AddSubelement(3,true);
    mesh3.GetFaces().at(3).GetSubelements().AddSubelement(6,true);
    mesh3.GetFaces().at(3).GetSubelements().AddSubelement(10,true);
    mesh3.GetFaces().at(3).GetSubelements().AddSubelement(7,true);
    mesh3.GetFaces().push_back(4);
    mesh3.GetFaces().at(4).GetSubelements().AddSubelement(2,true);
    mesh3.GetFaces().at(4).GetSubelements().AddSubelement(6,true);
    mesh3.GetFaces().at(4).GetSubelements().AddSubelement(9,true);
    mesh3.GetFaces().at(4).GetSubelements().AddSubelement(5,true);
    mesh3.GetFaces().push_back(5);
    mesh3.GetFaces().at(5).GetSubelements().AddSubelement(8,true);
    mesh3.GetFaces().at(5).GetSubelements().AddSubelement(9,true);
    mesh3.GetFaces().at(5).GetSubelements().AddSubelement(10,true);
    mesh3.GetFaces().at(5).GetSubelements().AddSubelement(11,true);
DBGCHECK;

    mesh3.GetFaces().at(0).SetNextBElem(1,0);
    mesh3.GetFaces().at(1).SetNextBElem(2,0);
    mesh3.GetFaces().at(2).SetNextBElem(3,0);
    mesh3.GetFaces().at(3).SetNextBElem(4,0);
    mesh3.GetFaces().at(4).SetNextBElem(5,0);
    mesh3.GetFaces().at(5).SetNextBElem(0,0);

    mesh3.GetCells().push_back(0);
    mesh3.GetCells().at(0).SetBoundaryElementIndex(3);
DBGCHECK;
    size_t tmp_face = mesh3.GetCells().at(0).GetBoundaryElementIndex();

    do {
        DBGVAR(tmp_face)
        for (auto& sube : mesh3.GetFaces().at(tmp_face).GetSubelements()) {
            DBGVAR(sube.index)
            if (sube.index < 12 ){
                DBGVAR(sube.index, mesh3.GetVertices().at(mesh3.GetEdges().at(sube.index).GetVertexAIndex()),mesh3.GetVertices().at(mesh3.GetEdges().at(sube.index).GetVertexBIndex()))
            }
        }

        tmp_face = mesh3.GetFaces().at(tmp_face).GetNextBElem(0);
    } while (tmp_face != mesh3.GetCells().at(0).GetBoundaryElementIndex());
//    mesh3.GetElements<0>().at(0).;


    DBGMSG("Iterator wrapper test");
    sit3::MeshElementWrap<2> elem(&mesh3, mesh3.GetFaces().at(0));
    for(auto i : elem.GetSubelements()){
        DBGVAR(i.index)
    }

    DBGMSG("3D cell vertices");
    CellsVertices<3,3,size_t, double, 6>::run(mesh3);



    DBGMSG("mesh conatiner test");
    MeshDataContainer<double, 3,2,1,0> cont(mesh3);

    //cont.GetDataDim<3>().resize(20);
    DBG(cont.GetDataPos<1>().size());

    DBG(cont.GetDataPos<0>().size());

    DBGMSG("faceCenters");
    MeshDataContainer<Vertex<3,double>, 3,2,1> centers(mesh3);
    _ComputeCenters<1,3>::compute<size_t, double, 6>(centers, mesh3);

    auto& faceCent = centers.GetDataDim<2>();
    for(auto& center : faceCent) {
        DBGVAR(center)
    }
    DBGMSG("cellCenter");
    DBGVAR(centers.GetDataDim<3>().at(0))

    //DBGVAR(centers.template GetDataDim<3>().at(0))

    //DBGVAR(centers.template GetDataDim<3>().at(0));
}

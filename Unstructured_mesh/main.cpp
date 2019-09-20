#include <iostream>
#include "../debug/debug.h"
#include "unstructuredmesh.h"
#include "mesh_functions.h"
using namespace std;




void cube(UnstructuredMesh<3, size_t, double, 6>& mesh3){
    DBGCHECK;
        mesh3.GetVertices().push_back({0, {0,0,0}});
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
        mesh3.GetEdges().push_back({6,3,7});
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
        mesh3.GetFaces().at(1).GetSubelements().AddSubelement(5,true);
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
}

void twoPrisms(UnstructuredMesh<3, size_t, double, 6>& mesh3){
    DBGCHECK;
        mesh3.GetVertices().push_back({0, {0,0,0}});
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
        mesh3.GetEdges().push_back({6,3,7});
        mesh3.GetEdges().push_back({7,2,6});
        mesh3.GetEdges().push_back({8,4,5});
        mesh3.GetEdges().push_back({9,5,7});
        mesh3.GetEdges().push_back({10,6,7});
        mesh3.GetEdges().push_back({11,4,6});
        mesh3.GetEdges().push_back({12,6,5});
        mesh3.GetEdges().push_back({13,2,1});
    DBGCHECK;
    size_t index = 0;
        mesh3.GetFaces().push_back(index);
        mesh3.GetFaces().at(index).GetSubelements().AddSubelement(0,true);
        mesh3.GetFaces().at(index).GetSubelements().AddSubelement(1,true);
        mesh3.GetFaces().at(index).GetSubelements().AddSubelement(13,true);
        mesh3.GetFaces().push_back(++index);
        mesh3.GetFaces().at(index).GetSubelements().AddSubelement(13,true);
        mesh3.GetFaces().at(index).GetSubelements().AddSubelement(2,true);
        mesh3.GetFaces().at(index).GetSubelements().AddSubelement(3,true);
        mesh3.GetFaces().push_back(++index);
        mesh3.GetFaces().at(index).GetSubelements().AddSubelement(0,true);
        mesh3.GetFaces().at(index).GetSubelements().AddSubelement(5,true);
        mesh3.GetFaces().at(index).GetSubelements().AddSubelement(8,true);
        mesh3.GetFaces().at(index).GetSubelements().AddSubelement(4,true);
        mesh3.GetFaces().push_back(++index);
        mesh3.GetFaces().at(index).GetSubelements().AddSubelement(1,true);
        mesh3.GetFaces().at(index).GetSubelements().AddSubelement(4,true);
        mesh3.GetFaces().at(index).GetSubelements().AddSubelement(11,true);
        mesh3.GetFaces().at(index).GetSubelements().AddSubelement(7,true);
        mesh3.GetFaces().push_back(++index);
        mesh3.GetFaces().at(index).GetSubelements().AddSubelement(3,true);
        mesh3.GetFaces().at(index).GetSubelements().AddSubelement(6,true);
        mesh3.GetFaces().at(index).GetSubelements().AddSubelement(10,true);
        mesh3.GetFaces().at(index).GetSubelements().AddSubelement(7,true);
        mesh3.GetFaces().push_back(++index);
        mesh3.GetFaces().at(index).GetSubelements().AddSubelement(2,true);
        mesh3.GetFaces().at(index).GetSubelements().AddSubelement(6,true);
        mesh3.GetFaces().at(index).GetSubelements().AddSubelement(9,true);
        mesh3.GetFaces().at(index).GetSubelements().AddSubelement(5,true);
        mesh3.GetFaces().push_back(++index);
        mesh3.GetFaces().at(index).GetSubelements().AddSubelement(8,true);
        mesh3.GetFaces().at(index).GetSubelements().AddSubelement(12,true);
        mesh3.GetFaces().at(index).GetSubelements().AddSubelement(11,true);
        mesh3.GetFaces().push_back(++index);
        mesh3.GetFaces().at(index).GetSubelements().AddSubelement(9,true);
        mesh3.GetFaces().at(index).GetSubelements().AddSubelement(10,true);
        mesh3.GetFaces().at(index).GetSubelements().AddSubelement(12,true);
        mesh3.GetFaces().push_back(++index);
        mesh3.GetFaces().at(index).GetSubelements().AddSubelement(12,true);
        mesh3.GetFaces().at(index).GetSubelements().AddSubelement(7,true);
        mesh3.GetFaces().at(index).GetSubelements().AddSubelement(13,true);
        mesh3.GetFaces().at(index).GetSubelements().AddSubelement(5,true);
    DBGCHECK;

        mesh3.GetFaces().at(0).SetNextBElem(2,0);
        mesh3.GetFaces().at(2).SetNextBElem(3,0);
        mesh3.GetFaces().at(3).SetNextBElem(8,0);
        mesh3.GetFaces().at(8).SetNextBElem(6,0);
        mesh3.GetFaces().at(6).SetNextBElem(0,0);
        mesh3.GetCells().push_back(0);
        mesh3.GetCells().at(0).SetBoundaryElementIndex(0);

        mesh3.GetFaces().at(1).SetNextBElem(5,1);
        mesh3.GetFaces().at(5).SetNextBElem(4,1);
        mesh3.GetFaces().at(4).SetNextBElem(8,1);
        mesh3.GetFaces().at(8).SetNextBElem(7,1);
        mesh3.GetFaces().at(7).SetNextBElem(1,1);
        mesh3.GetCells().push_back(1);
        mesh3.GetCells().at(1).SetBoundaryElementIndex(1);

    DBGCHECK;
}

void twoDeformedPrisms(UnstructuredMesh<3, size_t, double, 6>& mesh3){
    DBGCHECK;
        mesh3.GetVertices().push_back({0, {0,0,1}});
        mesh3.GetVertices().push_back({1, {1,0,0}});
        mesh3.GetVertices().push_back({2, {0,1,0}});
        mesh3.GetVertices().push_back({3, {1,1,0}});
        mesh3.GetVertices().push_back({4, {0,0,2}});
        mesh3.GetVertices().push_back({5, {1,0,1}});
        mesh3.GetVertices().push_back({6, {0,1,1}});
        mesh3.GetVertices().push_back({7, {1,1,2}});
    DBGCHECK;
        mesh3.GetEdges().push_back({0,0,1});
        mesh3.GetEdges().push_back({1,0,2});
        mesh3.GetEdges().push_back({2,1,3});
        mesh3.GetEdges().push_back({3,2,3});
        mesh3.GetEdges().push_back({4,0,4});
        mesh3.GetEdges().push_back({5,1,5});
        mesh3.GetEdges().push_back({6,3,7});
        mesh3.GetEdges().push_back({7,2,6});
        mesh3.GetEdges().push_back({8,4,5});
        mesh3.GetEdges().push_back({9,5,7});
        mesh3.GetEdges().push_back({10,6,7});
        mesh3.GetEdges().push_back({11,4,6});
        mesh3.GetEdges().push_back({12,6,5});
        mesh3.GetEdges().push_back({13,2,1});
    DBGCHECK;
    size_t index = 0;
        mesh3.GetFaces().push_back(index);
        mesh3.GetFaces().at(index).GetSubelements().AddSubelement(0,true);
        mesh3.GetFaces().at(index).GetSubelements().AddSubelement(1,true);
        mesh3.GetFaces().at(index).GetSubelements().AddSubelement(13,true);
        mesh3.GetFaces().push_back(++index);
        mesh3.GetFaces().at(index).GetSubelements().AddSubelement(13,true);
        mesh3.GetFaces().at(index).GetSubelements().AddSubelement(2,true);
        mesh3.GetFaces().at(index).GetSubelements().AddSubelement(3,true);
        mesh3.GetFaces().push_back(++index);
        mesh3.GetFaces().at(index).GetSubelements().AddSubelement(0,true);
        mesh3.GetFaces().at(index).GetSubelements().AddSubelement(5,true);
        mesh3.GetFaces().at(index).GetSubelements().AddSubelement(8,true);
        mesh3.GetFaces().at(index).GetSubelements().AddSubelement(4,true);
        mesh3.GetFaces().push_back(++index);
        mesh3.GetFaces().at(index).GetSubelements().AddSubelement(1,true);
        mesh3.GetFaces().at(index).GetSubelements().AddSubelement(4,true);
        mesh3.GetFaces().at(index).GetSubelements().AddSubelement(11,true);
        mesh3.GetFaces().at(index).GetSubelements().AddSubelement(7,true);
        mesh3.GetFaces().push_back(++index);
        mesh3.GetFaces().at(index).GetSubelements().AddSubelement(3,true);
        mesh3.GetFaces().at(index).GetSubelements().AddSubelement(6,true);
        mesh3.GetFaces().at(index).GetSubelements().AddSubelement(10,true);
        mesh3.GetFaces().at(index).GetSubelements().AddSubelement(7,true);
        mesh3.GetFaces().push_back(++index);
        mesh3.GetFaces().at(index).GetSubelements().AddSubelement(2,true);
        mesh3.GetFaces().at(index).GetSubelements().AddSubelement(6,true);
        mesh3.GetFaces().at(index).GetSubelements().AddSubelement(9,true);
        mesh3.GetFaces().at(index).GetSubelements().AddSubelement(5,true);
        mesh3.GetFaces().push_back(++index);
        mesh3.GetFaces().at(index).GetSubelements().AddSubelement(8,true);
        mesh3.GetFaces().at(index).GetSubelements().AddSubelement(12,true);
        mesh3.GetFaces().at(index).GetSubelements().AddSubelement(11,true);
        mesh3.GetFaces().push_back(++index);
        mesh3.GetFaces().at(index).GetSubelements().AddSubelement(9,true);
        mesh3.GetFaces().at(index).GetSubelements().AddSubelement(10,true);
        mesh3.GetFaces().at(index).GetSubelements().AddSubelement(12,true);
        mesh3.GetFaces().push_back(++index);
        mesh3.GetFaces().at(index).GetSubelements().AddSubelement(12,true);
        mesh3.GetFaces().at(index).GetSubelements().AddSubelement(7,true);
        mesh3.GetFaces().at(index).GetSubelements().AddSubelement(13,true);
        mesh3.GetFaces().at(index).GetSubelements().AddSubelement(5,true);
    DBGCHECK;

        mesh3.GetFaces().at(0).SetNextBElem(2,0);
        mesh3.GetFaces().at(2).SetNextBElem(3,0);
        mesh3.GetFaces().at(3).SetNextBElem(8,0);
        mesh3.GetFaces().at(8).SetNextBElem(6,0);
        mesh3.GetFaces().at(6).SetNextBElem(0,0);
        mesh3.GetCells().push_back(0);
        mesh3.GetCells().at(0).SetBoundaryElementIndex(0);

        mesh3.GetFaces().at(1).SetNextBElem(5,1);
        mesh3.GetFaces().at(5).SetNextBElem(4,1);
        mesh3.GetFaces().at(4).SetNextBElem(8,1);
        mesh3.GetFaces().at(8).SetNextBElem(7,1);
        mesh3.GetFaces().at(7).SetNextBElem(1,1);
        mesh3.GetCells().push_back(1);
        mesh3.GetCells().at(1).SetBoundaryElementIndex(1);

    DBGCHECK;
}

void testMesh2D() {
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

    auto centers = ComputeCenters(mesh);

    auto& faceCent = centers.GetDataDim<1>();
    for(auto& center : faceCent) {
        DBGVAR(center)
    }
    DBGMSG("cellCenter");
    for(sit::Cell& cell : mesh.GetCells()){
        DBGVAR(centers.GetDataDim<2>().at(cell.GetIndex()))
    }

/*
    DBGMSG("computing measures");

    auto measures = ComputeMeasures(mesh);


    for(double edgeM :measures.GetDataDim<1>()) {
        DBGVAR(edgeM)
    }

    for(double cellM :measures.GetDataDim<2>()) {
        DBGVAR(cellM)
    }
*/
}


void testMesh3D() {
    DBGMSG("3D test");

    using sit3 = UnstructuredMesh<3, size_t, double, 6>;
    UnstructuredMesh<3, size_t, double, 6> mesh3;
    twoPrisms(mesh3);
    size_t tmp_face = mesh3.GetCells().at(0).GetBoundaryElementIndex();

    do {
        DBGVAR(tmp_face)
        for (auto& sube : mesh3.GetFaces().at(tmp_face).GetSubelements()) {
            DBGVAR(sube.index)
            if (sube.index != INVALID_INDEX(size_t) ){
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
    DBGVAR(cont.GetDataPos<1>().size())

    DBGVAR(cont.GetDataPos<0>().size())

    DBGVAR(cont.GetDataDim<3>().size())


    DBGMSG("faceCenters");
    //MeshDataContainer<Vertex<3,double>, 3,2,1> centers(mesh3);


    //_ComputeCenters<1,3, 3,2,1>::compute<size_t, double, 6>(centers, mesh3);
    auto centers = ComputeCenters(mesh3);

    for(auto& face : mesh3.GetFaces()) {
        face.SetCenter(centers.template GetDataDim<2>().at(face.GetIndex()));
        DBGVAR(face.GetCenter())
    }
    DBGMSG("cellCenter");
    for(auto& cell : mesh3.GetCells()) {
        cell.SetCenter(centers.template GetDataDim<3>().at(cell.GetIndex()));
        DBGVAR(cell.GetCenter())
    }

    DBGMSG("measure computation");

    auto measures = ComputeMeasures(mesh3);
    for(double edgeM : measures.GetDataDim<1>()) {
        DBGVAR(edgeM)
    }
    for(double faceM : measures.GetDataDim<2>()) {
        DBGVAR(faceM)
    }

    for(double cellM : measures.GetDataDim<3>()) {
        DBGVAR(cellM)
    }

    //DBGVAR(centers.template GetDataDim<3>().at(0))

    //DBGVAR(centers.template GetDataDim<3>().at(0));
}



void test3DMeshDeformedPrisms() {
    UnstructuredMesh<3, size_t, double, 6> mesh3;
    twoDeformedPrisms(mesh3);

    //_ComputeCenters<1,3, 3,2,1>::compute<size_t, double, 6>(centers, mesh3);
    auto centers = ComputeCenters(mesh3);

    for(auto& face : mesh3.GetFaces()) {
        face.SetCenter(centers[face]);
        DBGVAR(face.GetCenter())
    }
    DBGMSG("cellCenter");
    for(auto& cell : mesh3.GetCells()) {
        cell.SetCenter(centers[cell]);
        DBGVAR(cell.GetCenter())
    }

    DBGMSG("measure computation");

    auto measures = ComputeMeasures(mesh3);
    for(double edgeM : measures.GetDataDim<1>()) {
        DBGVAR(edgeM)
    }
    for(double faceM : measures.GetDataDim<2>()) {
        DBGVAR(faceM)
    }

    for(double cellM : measures.GetDataDim<3>()) {
        DBGVAR(cellM)
    }
}








void testMeshDataContainer() {
    UnstructuredMesh<3, size_t, double, 6> mesh3;
    twoDeformedPrisms(mesh3);


    MeshDataContainer<std::tuple<int, double, char, double>, 3,2,0> container(mesh3);



    for(auto& c : container.GetDataDim<0>()) {
        c=42;
    }

    for(auto& v : mesh3.GetVertices()){
        DBGVAR(container.at(v))
    }
}

template <unsigned int ... Is>
class ClassA
 {
   public:
      ClassA (std::integer_sequence<unsigned int,Is...>)
      {DBGVAR(sizeof... (Is), std::get<0>(std::array<size_t, sizeof...(Is)>{Is...})) }

      static void fun (std::index_sequence<Is...>)
      {DBGVAR(sizeof... (Is), std::get<0>(std::array<size_t, sizeof...(Is)>{Is...})) }

 };


//template <typename ... t> class ClassB{};

template<typename  Tuple, unsigned int ... Is>
class ClassB
 {
   public:
      ClassB (const std::integer_sequence<unsigned int, Is...>, Tuple t)
      {std::tuple_element_t<0,Tuple> typ = 0;
          DBGVAR(sizeof... (Is), std::get<0>(std::array<size_t, sizeof...(Is)>{Is...}), std::get<0>(t), typ) }

 };


template <typename ... t> class ClassC{};

template<unsigned int ... Is, typename ...Types>
class ClassC<std::integer_sequence<unsigned int,Is...>, std::tuple<Types...>>
 {
   public:
      ClassC (const std::integer_sequence<unsigned int, Is...>, std::tuple<Types...>)
      {std::tuple_element_t<0,std::tuple<Types...>> typ = 0;
          DBGVAR(sizeof... (Is), std::get<0>(std::array<size_t, sizeof...(Is)>{Is...}), typ) }

      ClassC () {
          std::tuple_element_t<0,std::tuple<Types...>> typ = 42.15;
          std::tuple_element_t<1,std::tuple<Types...>> typ2 = 42.15;
          DBGVAR(sizeof... (Is), std::get<0>(std::array<size_t, sizeof...(Is)>{Is...}), typ, typ2)
      }
 };



void testTemplate() {
    ClassA n(std::make_integer_sequence<unsigned int, 3>{});
    UnstructuredMesh<3,size_t, double,6> mesh3;
    //MeshDataContainer<Vertex<3, double>, 0,1,2> centers2(mesh3,std::make_integer_sequence<unsigned int, 3>{}, Vertex<3, double>{});
    //ComputeCenters(mesh3);

    ClassA p(make_custom_integer_sequence<unsigned int, 10, 0, -2>{});
    std::tuple<double, char> t{};
    t={1,2};
    ClassB u(make_custom_integer_sequence<unsigned int, 2, 0, -2>{}, t);

    ClassC<std::integer_sequence<unsigned int, 2,0>, std::tuple<double, char>> c(make_custom_integer_sequence<unsigned int, 2, 0, -2>{}, std::tuple<double, char>{});
    ClassC<std::integer_sequence<unsigned int, 2,0>, std::tuple<double, char>> cc;
    ClassC<std::integer_sequence<unsigned int, 2,0>, decltype(std::make_tuple(1.0, 'a'))> ccc;
}





int main()
{
    //testMesh2D();
    //testMesh3D();
    //test3DMeshDeformedPrisms();
    testMeshDataContainer();
    //testTemplate();
}

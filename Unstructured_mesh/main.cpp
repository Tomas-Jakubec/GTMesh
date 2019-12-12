#include <iostream>
#include <string>
//#define UNDEBUG
#define CONSOLE_COLOURED_OUTPUT
#include "../src/Debug/Debug.h"
#include "../src/UnstructuredMesh/UnstructuredMesh.h"
#include "../src/UnstructuredMesh/MeshFunctions/MeshFunctions.h"
#include "../src/UnstructuredMesh/MeshIO/MeshReader/VTKMeshReader.h"
#include "../src/UnstructuredMesh/MeshIO/MeshWriter/VTKMeshWriter.h"
#include "../src/UnstructuredMesh/MeshDataContainer/MeshDataIO/VTKMeshDataWriter.h"
#include "../src/UnstructuredMesh/MeshDataContainer/MeshDataIO/VTKMeshDataReader.h"
#include "../src/UnstructuredMesh/MeshIO/MeshReader/FPMAMeshReader.h"
#include "../src/UnstructuredMesh/MeshIO/MeshWriter/FPMAMeshWriter.h"

#include "../src/Traits/MemberApproach/MemberApproach.h"
#include <fstream>
#include <list>
using namespace std;




void cube(UnstructuredMesh<3, size_t, double, 6>& mesh3){
    DBGCHECK;
        mesh3.getVertices().push_back({0, {0,0,0}});
        mesh3.getVertices().push_back({1, {1,0,0}});
        mesh3.getVertices().push_back({2, {0,1,0}});
        mesh3.getVertices().push_back({3, {1,1,0}});
        mesh3.getVertices().push_back({4, {0,0,1}});
        mesh3.getVertices().push_back({5, {1,0,1}});
        mesh3.getVertices().push_back({6, {0,1,1}});
        mesh3.getVertices().push_back({7, {1,1,1}});
    DBGCHECK;
        mesh3.getEdges().push_back({0,0,1});
        mesh3.getEdges().push_back({1,0,2});
        mesh3.getEdges().push_back({2,1,3});
        mesh3.getEdges().push_back({3,2,3});
        mesh3.getEdges().push_back({4,0,4});
        mesh3.getEdges().push_back({5,1,5});
        mesh3.getEdges().push_back({6,3,7});
        mesh3.getEdges().push_back({7,2,6});
        mesh3.getEdges().push_back({8,4,5});
        mesh3.getEdges().push_back({9,5,7});
        mesh3.getEdges().push_back({10,6,7});
        mesh3.getEdges().push_back({11,4,6});
    DBGCHECK;
        mesh3.getFaces().push_back(0);
        mesh3.getFaces().at(0).getSubelements().addSubelement(0);
        mesh3.getFaces().at(0).getSubelements().addSubelement(1);
        mesh3.getFaces().at(0).getSubelements().addSubelement(2);
        mesh3.getFaces().at(0).getSubelements().addSubelement(3);
        mesh3.getFaces().push_back(1);
        mesh3.getFaces().at(1).getSubelements().addSubelement(0);
        mesh3.getFaces().at(1).getSubelements().addSubelement(5);
        mesh3.getFaces().at(1).getSubelements().addSubelement(8);
        mesh3.getFaces().at(1).getSubelements().addSubelement(4);
        mesh3.getFaces().push_back(2);
        mesh3.getFaces().at(2).getSubelements().addSubelement(1);
        mesh3.getFaces().at(2).getSubelements().addSubelement(4);
        mesh3.getFaces().at(2).getSubelements().addSubelement(11);
        mesh3.getFaces().at(2).getSubelements().addSubelement(7);
        mesh3.getFaces().push_back(3);
        mesh3.getFaces().at(3).getSubelements().addSubelement(3);
        mesh3.getFaces().at(3).getSubelements().addSubelement(6);
        mesh3.getFaces().at(3).getSubelements().addSubelement(10);
        mesh3.getFaces().at(3).getSubelements().addSubelement(7);
        mesh3.getFaces().push_back(4);
        mesh3.getFaces().at(4).getSubelements().addSubelement(2);
        mesh3.getFaces().at(4).getSubelements().addSubelement(6);
        mesh3.getFaces().at(4).getSubelements().addSubelement(9);
        mesh3.getFaces().at(4).getSubelements().addSubelement(5);
        mesh3.getFaces().push_back(5);
        mesh3.getFaces().at(5).getSubelements().addSubelement(8);
        mesh3.getFaces().at(5).getSubelements().addSubelement(9);
        mesh3.getFaces().at(5).getSubelements().addSubelement(10);
        mesh3.getFaces().at(5).getSubelements().addSubelement(11);
    DBGCHECK;

        mesh3.getFaces().at(0).setNextBElem(1,0);
        mesh3.getFaces().at(1).setNextBElem(2,0);
        mesh3.getFaces().at(2).setNextBElem(3,0);
        mesh3.getFaces().at(3).setNextBElem(4,0);
        mesh3.getFaces().at(4).setNextBElem(5,0);
        mesh3.getFaces().at(5).setNextBElem(0,0);

        mesh3.getCells().push_back(0);
        mesh3.getCells().at(0).setBoundaryElementIndex(3);
        mesh3.updateSignature();
    DBGCHECK;
}

void twoPrisms(UnstructuredMesh<3, size_t, double, 6>& mesh3){
    DBGCHECK;
        mesh3.getVertices().push_back({0, {0,0,0}});
        mesh3.getVertices().push_back({1, {1,0,0}});
        mesh3.getVertices().push_back({2, {0,1,0}});
        mesh3.getVertices().push_back({3, {1,1,0}});
        mesh3.getVertices().push_back({4, {0,0,1}});
        mesh3.getVertices().push_back({5, {1,0,1}});
        mesh3.getVertices().push_back({6, {0,1,1}});
        mesh3.getVertices().push_back({7, {1,1,1}});
    DBGCHECK;
        mesh3.getEdges().push_back({0,0,1});
        mesh3.getEdges().push_back({1,0,2});
        mesh3.getEdges().push_back({2,1,3});
        mesh3.getEdges().push_back({3,2,3});
        mesh3.getEdges().push_back({4,0,4});
        mesh3.getEdges().push_back({5,1,5});
        mesh3.getEdges().push_back({6,3,7});
        mesh3.getEdges().push_back({7,2,6});
        mesh3.getEdges().push_back({8,4,5});
        mesh3.getEdges().push_back({9,5,7});
        mesh3.getEdges().push_back({10,6,7});
        mesh3.getEdges().push_back({11,4,6});
        mesh3.getEdges().push_back({12,6,5});
        mesh3.getEdges().push_back({13,2,1});
    DBGCHECK;
    size_t index = 0;
        mesh3.getFaces().push_back(index);
        mesh3.getFaces().at(index).getSubelements().addSubelement(0);
        mesh3.getFaces().at(index).getSubelements().addSubelement(1);
        mesh3.getFaces().at(index).getSubelements().addSubelement(13);
        mesh3.getFaces().push_back(++index);
        mesh3.getFaces().at(index).getSubelements().addSubelement(13);
        mesh3.getFaces().at(index).getSubelements().addSubelement(2);
        mesh3.getFaces().at(index).getSubelements().addSubelement(3);
        mesh3.getFaces().push_back(++index);
        mesh3.getFaces().at(index).getSubelements().addSubelement(0);
        mesh3.getFaces().at(index).getSubelements().addSubelement(5);
        mesh3.getFaces().at(index).getSubelements().addSubelement(8);
        mesh3.getFaces().at(index).getSubelements().addSubelement(4);
        mesh3.getFaces().push_back(++index);
        mesh3.getFaces().at(index).getSubelements().addSubelement(1);
        mesh3.getFaces().at(index).getSubelements().addSubelement(4);
        mesh3.getFaces().at(index).getSubelements().addSubelement(11);
        mesh3.getFaces().at(index).getSubelements().addSubelement(7);
        mesh3.getFaces().push_back(++index);
        mesh3.getFaces().at(index).getSubelements().addSubelement(3);
        mesh3.getFaces().at(index).getSubelements().addSubelement(6);
        mesh3.getFaces().at(index).getSubelements().addSubelement(10);
        mesh3.getFaces().at(index).getSubelements().addSubelement(7);
        mesh3.getFaces().push_back(++index);
        mesh3.getFaces().at(index).getSubelements().addSubelement(2);
        mesh3.getFaces().at(index).getSubelements().addSubelement(6);
        mesh3.getFaces().at(index).getSubelements().addSubelement(9);
        mesh3.getFaces().at(index).getSubelements().addSubelement(5);
        mesh3.getFaces().push_back(++index);
        mesh3.getFaces().at(index).getSubelements().addSubelement(8);
        mesh3.getFaces().at(index).getSubelements().addSubelement(12);
        mesh3.getFaces().at(index).getSubelements().addSubelement(11);
        mesh3.getFaces().push_back(++index);
        mesh3.getFaces().at(index).getSubelements().addSubelement(9);
        mesh3.getFaces().at(index).getSubelements().addSubelement(10);
        mesh3.getFaces().at(index).getSubelements().addSubelement(12);
        mesh3.getFaces().push_back(++index);
        mesh3.getFaces().at(index).getSubelements().addSubelement(12);
        mesh3.getFaces().at(index).getSubelements().addSubelement(7);
        mesh3.getFaces().at(index).getSubelements().addSubelement(13);
        mesh3.getFaces().at(index).getSubelements().addSubelement(5);
    DBGCHECK;

        mesh3.getFaces().at(0).setNextBElem(2,0);
        mesh3.getFaces().at(2).setNextBElem(3,0);
        mesh3.getFaces().at(3).setNextBElem(8,0);
        mesh3.getFaces().at(8).setNextBElem(6,0);
        mesh3.getFaces().at(6).setNextBElem(0,0);
        mesh3.getCells().push_back(0);
        mesh3.getCells().at(0).setBoundaryElementIndex(0);

        mesh3.getFaces().at(1).setNextBElem(5,1);
        mesh3.getFaces().at(5).setNextBElem(4,1);
        mesh3.getFaces().at(4).setNextBElem(8,1);
        mesh3.getFaces().at(8).setNextBElem(7,1);
        mesh3.getFaces().at(7).setNextBElem(1,1);
        mesh3.getCells().push_back(1);
        mesh3.getCells().at(1).setBoundaryElementIndex(1);

        mesh3.updateSignature();
    DBGCHECK;
}

void twoDeformedPrisms(UnstructuredMesh<3, size_t, double, 6>& mesh3){
    DBGCHECK;
        mesh3.getVertices().push_back({0, {0,0,1}});
        mesh3.getVertices().push_back({1, {1,0,0}});
        mesh3.getVertices().push_back({2, {0,1,0}});
        mesh3.getVertices().push_back({3, {1,1,0}});
        mesh3.getVertices().push_back({4, {0,0,2}});
        mesh3.getVertices().push_back({5, {1,0,1}});
        mesh3.getVertices().push_back({6, {0,1,1}});
        mesh3.getVertices().push_back({7, {1,1,2}});
    DBGCHECK;
        mesh3.getEdges().push_back({0,0,1});
        mesh3.getEdges().push_back({1,0,2});
        mesh3.getEdges().push_back({2,1,3});
        mesh3.getEdges().push_back({3,2,3});
        mesh3.getEdges().push_back({4,0,4});
        mesh3.getEdges().push_back({5,1,5});
        mesh3.getEdges().push_back({6,3,7});
        mesh3.getEdges().push_back({7,2,6});
        mesh3.getEdges().push_back({8,4,5});
        mesh3.getEdges().push_back({9,5,7});
        mesh3.getEdges().push_back({10,6,7});
        mesh3.getEdges().push_back({11,4,6});
        mesh3.getEdges().push_back({12,6,5});
        mesh3.getEdges().push_back({13,2,1});
    DBGCHECK;
    size_t index = 0;
        mesh3.getFaces().push_back(index);
        mesh3.getFaces().at(index).getSubelements().addSubelement(0);
        mesh3.getFaces().at(index).getSubelements().addSubelement(1);
        mesh3.getFaces().at(index).getSubelements().addSubelement(13);
        mesh3.getFaces().push_back(++index);
        mesh3.getFaces().at(index).getSubelements().addSubelement(13);
        mesh3.getFaces().at(index).getSubelements().addSubelement(2);
        mesh3.getFaces().at(index).getSubelements().addSubelement(3);
        mesh3.getFaces().push_back(++index);
        mesh3.getFaces().at(index).getSubelements().addSubelement(0);
        mesh3.getFaces().at(index).getSubelements().addSubelement(5);
        mesh3.getFaces().at(index).getSubelements().addSubelement(8);
        mesh3.getFaces().at(index).getSubelements().addSubelement(4);
        mesh3.getFaces().push_back(++index);
        mesh3.getFaces().at(index).getSubelements().addSubelement(1);
        mesh3.getFaces().at(index).getSubelements().addSubelement(4);
        mesh3.getFaces().at(index).getSubelements().addSubelement(11);
        mesh3.getFaces().at(index).getSubelements().addSubelement(7);
        mesh3.getFaces().push_back(++index);
        mesh3.getFaces().at(index).getSubelements().addSubelement(3);
        mesh3.getFaces().at(index).getSubelements().addSubelement(6);
        mesh3.getFaces().at(index).getSubelements().addSubelement(10);
        mesh3.getFaces().at(index).getSubelements().addSubelement(7);
        mesh3.getFaces().push_back(++index);
        mesh3.getFaces().at(index).getSubelements().addSubelement(2);
        mesh3.getFaces().at(index).getSubelements().addSubelement(6);
        mesh3.getFaces().at(index).getSubelements().addSubelement(9);
        mesh3.getFaces().at(index).getSubelements().addSubelement(5);
        mesh3.getFaces().push_back(++index);
        mesh3.getFaces().at(index).getSubelements().addSubelement(8);
        mesh3.getFaces().at(index).getSubelements().addSubelement(12);
        mesh3.getFaces().at(index).getSubelements().addSubelement(11);
        mesh3.getFaces().push_back(++index);
        mesh3.getFaces().at(index).getSubelements().addSubelement(9);
        mesh3.getFaces().at(index).getSubelements().addSubelement(10);
        mesh3.getFaces().at(index).getSubelements().addSubelement(12);
        mesh3.getFaces().push_back(++index);
        mesh3.getFaces().at(index).getSubelements().addSubelement(12);
        mesh3.getFaces().at(index).getSubelements().addSubelement(7);
        mesh3.getFaces().at(index).getSubelements().addSubelement(13);
        mesh3.getFaces().at(index).getSubelements().addSubelement(5);
    DBGCHECK;

        mesh3.getFaces().at(0).setNextBElem(2,0);
        mesh3.getFaces().at(2).setNextBElem(3,0);
        mesh3.getFaces().at(3).setNextBElem(8,0);
        mesh3.getFaces().at(8).setNextBElem(6,0);
        mesh3.getFaces().at(6).setNextBElem(0,0);
        mesh3.getCells().push_back(0);
        mesh3.getCells().at(0).setBoundaryElementIndex(0);

        mesh3.getFaces().at(1).setNextBElem(5,1);
        mesh3.getFaces().at(5).setNextBElem(4,1);
        mesh3.getFaces().at(4).setNextBElem(8,1);
        mesh3.getFaces().at(8).setNextBElem(7,1);
        mesh3.getFaces().at(7).setNextBElem(1,1);
        mesh3.getCells().push_back(1);
        mesh3.getCells().at(1).setBoundaryElementIndex(1);

        mesh3.updateSignature();
    DBGCHECK;
}

void testMesh2D() {
    using sit = UnstructuredMesh<2, size_t, double>;
    sit mesh;

    mesh.getVertices().resize(4);

    mesh.getVertices().at(0).setIndex(0);
    mesh.getVertices().at(0) = {0,0};
    mesh.getVertices().at(1).setIndex(1);
    mesh.getVertices().at(1) = {0,1};
    mesh.getVertices().at(2).setIndex(2);
    mesh.getVertices().at(2) = {1,0};
    mesh.getVertices().at(3).setIndex(3);
    mesh.getVertices().at(3) = {1,1};

    mesh.getFaces().resize(5);
    mesh.getFaces().at(0).vertexAIndex = 1;
    mesh.getFaces().at(0).vertexBIndex = 0;
    mesh.getFaces().at(1).vertexAIndex = 0;
    mesh.getFaces().at(1).vertexBIndex = 2;
    mesh.getFaces().at(2).vertexAIndex = 2;
    mesh.getFaces().at(2).vertexBIndex = 1;
    mesh.getFaces().at(3).vertexAIndex = 3;
    mesh.getFaces().at(3).vertexBIndex = 1;
    mesh.getFaces().at(4).vertexAIndex = 2;
    mesh.getFaces().at(4).vertexBIndex = 3;
    for(size_t i = 0; i < 5; i++)
        mesh.getFaces().at(i).setIndex(i);

    mesh.getFaces().at(0).setNextBElem(1,0);
    mesh.getFaces().at(1).setNextBElem(2,0);
    mesh.getFaces().at(2).setNextBElem(0,0);

    mesh.getFaces().at(2).setNextBElem(3,1);
    mesh.getFaces().at(3).setNextBElem(4,1);
    mesh.getFaces().at(4).setNextBElem(2,1);


    mesh.getCells().resize(2);
    mesh.getCells().at(0).setBoundaryElementIndex(0);
    mesh.getCells().at(0).setIndex(0);
    mesh.getCells().at(1).setBoundaryElementIndex(2);
    mesh.getCells().at(1).setIndex(1);

    mesh.updateSignature();

    DBGVAR(mesh.getSignature());
    mesh.initializeCenters();



    sit::MeshElementWrap<1> ele(&mesh, mesh.getEdges().at(0));
    sit::MeshElementWrap<2> cell(&mesh, mesh.getCells().at(0));

    DBGMSG("cell boundary iterator test");

    for(auto i : cell.getSubelements()){
        DBGVAR(i);
    }


    for(auto i : mesh.getElement<2>(1).getSubelements()){
        DBGVAR(i);
    }

    DBGMSG("cell vertices using iterator");
    for(size_t i = 0; i < mesh.getCells().size(); i++){
        for(auto j : mesh.getElement<2>(i).getSubelements()){
            auto e = mesh.getElement<1>(j);
            DBGVAR(e.getElement().getVertexAIndex(), e.getElement().getVertexBIndex());
        }
    }

    MeshDataContainer<std::set<size_t>, 2> vertices(mesh);

    DBGMSG("vertices of cells in 2D");
    for (auto& cell : mesh.getCells()){
        std::set<size_t>& _set = vertices.at(cell);
        DBGVAR(cell.getIndex());
        for (size_t index: _set){
            DBGVAR(index);
        }
    }



    auto centers = ComputeCenters<DEFAULT>(mesh);

    auto& faceCent = centers.getDataByDim<1>();
    for(auto& center : faceCent) {
        DBGVAR(center);
    }
    DBGMSG("cellCenter");
    for(sit::Cell& cell : mesh.getCells()){
        DBGVAR(centers.getDataByDim<2>().at(cell.getIndex()));
    }

    DBGMSG("computing measures");

    auto measures = mesh.computeElementMeasures();


    for(double edgeM :measures.getDataByDim<1>()) {
        DBGVAR(edgeM);
    }

    for(double cellM :measures.getDataByDim<2>()) {
        DBGVAR(cellM);
    }



    DBGMSG("2D normals test");

    auto normals = ComputeFaceNormals<DEFAULT>(mesh);
    for(auto& edge : mesh.getEdges()){
        DBGVAR(edge.getIndex(),normals.at(edge));
    }

    DBGMSG("2D cells distances");

    auto distances = ComputeCellsDistance(mesh);
    for(auto& edge : mesh.getEdges()){
        DBGVAR(edge.getIndex(),distances.at(edge));
    }

}



void testMesh2DLoadAndWrite(){
    using Mesh = UnstructuredMesh<2, size_t, double>;
    Mesh mesh;
    DBGMSG("load from vtk file test");
    VTKMeshReader<2> reader;
    ifstream ifst("Test_obdelnik.vtk");
    DBGVAR(bool(ifst));
    reader.loadFromStream(ifst, mesh);

    DBGVAR(mesh.getVertices().size(), mesh.getVertices().at(4),mesh.getCells().size());


    DBGMSG("mesh apply test");
    MeshRun<2, 2, 0, 2,false, true>::run(mesh,size_t(4), size_t(4), [](size_t ori, size_t i){
        DBGVAR(ori,i);
    });

    mesh.initializeCenters();
    auto normals = mesh.computeFaceNormals();
    auto measures = mesh.computeElementMeasures();
    DBGVAR(normals.getDataByPos<0>().at(0), measures.getDataByDim<2>().at(100), measures.getDataByDim<1>().at(100), mesh.getCells().at(100).getCenter());

    VTKMeshWriter<2,size_t, double> writer;
    ofstream ofst("test_mesh_export.vtk");
    writer.writeHeader(ofst, "superData");
    writer.writeToStream(ofst, mesh, reader.getCellTypes());
}

void testMesh3D() {
    DBGMSG("3D test");

    using sit3 = UnstructuredMesh<3, size_t, double, 6>;
    UnstructuredMesh<3, size_t, double, 6> mesh3;
    twoPrisms(mesh3);
    mesh3.setupBoundaryCells();
    size_t tmp_face = mesh3.getCells().at(0).getBoundaryElementIndex();

    do {
        DBGVAR(tmp_face);
        for (auto& sube : mesh3.getFaces().at(tmp_face).getSubelements()) {
            DBGVAR(sube.index);
            if (sube.index != INVALID_INDEX(size_t) ){
                DBGVAR(sube.index, mesh3.getVertices().at(mesh3.getEdges().at(sube.index).getVertexAIndex()),mesh3.getVertices().at(mesh3.getEdges().at(sube.index).getVertexBIndex()));
            }
        }

        tmp_face = mesh3.getFaces().at(tmp_face).getNextBElem(0);
    } while (tmp_face != mesh3.getCells().at(0).getBoundaryElementIndex());
//    mesh3.getElements<0>().at(0).;


    DBGMSG("Iterator wrapper test");
    sit3::MeshElementWrap<2> elem(&mesh3, mesh3.getFaces().at(0));
    for(auto i : elem.getSubelements()){
        DBGVAR(i.index);
    }


    DBGVAR(mesh3.getSignature());

    DBGMSG("mesh conatiner test");
    MeshDataContainer<double, 3,2,1,0> cont(mesh3);

    MakeMeshDataContainer_t<double, make_custom_integer_sequence_t<unsigned int, 3,0,-1>> cont1(mesh3);

    cont = cont1;

    //cont.getDataByDim<3>().resize(20);
    DBGVAR(cont.getDataByPos<1>().size(), cont.getDataByPos<1>().getMappedDimension());

    DBGVAR(cont.getDataByPos<0>().size());

    DBGVAR(cont.getDataByDim<3>().size());


    DBGMSG("faceCenters");
    //MeshDataContainer<Vertex<3,double>, 3,2,1> centers(mesh3);


    //_ComputeCenters<1,3, 3,2,1>::compute<size_t, double, 6>(centers, mesh3);
    auto centers = ComputeCenters<DEFAULT>(mesh3);



    for(auto& face : mesh3.getFaces()) {
        face.setCenter(centers.template getDataByDim<2>().at(face.getIndex()));
        DBGVAR(face.getCenter());
    }
    DBGMSG("cellCenter - default method");
    for(auto& cell : mesh3.getCells()) {
        cell.setCenter(centers.template getDataByDim<3>().at(cell.getIndex()));
        DBGVAR(cell.getCenter());
    }

    DBGMSG("centers - tessellated faces");
    auto centers1 = ComputeCenters<TESSELLATED>(mesh3);

    for(auto& face : mesh3.getFaces()) {
        face.setCenter(centers1.template getDataByDim<2>().at(face.getIndex()));
        DBGVAR(face.getCenter());
    }
    DBGMSG("cellCenter");
    for(auto& cell : mesh3.getCells()) {
        cell.setCenter(centers1.template getDataByDim<3>().at(cell.getIndex()));
        DBGVAR(cell.getCenter());
    }

    DBGMSG("measure computation");

DBGMSG("tessellated cell volume");
    auto measures1 = ComputeMeasures<TESSELLATED>(mesh3);

    DBGVAR(measures1.getDataByDim<3>());


    auto measures = ComputeMeasures<DEFAULT>(mesh3);
    for(double edgeM : measures.getDataByDim<1>()) {
        DBGVAR(edgeM);
    }
    for(double faceM : measures.getDataByDim<2>()) {
        DBGVAR(faceM);
    }

    for(double cellM : measures.getDataByDim<3>()) {
        DBGVAR(cellM);
    }


    DBGMSG("3D normals test");

    auto normals = mesh3.computeFaceNormals();
    for(auto& face : mesh3.getFaces()){
        DBGVAR(face.getIndex(),normals.at(face));
    }

    DBGMSG("3D normals test");

    auto normalsTess = mesh3.computeFaceNormals<TESSELLATED>();
    for(auto& face : mesh3.getFaces()){
        DBGVAR(face.getIndex(),normalsTess.at(face));
    }

    DBGMSG("mesh apply test");
    MeshApply<3, 2>::apply(mesh3, [](size_t ori, size_t i){
        DBGVAR(ori,i);
    });
    DBGMSG("mesh apply test");
    MeshApply<2, 3>::apply(mesh3,[](size_t ori, size_t i){
        DBGVAR(ori,i);
    });

    DBGMSG("3D edge orientation");
    MeshApply<2, 1>::apply(mesh3,[&mesh3](size_t faceIndex, size_t edgeIndex){
        size_t iA = mesh3.getEdges().at(edgeIndex).getVertexAIndex(), iB =mesh3.getEdges().at(edgeIndex).getVertexBIndex();
        DBGVAR(faceIndex,
               edgeIndex,
               iA, iB,
               edgeIsLeft(mesh3,faceIndex, edgeIndex));
    });


    auto orientation =  edgesOrientation(mesh3);
    for(auto & face : mesh3.getFaces()){
        DBGVAR(face.getIndex(),orientation[face]);
    }


    DBGMSG("connection test");
    auto con = MeshConnections<3,0>::connections(mesh3);
    for (auto& cell : mesh3.getCells()){
            DBGVAR(cell.getIndex(), con[cell]);
    }

    DBGMSG("connections original order");
    auto conOrig = MeshConnections<3,0,Order::ORDER_ORIGINAL>::connections(mesh3);
    for (auto& cell : mesh3.getCells()){
            DBGVAR(cell.getIndex(), conOrig[cell]);
    }

    DBGMSG("connection test oposite");
    auto con1 = MeshConnections<0,3>::connections(mesh3);
    for (auto& vert : mesh3.getVertices()){
        DBGVAR(vert.getIndex(), con1[vert]);
    }

    DBGMSG("face to vertex colouring");
    auto colours = ColorMesh<2,0>::color(mesh3);
    for (auto& face : mesh3.getFaces()){
        DBGVAR(face.getIndex(), colours.at(face));
    }

    DBGMSG("vertex to face colouring");
    auto colours1 = ColorMesh<0,2>::color(mesh3);
    for (auto& vert : mesh3.getVertices()){
        DBGVAR(vert.getIndex(), colours1.at(vert));
    }

    MeshDataContainer<MeshNativeType<3>::ElementType,3> types(mesh3, MeshNativeType<3>::ElementType::WEDGE);

    VTKMeshWriter<3, size_t, double> writer;
    writer.indexMesh(mesh3, types);

    DBGVAR(writer.cellVert.getDataByPos<0>());
    ofstream out3D("3D_test_mesh_two_prisms.vtk");
    writer.writeHeader(out3D, "test data");
    writer.writeToStream(out3D, mesh3, types);


    MeshDataContainer<MeshNativeType<3>::ElementType,3> types1(mesh3);

    types1.getDataByPos<0>().at(0) = MeshNativeType<3>::ElementType::WEDGE;

    types1.getDataByPos<0>().at(1) = MeshNativeType<3>::ElementType::POLYHEDRON;

    VTKMeshWriter<3, size_t, double> writer1;
    ofstream out3D1("3D_test_mesh_two_prisms_split.vtk");
    writer1.writeHeader(out3D1, "test data");
    writer1.writeToStream(out3D1, mesh3, types1);
    DBGVAR(writer1.backwardCellIndexMapping);

}


struct colourData {
    unsigned int colour;
    Vector<3, double> firstEdgeNormal;
};

MAKE_ATTRIBUTE_TRAIT(colourData, colour, firstEdgeNormal);

void testMeshRefine() {
    UnstructuredMesh<3, size_t, double, 6> mesh;
    twoPrisms(mesh);
    mesh.initializeCenters();

    MeshDataContainer<MeshNativeType<3>::ElementType,3> types(mesh, MeshNativeType<3>::WEDGE);
    VTKMeshWriter<3, size_t, double> writer;
    ofstream out3D;
    out3D.open("mesh_refine_0.vtk");
    writer.writeHeader(out3D, "test data");
    writer.writeToStream(out3D, mesh, types);

    auto colours = MeshColoring<3,0>::color(mesh);

    out3D << "CELL_DATA " << mesh.getCells().size() << endl;
    out3D << "SCALARS cell_wrt_vertex_colour double 1\nLOOKUP_TABLE default" << endl;
    for (auto colour : colours.getDataByPos<0>()) {
        out3D << colour << ' ';
    }
    out3D.close();

    MeshDataContainer<MeshNativeType<3>::ElementType,3> types1(mesh, MeshNativeType<3>::POLYHEDRON);

    VTKMeshWriter<3, size_t, double> writer1;
    out3D.open("mesh_refine_1.vtk");
    writer1.writeHeader(out3D, "test data");
    writer1.writeToStream(out3D, mesh, types1);
    auto colours1 = MeshColoring<3,0>::color(mesh);

    MeshDataContainer<colourData, 3> cd(mesh);
    auto normals = mesh.computeFaceNormals();

    for(auto& cell : mesh.getCells()){
        cd.at(cell).colour = colours1.at(cell);
        cd.at(cell).firstEdgeNormal = normals.getDataByDim<2>().at(mesh.getFaces().at(cell.getBoundaryElementIndex()).getNextBElem(cell.getIndex()));
    }
    DBGVAR(cd.getDataByDim<3>());

    VTKMeshDataWriter<3> dataWriter;

    dataWriter.writeToStream(out3D, cd, writer1);

    //out3D << "CELL_DATA " << writer1.cellVert.getDataByPos<0>().size() << endl;
    out3D << "SCALARS cell_wrt_vertex_colour double 1\nLOOKUP_TABLE default" << endl;
    size_t realIndex = 0;
    for (size_t i = 0; i < writer1.cellVert.getDataByPos<0>().size(); i++) {
        auto iterator = writer1.backwardCellIndexMapping.find(i);
        if (iterator == writer1.backwardCellIndexMapping.end()){
            out3D << colours1.getDataByPos<0>().at(realIndex) << ' ';
            realIndex++;
        } else {
            out3D << colours1.getDataByPos<0>().at(iterator->second) << ' ';
            realIndex = iterator->second;
        }
    }
    out3D.close();


    ifstream in3D;
    in3D.open("mesh_refine_1.vtk");
    VTKMeshReader<3> reader;
    reader.loadFromStream(in3D, mesh);

    MeshDataContainer<colourData, 3> cd1(mesh);
    VTKMeshDataReader<3, size_t>::readData(in3D, cd1.getDataByPos<0>());

    DBGVAR(cd1.getDataByDim<3>());
    in3D.close();

    mesh.initializeCenters();

    MeshDataContainer<MeshNativeType<3>::ElementType,3> types2(mesh, MeshNativeType<3>::POLYHEDRON);
    out3D.open("mesh_refine_2.vtk");
    writer1.writeHeader(out3D, "test data");
    writer1.writeToStream(out3D, mesh, types2);

    auto colours2 = MeshColoring<3,0>::color(mesh);

    out3D << "CELL_DATA " << writer1.cellVert.getDataByPos<0>().size() << endl;
    out3D << "SCALARS cell_wrt_vertex_colour double 1\nLOOKUP_TABLE default" << endl;
    realIndex = 0;
    for (size_t i = 0; i < writer1.cellVert.getDataByPos<0>().size(); i++) {
        auto iterator = writer1.backwardCellIndexMapping.find(i);
        if (iterator == writer1.backwardCellIndexMapping.end()){
            out3D << colours2.getDataByPos<0>().at(realIndex) << ' ';
            realIndex++;
        } else {
            out3D << colours2.getDataByPos<0>().at(iterator->second) << ' ';
            realIndex = iterator->second;
        }
    }
    out3D.close();


    in3D.open("mesh_refine_2.vtk");
    reader.loadFromStream(in3D, mesh);
    in3D.close();

    mesh.initializeCenters();

    MeshDataContainer<MeshNativeType<3>::ElementType,3> types3(mesh, MeshNativeType<3>::POLYHEDRON);
    out3D.open("mesh_refine_3.vtk");
    writer1.writeHeader(out3D, "test data");
    writer1.writeToStream(out3D, mesh, types3);
    auto colours3 = MeshColoring<3,0>::color(mesh);

    out3D << "CELL_DATA " << writer1.cellVert.getDataByPos<0>().size() << endl;
    out3D << "SCALARS cell_wrt_vertex_colour double 1\nLOOKUP_TABLE default" << endl;
    realIndex = 0;
    for (size_t i = 0; i < writer1.cellVert.getDataByPos<0>().size(); i++) {
        auto iterator = writer1.backwardCellIndexMapping.find(i);
        if (iterator == writer1.backwardCellIndexMapping.end()){
            out3D << colours3.getDataByPos<0>().at(realIndex) << ' ';
            realIndex++;
        } else {
            out3D << colours3.getDataByPos<0>().at(iterator->second) << ' ';
            realIndex = iterator->second;
        }
    }
    out3D.close();
}





void test3DMeshDeformedPrisms() {
    UnstructuredMesh<3, size_t, double, 6> mesh3;
    twoDeformedPrisms(mesh3);

    DBGVAR(mesh3.getSignature());
    //_ComputeCenters<1,3, 3,2,1>::compute<size_t, double, 6>(centers, mesh3);
    auto centers = ComputeCenters<DEFAULT>(mesh3);

    for(auto& face : mesh3.getFaces()) {
        face.setCenter(centers[face]);
        DBGVAR(face.getCenter());
    }
    DBGMSG("cellCenter");
    for(auto& cell : mesh3.getCells()) {
        cell.setCenter(centers[cell]);
        DBGVAR(cell.getCenter());
    }

    DBGMSG("measure computation");

    auto measures = ComputeMeasures<DEFAULT>(mesh3);
    for(double edgeM : measures.getDataByDim<1>()) {
        DBGVAR(edgeM);
    }
    for(double faceM : measures.getDataByDim<2>()) {
        DBGVAR(faceM);
    }

    for(double cellM : measures.getDataByDim<3>()) {
        DBGVAR(cellM);
    }


    MeshDataContainer<MeshNativeType<3>::ElementType,3> types(mesh3, MeshNativeType<3>::ElementType::WEDGE);

    VTKMeshWriter<3, size_t, double> writer;
    ofstream out3D("3D_test_mesh_two_deformed_prisms.vtk");
    writer.writeHeader(out3D, "test data");
    writer.writeToStream(out3D, mesh3, types);

}






void testMeshDataContainer() {
    UnstructuredMesh<3, size_t, double, 6> mesh3;
    twoDeformedPrisms(mesh3);


    MeshDataContainer<std::tuple<int, double, char, double>, 3,2,0> container(mesh3);


    MeshDataContainer<double, 3,3,1> cont(mesh3, 52.2);
    MeshDataContainer<double, 1,1> contAlloc;
    contAlloc.allocateData(cont, 42.2);

    DBGVAR(cont.getDataByDim<1>(), contAlloc.getDataByDim<1>());


    for(auto& c : container.getDataByDim<0>()) {
        c=42;
    }

    for(auto& v : mesh3.getVertices()){
        DBGVAR(container.at(v));
    }

    MeshDataContainer<std::tuple<int, double, char, int>, 3, 2, 0, 1> containerIni(mesh3,3, 42.15, 'a', 15);


    for (auto& val : containerIni.getDataByPos<0>()){
        DBGVAR(val);
    }

    for (auto& val : containerIni.getDataByPos<1>()){
        DBGVAR(val);
    }

    for (auto& val : containerIni.getDataByPos<2>()){
        DBGVAR(val);
    }

    for (auto& val : containerIni.getDataByPos<3>()){
        DBGVAR(val);
    }

    DBGMSG("assign test");
    MeshDataContainer<std::tuple<int, double, char, int>, 3, 2, 0, 1> containerAssign;

    containerAssign = containerIni;

    contAlloc.getDataByDim<1>().clear();
    contAlloc.allocateData(containerAssign, 87);


    DBGVAR(contAlloc.getDataByDim<1>(), contAlloc.getDataByPos<1>(),containerAssign.getDataByPos<0>());
}




void test3DMeshLoad() {
    UnstructuredMesh<3, size_t, double, 6> mesh;
    VTKMeshReader<3> reader;

    ifstream file("test_3Dmesh.vtk");
    reader.loadFromStream(file, mesh);

DBGVAR(mesh.getVertices().size(),mesh.getEdges().size(), mesh.getFaces().size(), mesh.getCells().size());
    DBGVAR(mesh.getVertices());


    DBGMSG("connection test");
    auto con2 = MeshConnections<1,0>::connections(mesh);
    for (auto& edge : mesh.getEdges()){
            DBGVAR(edge.getIndex(), con2[edge]);
    }

    for (auto& face : mesh.getFaces()) {
        DBGVAR(face.getIndex());
        for (auto& sube : face.getSubelements()){
            DBGVAR(sube.index);
        }
    }

    DBGMSG("connection test");
    auto con1 = MeshConnections<2,1>::connections(mesh);
    for (auto& face : mesh.getFaces()){
            DBGVAR(face.getIndex(), con1[face]);
    }

    DBGMSG("connection test");
    auto con = MeshConnections<3,2>::connections(mesh);
    for (auto& cell : mesh.getCells()){
            DBGVAR(cell.getIndex(), con[cell]);
    }

    mesh.initializeCenters();
    DBGVAR(mesh.computeElementMeasures().getDataByDim<3>(),ComputeCenters<DEFAULT>(mesh).getDataByDim<2>(),mesh.computeFaceNormals().getDataByPos<0>());



    VTKMeshWriter<3, size_t, double> writer;
    ofstream out3D("3D_test_mesh_output.vtk");
    writer.writeHeader(out3D, "test data");
    writer.writeToStream(out3D, mesh, reader.getCellTypes());

}


void testFPMARW(){
    UnstructuredMesh<3, size_t, double, 8> mesh;
    FPMAMeshReader<3> reader;
    ifstream file("Poly_simple.fpma");
    reader.loadFromStream(file, mesh);

    DBGVAR(mesh.getCells().size(), mesh.getVertices().size());

    auto faceVert = MeshConnections<2,0,ORDER_ORIGINAL>::connections(mesh);
    DBGVAR(faceVert.getDataByPos<0>());

    mesh.initializeCenters();
    VTKMeshWriter<3, size_t, double> writer;
    ofstream ofile("Poly_simple.vtk");
    writer.writeHeader(ofile, "fpma_output_test");
    writer.writeToStream(ofile, mesh, MeshDataContainer<MeshNativeType<3>::ElementType, 3>(mesh, MeshNativeType<3>::POLYHEDRON));

    ofile << "CELL_DATA " << writer.cellVert.getDataByPos<0>().size() << endl;
    ofile << "SCALARS cell_wrt_vertex_colour double 1\nLOOKUP_TABLE default" << endl;
    size_t realIndex = 0;
    for (size_t i = 0; i < writer.cellVert.getDataByPos<0>().size(); i++) {
        auto iterator = writer.backwardCellIndexMapping.find(i);
        if (iterator == writer.backwardCellIndexMapping.end()){
            ofile << realIndex << ' ';
            realIndex++;
        } else {
            ofile << iterator->second << ' ';
            realIndex = iterator->second;
        }
    }
    ofile.close();

    FPMAMeshWriter<3,size_t, double> writerFPMA;
    ofstream ofile1("Poly_simple_test.fpma");
    writerFPMA.writeToStream(ofile1, mesh);


}


void testFPMA_poly(){
    UnstructuredMesh<3, size_t, double, 0> mesh;
    FPMAMeshReader<3> reader;
    ifstream file("Spark_mesh.fpma");
    reader.loadFromStream(file, mesh);

    DBGVAR_CSV(mesh.getCells().size(), mesh.getVertices().size());

    mesh.initializeCenters();
    VTKMeshWriter<3, size_t, double> writer;
    ofstream ofile("Spark_mesh.vtk");
    writer.writeHeader(ofile, "fpma_output_test");
    writer.writeToStream(ofile, mesh, MeshDataContainer<MeshNativeType<3>::ElementType, 3>(mesh, MeshNativeType<3>::POLYHEDRON));

    ofile << "CELL_DATA " << writer.cellVert.getDataByPos<0>().size() << endl;
    ofile << "SCALARS cell_wrt_vertex_colour double 1\nLOOKUP_TABLE default" << endl;
    size_t realIndex = 0;
    for (size_t i = 0; i < writer.cellVert.getDataByPos<0>().size(); i++) {
        auto iterator = writer.backwardCellIndexMapping.find(i);
        if (iterator == writer.backwardCellIndexMapping.end()){
            ofile << realIndex << ' ';
            realIndex++;
        } else {
            ofile << iterator->second << ' ';
            realIndex = iterator->second;
        }
    }
    ofile.close();

}


int main()
{

    testMesh2D();
    //testMesh2DLoadAndWrite();
    testMesh3D();
    test3DMeshDeformedPrisms();
    testMeshRefine();
    testMeshDataContainer();
    //UnstructuredMesh<5, size_t, double, 6,5,4> m;
    //m.ComputeElementMeasures();
    //test3DMeshLoad();

    //testFPMA_poly();
    DBGVAR(vector<double>());

}

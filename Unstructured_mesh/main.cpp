#include <iostream>
//#define UNDEBUG
#include "../debug/Debug.h"
#include "UnstructuredMesh.h"
#include "MeshFunctions.h"
#include "VTKMeshReader.h"
#include "VTKMeshWriter.h"
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
        mesh3.getFaces().at(0).getSubelements().addSubelement(0,true);
        mesh3.getFaces().at(0).getSubelements().addSubelement(1,true);
        mesh3.getFaces().at(0).getSubelements().addSubelement(2,true);
        mesh3.getFaces().at(0).getSubelements().addSubelement(3,true);
        mesh3.getFaces().push_back(1);
        mesh3.getFaces().at(1).getSubelements().addSubelement(0,true);
        mesh3.getFaces().at(1).getSubelements().addSubelement(5,true);
        mesh3.getFaces().at(1).getSubelements().addSubelement(8,true);
        mesh3.getFaces().at(1).getSubelements().addSubelement(4,true);
        mesh3.getFaces().push_back(2);
        mesh3.getFaces().at(2).getSubelements().addSubelement(1,true);
        mesh3.getFaces().at(2).getSubelements().addSubelement(4,true);
        mesh3.getFaces().at(2).getSubelements().addSubelement(11,true);
        mesh3.getFaces().at(2).getSubelements().addSubelement(7,true);
        mesh3.getFaces().push_back(3);
        mesh3.getFaces().at(3).getSubelements().addSubelement(3,true);
        mesh3.getFaces().at(3).getSubelements().addSubelement(6,true);
        mesh3.getFaces().at(3).getSubelements().addSubelement(10,true);
        mesh3.getFaces().at(3).getSubelements().addSubelement(7,true);
        mesh3.getFaces().push_back(4);
        mesh3.getFaces().at(4).getSubelements().addSubelement(2,true);
        mesh3.getFaces().at(4).getSubelements().addSubelement(6,true);
        mesh3.getFaces().at(4).getSubelements().addSubelement(9,true);
        mesh3.getFaces().at(4).getSubelements().addSubelement(5,true);
        mesh3.getFaces().push_back(5);
        mesh3.getFaces().at(5).getSubelements().addSubelement(8,true);
        mesh3.getFaces().at(5).getSubelements().addSubelement(9,true);
        mesh3.getFaces().at(5).getSubelements().addSubelement(10,true);
        mesh3.getFaces().at(5).getSubelements().addSubelement(11,true);
    DBGCHECK;

        mesh3.getFaces().at(0).setNextBElem(1,0);
        mesh3.getFaces().at(1).setNextBElem(2,0);
        mesh3.getFaces().at(2).setNextBElem(3,0);
        mesh3.getFaces().at(3).setNextBElem(4,0);
        mesh3.getFaces().at(4).setNextBElem(5,0);
        mesh3.getFaces().at(5).setNextBElem(0,0);

        mesh3.getCells().push_back(0);
        mesh3.getCells().at(0).setBoundaryElementIndex(3);
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
        mesh3.getFaces().at(index).getSubelements().addSubelement(0,true);
        mesh3.getFaces().at(index).getSubelements().addSubelement(1,true);
        mesh3.getFaces().at(index).getSubelements().addSubelement(13,true);
        mesh3.getFaces().push_back(++index);
        mesh3.getFaces().at(index).getSubelements().addSubelement(13,true);
        mesh3.getFaces().at(index).getSubelements().addSubelement(2,true);
        mesh3.getFaces().at(index).getSubelements().addSubelement(3,true);
        mesh3.getFaces().push_back(++index);
        mesh3.getFaces().at(index).getSubelements().addSubelement(0,true);
        mesh3.getFaces().at(index).getSubelements().addSubelement(5,true);
        mesh3.getFaces().at(index).getSubelements().addSubelement(8,true);
        mesh3.getFaces().at(index).getSubelements().addSubelement(4,true);
        mesh3.getFaces().push_back(++index);
        mesh3.getFaces().at(index).getSubelements().addSubelement(1,true);
        mesh3.getFaces().at(index).getSubelements().addSubelement(4,true);
        mesh3.getFaces().at(index).getSubelements().addSubelement(11,true);
        mesh3.getFaces().at(index).getSubelements().addSubelement(7,true);
        mesh3.getFaces().push_back(++index);
        mesh3.getFaces().at(index).getSubelements().addSubelement(3,true);
        mesh3.getFaces().at(index).getSubelements().addSubelement(6,true);
        mesh3.getFaces().at(index).getSubelements().addSubelement(10,true);
        mesh3.getFaces().at(index).getSubelements().addSubelement(7,true);
        mesh3.getFaces().push_back(++index);
        mesh3.getFaces().at(index).getSubelements().addSubelement(2,true);
        mesh3.getFaces().at(index).getSubelements().addSubelement(6,true);
        mesh3.getFaces().at(index).getSubelements().addSubelement(9,true);
        mesh3.getFaces().at(index).getSubelements().addSubelement(5,true);
        mesh3.getFaces().push_back(++index);
        mesh3.getFaces().at(index).getSubelements().addSubelement(8,true);
        mesh3.getFaces().at(index).getSubelements().addSubelement(12,true);
        mesh3.getFaces().at(index).getSubelements().addSubelement(11,true);
        mesh3.getFaces().push_back(++index);
        mesh3.getFaces().at(index).getSubelements().addSubelement(9,true);
        mesh3.getFaces().at(index).getSubelements().addSubelement(10,true);
        mesh3.getFaces().at(index).getSubelements().addSubelement(12,true);
        mesh3.getFaces().push_back(++index);
        mesh3.getFaces().at(index).getSubelements().addSubelement(12,true);
        mesh3.getFaces().at(index).getSubelements().addSubelement(7,true);
        mesh3.getFaces().at(index).getSubelements().addSubelement(13,true);
        mesh3.getFaces().at(index).getSubelements().addSubelement(5,true);
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
        mesh3.getFaces().at(index).getSubelements().addSubelement(0,true);
        mesh3.getFaces().at(index).getSubelements().addSubelement(1,true);
        mesh3.getFaces().at(index).getSubelements().addSubelement(13,true);
        mesh3.getFaces().push_back(++index);
        mesh3.getFaces().at(index).getSubelements().addSubelement(13,true);
        mesh3.getFaces().at(index).getSubelements().addSubelement(2,true);
        mesh3.getFaces().at(index).getSubelements().addSubelement(3,true);
        mesh3.getFaces().push_back(++index);
        mesh3.getFaces().at(index).getSubelements().addSubelement(0,true);
        mesh3.getFaces().at(index).getSubelements().addSubelement(5,true);
        mesh3.getFaces().at(index).getSubelements().addSubelement(8,true);
        mesh3.getFaces().at(index).getSubelements().addSubelement(4,true);
        mesh3.getFaces().push_back(++index);
        mesh3.getFaces().at(index).getSubelements().addSubelement(1,true);
        mesh3.getFaces().at(index).getSubelements().addSubelement(4,true);
        mesh3.getFaces().at(index).getSubelements().addSubelement(11,true);
        mesh3.getFaces().at(index).getSubelements().addSubelement(7,true);
        mesh3.getFaces().push_back(++index);
        mesh3.getFaces().at(index).getSubelements().addSubelement(3,true);
        mesh3.getFaces().at(index).getSubelements().addSubelement(6,true);
        mesh3.getFaces().at(index).getSubelements().addSubelement(10,true);
        mesh3.getFaces().at(index).getSubelements().addSubelement(7,true);
        mesh3.getFaces().push_back(++index);
        mesh3.getFaces().at(index).getSubelements().addSubelement(2,true);
        mesh3.getFaces().at(index).getSubelements().addSubelement(6,true);
        mesh3.getFaces().at(index).getSubelements().addSubelement(9,true);
        mesh3.getFaces().at(index).getSubelements().addSubelement(5,true);
        mesh3.getFaces().push_back(++index);
        mesh3.getFaces().at(index).getSubelements().addSubelement(8,true);
        mesh3.getFaces().at(index).getSubelements().addSubelement(12,true);
        mesh3.getFaces().at(index).getSubelements().addSubelement(11,true);
        mesh3.getFaces().push_back(++index);
        mesh3.getFaces().at(index).getSubelements().addSubelement(9,true);
        mesh3.getFaces().at(index).getSubelements().addSubelement(10,true);
        mesh3.getFaces().at(index).getSubelements().addSubelement(12,true);
        mesh3.getFaces().push_back(++index);
        mesh3.getFaces().at(index).getSubelements().addSubelement(12,true);
        mesh3.getFaces().at(index).getSubelements().addSubelement(7,true);
        mesh3.getFaces().at(index).getSubelements().addSubelement(13,true);
        mesh3.getFaces().at(index).getSubelements().addSubelement(5,true);
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



    auto centers = ComputeCenters(mesh);

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

    auto normals = ComputeFaceNormals(mesh);
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
    VTKMeshReader<2, size_t, double> reader;
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




    DBGMSG("mesh conatiner test");
    MeshDataContainer<double, 3,2,1,0> cont(mesh3);

    MakeMeshDataContainer_t<double, make_custom_integer_sequence_t<unsigned int, 0,3>> cont1(mesh3);


    //cont.getDataByDim<3>().resize(20);
    DBGVAR(cont.getDataByPos<1>().size(), cont.getDataByPos<1>().getMappedDimension());

    DBGVAR(cont.getDataByPos<0>().size());

    DBGVAR(cont.getDataByDim<3>().size());


    DBGMSG("faceCenters");
    //MeshDataContainer<Vertex<3,double>, 3,2,1> centers(mesh3);


    //_ComputeCenters<1,3, 3,2,1>::compute<size_t, double, 6>(centers, mesh3);
    auto centers = ComputeCenters(mesh3);

    for(auto& face : mesh3.getFaces()) {
        face.setCenter(centers.template getDataByDim<2>().at(face.getIndex()));
        DBGVAR(face.getCenter());
    }
    DBGMSG("cellCenter");
    for(auto& cell : mesh3.getCells()) {
        cell.setCenter(centers.template getDataByDim<3>().at(cell.getIndex()));
        DBGVAR(cell.getCenter());
    }

    DBGMSG("measure computation");

    auto measures = ComputeMeasures(mesh3);
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

    DBGMSG("mesh apply test");
    MeshApply<3, 2, 3>::apply(mesh3, [](size_t ori, size_t i){
        DBGVAR(ori,i);
    });
    DBGMSG("mesh apply test");
    MeshApply<2, 3, 3>::apply(mesh3,[](size_t ori, size_t i){
        DBGVAR(ori,i);
    });

    DBGMSG("3D edge orientation");
    MeshApply<2, 1, 3>::apply(mesh3,[&mesh3](size_t faceIndex, size_t edgeIndex){
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


    DBGMSG("connection test oposite");
    auto con1 = MeshConnections<0,3>::connections(mesh3);
    for (auto& vert : mesh3.getVertices()){
        DBGVAR(vert.getIndex(), con1[vert]);
    }

    DBGMSG("face to vertex colouring");
    auto colours = ColourMesh<2,0>::colour(mesh3);
    for (auto& face : mesh3.getFaces()){
        DBGVAR(face.getIndex(), colours.at(face));
    }

    DBGMSG("vertex to face colouring");
    auto colours1 = ColourMesh<0,2>::colour(mesh3);
    for (auto& vert : mesh3.getVertices()){
        DBGVAR(vert.getIndex(), colours1.at(vert));
    }

    MeshDataContainer<MeshNativeType<3>::ElementType,3> types(mesh3, MeshNativeType<3>::ElementType::WEDGE);

    VTKMeshWriter<3, size_t, double, 6> writer;
    ofstream out3D("3D_test_mesh_two_prisms.vtk");
    writer.writeHeader(out3D, "test data");
    writer.writeToStream(out3D, mesh3, types);
}




void test3DMeshDeformedPrisms() {
    UnstructuredMesh<3, size_t, double, 6> mesh3;
    twoDeformedPrisms(mesh3);

    //_ComputeCenters<1,3, 3,2,1>::compute<size_t, double, 6>(centers, mesh3);
    auto centers = ComputeCenters(mesh3);

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

    auto measures = ComputeMeasures(mesh3);
    for(double edgeM : measures.getDataByDim<1>()) {
        DBGVAR(edgeM);
    }
    for(double faceM : measures.getDataByDim<2>()) {
        DBGVAR(faceM);
    }

    for(double cellM : measures.getDataByDim<3>()) {
        DBGVAR(cellM);
    }
}








void testMeshDataContainer() {
    UnstructuredMesh<3, size_t, double, 6> mesh3;
    twoDeformedPrisms(mesh3);


    MeshDataContainer<std::tuple<int, double, char, double>, 3,2,0> container(mesh3);



    for(auto& c : container.getDataByDim<0>()) {
        c=42;
    }

    for(auto& v : mesh3.getVertices()){
        DBGVAR(container.at(v));
    }

    MeshDataContainer<std::tuple<int, double, char, int>, 3, 2, 0, 2> containerIni(mesh3,3, 42.15, 'a', 15);

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
}




void test3DMeshLoad() {
    UnstructuredMesh<3, size_t, double, 6> mesh;
    VTKMeshReader<3, size_t, double, 6> reader;

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
    DBGVAR(mesh.computeElementMeasures().getDataByDim<3>(),ComputeCenters(mesh).getDataByDim<2>(),mesh.computeFaceNormals().getDataByPos<0>());

}







int main()
{
    testMesh2D();
    //testMesh2DLoadAndWrite();
    testMesh3D();
    //test3DMeshDeformedPrisms();
    //testMeshDataContainer();
    //UnstructuredMesh<5, size_t, double, 6,5,4> m;
    //m.ComputeElementMeasures();
    //test3DMeshLoad();


}

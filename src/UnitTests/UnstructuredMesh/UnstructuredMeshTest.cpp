// Test of Traits class
#ifdef HAVE_GTEST
#include <gtest/gtest.h>
#else
#define TEST(_1,_2) void _1()
#define EXPECT_TRUE(_1) DBGCHECK;(void)(_1)
#define EXPECT_FALSE(_1) DBGCHECK;(void)(_1)
#define EXPECT_EQ(_1,_2) DBGCHECK;(void)(_1 == _2)
#define EXPECT_ANY_THROW(_1) DBGCHECK;try{(_1);}catch(...){}
#endif
#include "GTMesh/Debug/Debug.h"
#include <list>
#include <map>
#include <array>
#include <string>
#include <fstream>
#include <ios>
#include "GTMesh/UnstructuredMesh/UnstructuredMesh.h"
#include "GTMesh/UnstructuredMesh/MeshDataContainer/MeshDataIO/VTKMeshDataReader.h"
#include "GTMesh/UnstructuredMesh/MeshDataContainer/MeshDataIO/VTKMeshDataWriter.h"
#include "MeshSetup.h"


bool floatArrayCompare(const double& _1, const double& _2, double treshold = 1e-8){
    return fabs(_1 - _2) < treshold;
}

template<typename ArrayT, typename ArrayT2, typename ..., std::enable_if_t<IsIndexable<ArrayT>::value && IsIndexable<ArrayT2>::value,bool> = true>
bool floatArrayCompare(const ArrayT& _1, const ArrayT2& _2, double treshold = 1e-8){
    for (decltype (_1.size()) i = 0; i < _1.size(); i++) {
        if (!floatArrayCompare(_1[i], _2[i], treshold)){
            return false;
        }
    }
    return true;
}


template<unsigned int ColoredDim, unsigned int ConnectingDim, typename Mesh>
bool testProperColoring(const Mesh& mesh, const MeshDataContainer<unsigned int, ColoredDim>& coloring) {
    auto neighbors = mesh.template neighborhood<ColoredDim, ConnectingDim>();
    for (auto element : mesh.template getElements<ColoredDim>()){
        for (auto neighborIndex : neighbors[element]){
            if (coloring[element] == coloring.template getDataByPos<0>()[neighborIndex]){
                return false;
            }
        }
    }
    return true;
}

TEST( UnstructuredMesh2D_Functions_Test, basicTest )
{
    using MeshType = UnstructuredMesh<2, size_t, double>;
    MeshType mesh;

    square(mesh); // square mesh made of 2 triangles

    EXPECT_TRUE(mesh.updateSignature() != 0);
    mesh.initializeCenters();


    auto& face = mesh.getFaces()[0];
    EXPECT_TRUE(isInvalidIndex(face.getCellRightIndex()));

    mesh.setupBoundaryCells();
    mesh.setupBoundaryCellsCenters();
    face = mesh.getFaces()[0];
    EXPECT_TRUE(floatArrayCompare( face.getCenter(), (Vector<2, double>{0, 0.5})));
    EXPECT_TRUE(isBoundaryIndex(face.getCellRightIndex()));
    EXPECT_EQ(extractBoundaryIndex(face.getCellRightIndex()), 0);
    EXPECT_TRUE(face.getCellLeftIndex() == 0);

    EXPECT_EQ(face.getOtherCellIndex(makeBoundaryIndex<size_t>(0)), 0);
    EXPECT_EQ(face.getOtherCellIndex(0), makeBoundaryIndex<size_t>(0));

    EXPECT_EQ(face.getNextBElem(face.getCellLeftIndex()), 1);
    EXPECT_EQ(face.getNextBElem(face.getCellRightIndex()), INVALID_INDEX(size_t));

    EXPECT_ANY_THROW(face.getOtherCellIndex(1));
    EXPECT_ANY_THROW(face.getNextBElem(1));

    // compute centers
    auto centers = computeCenters<METHOD_DEFAULT>(mesh);
    std::vector<Vertex<2, double>> expectCenter = { {0.333333, 0.333333}, {0.666667, 0.666667} };
    EXPECT_EQ((centers.getDataByDim<2>()), expectCenter);

    // measures test
    auto measures = mesh.computeElementMeasures();
    std::vector<double> expectEdgeM = { 1.0, 1.0, 1.41421, 1.0, 1.0 };
    EXPECT_EQ(measures.getDataByDim<1>(), expectEdgeM);
    std::vector<double> expectCellM = { 0.5, 0.5 };
    EXPECT_EQ(measures.getDataByDim<2>(), expectCellM);

    // face normals test
    auto normals = mesh.computeFaceNormals();
    std::vector<Vector<2, double>> expectNormals = { {-1, 0}, {0, -1}, { 0.707107, 0.707107}, {0, 1}, {1, 0} };
    EXPECT_EQ(normals.getDataByPos<0>(), expectNormals);

    // center distance
    auto dist = computeCellsDistance(mesh);
    std::vector<double> expectCellDist = { 0.372678, 0.372678, 0.471405, 0.372678, 0.372678 };
    EXPECT_EQ(dist.getDataByPos<0>(), expectCellDist);

    // Mesh connections
    std::vector<std::vector<size_t>> expCon20 = { { 0, 1, 2 }, { 1, 2, 3 } };
    std::vector<std::vector<size_t>> expCon21 = { { 0, 1, 2 }, { 2, 3, 4 } };
    std::vector<std::vector<size_t>> expCon12 = { { 0 }, { 0 }, { 0, 1 }, { 1 }, { 1 } };
    std::vector<std::vector<size_t>> expCon02 = { { 0 }, { 0, 1 }, { 0, 1 }, { 1 } };

    EXPECT_EQ((mesh.connections<2,0>().getDataByPos<0>()), expCon20);
    EXPECT_EQ((mesh.connections<2,1>().getDataByPos<0>()), expCon21);
    EXPECT_EQ((mesh.connections<1,2>().getDataByPos<0>()), expCon12);
    EXPECT_EQ((mesh.connections<0,2>().getDataByPos<0>()), expCon02);


    std::vector<std::vector<size_t>> expNeighbors20 = { { 1 }, { 0 } };
    std::vector<std::vector<size_t>> expNeighbors210 = { { 0, 1, 2 }, { 1, 2, 3 } };
    std::vector<std::vector<size_t>> expNeighbors120 = { { 0, 1, 2 }, { 0, 1, 2 }, { 0, 1, 2, 3 }, { 1, 2, 3 }, { 1, 2, 3 } };
    std::vector<std::vector<size_t>> expNeighbors02 = { { 1, 2 }, { 0, 2, 3 }, { 0, 1, 3 }, { 1, 2 } };

    EXPECT_EQ((mesh.neighborhood<2,0>().getDataByPos<0>()),   expNeighbors20);
    EXPECT_EQ((mesh.neighborhood<2,1,0>().getDataByPos<0>()), expNeighbors210);
    EXPECT_EQ((mesh.neighborhood<1,2,0>().getDataByPos<0>()), expNeighbors120);
    EXPECT_EQ((mesh.neighborhood<0,2>().getDataByPos<0>()),   expNeighbors02);


    auto greedyColor10 = mesh.coloring< 1, 0 >();
    EXPECT_TRUE((testProperColoring<1, 0>(mesh, greedyColor10)));
    auto greedyColor01 = mesh.coloring< 0, 1 >();
    EXPECT_TRUE((testProperColoring<0, 1>(mesh, greedyColor01)));

    auto randomColor01 = mesh.coloring< 1, 0, METHOD_RANDOM >();
    EXPECT_TRUE((testProperColoring<1, 0>(mesh, randomColor01)));

}

template<unsigned int Dim>
struct CellData {
    unsigned int color;
    Vertex<Dim, double> center;
};

MAKE_ATTRIBUTE_TEMPLATE_TRAIT((CellData<Dim>), (unsigned int Dim), center, color);

TEST( UnstructuredMesh2DReadWrite, basicTest )
{
    UnstructuredMesh<2, size_t, double> mesh;
    square(mesh);
    mesh.initializeCenters();

    auto col = mesh.coloring<2,1>();

    MeshDataContainer<CellData<2>, 2> cellData(mesh);
    for (auto& cell : mesh.getCells()){
        cellData[cell].color = col[cell];
        cellData[cell].center = cell.getCenter();
    }

    MeshDataContainer<MeshNativeType<2>::ElementType,2> types(mesh, MeshNativeType<2>::TRIANGLE);
    VTKMeshWriter<2, size_t, double> writer;
    std::ofstream out2D;
    out2D.open("test_mesh_2D.vtk");
    EXPECT_TRUE(bool(out2D));
    writer.writeHeader(out2D, "test data");
    writer.writeToStream(out2D, mesh, types);
    VTKMeshDataWriter<2>::writeToStream(out2D, cellData, writer);
    out2D.close();

    mesh.clear();


    VTKMeshReader<2> reader;
    std::ifstream ifst("test_mesh_2D.vtk", std::ios::binary | std::ios::in);
    EXPECT_TRUE(bool(ifst));
    reader.loadFromStream(ifst, mesh);
    MeshDataContainer<CellData<2>, 2> cellDataLoad(mesh);
    VTKMeshDataReader<2, size_t>::readFromStream(ifst, cellDataLoad);
    ifst.close();

    mesh.initializeCenters();
    auto measures = mesh.computeElementMeasures();
    std::vector<double> expectEdgeM = { 1.0, 1.0, 1.41421, 1.0, 1.0 };
    EXPECT_EQ(measures.getDataByDim<1>(), expectEdgeM);
    std::vector<double> expectCellM = { 0.5, 0.5 };
    EXPECT_EQ(measures.getDataByDim<2>(), expectCellM);

    for (auto& cell : mesh.getCells()){
        EXPECT_EQ(cellDataLoad[cell].color, cellData[cell].color);
        EXPECT_TRUE(floatArrayCompare( cellDataLoad[cell].center, cellData[cell].center, 1e-4 ));
    }
}




TEST( UnstructuredMesh3D_Functions_Test, 3DMeshTest )
{
    using MeshType = UnstructuredMesh<3, size_t, double>;
    MeshType mesh;

    twoPrisms(mesh); // block made of 2 prisms

    EXPECT_TRUE(mesh.updateSignature() != 0);
    mesh.initializeCenters();


    auto& face = mesh.getFaces()[0];
    EXPECT_TRUE(isInvalidIndex(face.getCellRightIndex()));

    mesh.setupBoundaryCells();
    mesh.setupBoundaryCellsCenters();
    EXPECT_TRUE(isBoundaryIndex(face.getCellRightIndex()));
    EXPECT_EQ(extractBoundaryIndex(face.getCellRightIndex()), 0);
    EXPECT_TRUE(face.getCellLeftIndex() == 0);

    EXPECT_EQ(face.getOtherCellIndex(makeBoundaryIndex<size_t>(0)), 0);
    EXPECT_EQ(face.getOtherCellIndex(0), makeBoundaryIndex(0));

    EXPECT_EQ(face.getNextBElem(face.getCellLeftIndex()), 1);
    EXPECT_EQ(face.getNextBElem(face.getCellRightIndex()), INVALID_INDEX(size_t));

    EXPECT_ANY_THROW(face.getOtherCellIndex(1));
    EXPECT_ANY_THROW(face.getNextBElem(1));


    // compute centers
    auto centers = computeCenters<METHOD_DEFAULT>(mesh);
    std::vector<Vertex<3, double>> expectCenter = { {0.333333, 0.333333, 0.5}, {0.666667, 0.666667, 0.5} };
    EXPECT_EQ((centers.getDataByDim<3>()), expectCenter);

    std::vector<Vertex<3, double>> expectCenterFace = {{ 0.333333, 0.333333, 0 }, { 0.666667, 0.666667, 0 }, { 0.5, 0, 0.5 }, { 0, 0.5, 0.5 }, { 0.5, 1, 0.5 }, { 1, 0.5, 0.5 }, { 0.333333, 0.333333, 1 }, { 0.666667, 0.666667, 1 }, { 0.5, 0.5, 0.5 }};
    EXPECT_TRUE(floatArrayCompare((centers.getDataByDim<3>()), expectCenterFace));


    // measures test
    auto measures = mesh.computeElementMeasures();
    auto measuresTess = mesh.computeElementMeasures<METHOD_TESSELLATED>();

    std::vector<double> expectEdgeM = { 1, 1.41421, 1, 1, 1, 1, 1, 1.41421, 1, 1, 1, 1, 1, 1 };
    EXPECT_EQ(measures.getDataByDim<1>(), expectEdgeM);
    EXPECT_EQ(measuresTess.getDataByDim<1>(), expectEdgeM);
    std::vector<double> expectFaceM = {  0.5, 1, 1.41421, 1, 0.5, 0.5, 1, 1, 0.5 };
    EXPECT_EQ(measures.getDataByDim<2>(), expectFaceM);
    EXPECT_EQ(measuresTess.getDataByDim<2>(), expectFaceM);
    std::vector<double> expectCellM = {  0.5, 0.5 };
    EXPECT_TRUE(floatArrayCompare(measures.getDataByDim<3>(), expectCellM));
    EXPECT_TRUE(floatArrayCompare(measuresTess.getDataByDim<3>(), expectCellM));
    // face normals test
    auto normals = mesh.computeFaceNormals();
    std::vector<Vector<3, double>> expectNormals = {{ 0, 0, -1 }, { 0, 0, -1 }, { 0, -1, 0 }, { -1, 0, 0 }, { 0, 1, 0 }, { 1, 0, 0 }, { 0, 0, 1 }, { 0, 0, 1 }, { 0.707107, 0.707107, 0 }};
    EXPECT_TRUE(floatArrayCompare( normals.getDataByPos<0>(), expectNormals ));
    EXPECT_TRUE(floatArrayCompare( mesh.computeFaceNormals<METHOD_TESSELLATED>().getDataByPos<0>(), expectNormals ));

    // center distance
    auto dist = computeCellsDistance(mesh);
    std::vector<double> expectCellDist = { 0.687184, 1.06719, 0.687184, 0.687184, 1.06719, 1.06719, 0.687184, 1.06719, 0.471405 };
    EXPECT_EQ(dist.getDataByPos<0>(), expectCellDist);

    // Mesh connections
    std::vector<std::vector<size_t>> expCon20 = { { 0, 1, 2 }, { 1, 2, 3 }, { 0, 1, 4, 5 }, { 0, 2, 4, 6 }, { 2, 3, 6, 7 }, { 1, 3, 5, 7 }, { 4, 5, 6 }, { 5, 6, 7 }, { 1, 2, 5, 6 } };
    std::vector<std::vector<size_t>> expCon30 = { { 0, 1, 2, 4, 5, 6 }, { 1, 2, 3, 5, 6, 7 } };
    std::vector<std::vector<size_t>> expCon03 = { { 0 }, { 0, 1 }, { 0, 1 }, { 1 }, { 0 }, { 0, 1 }, { 0, 1 }, { 1 } };
    std::vector<std::vector<size_t>> expCon13 = { { 0 }, { 0 }, { 1 }, { 1 }, { 0 }, { 0, 1 }, { 1 }, { 0, 1 }, { 0 }, { 1 }, { 1 }, { 0 }, { 0, 1 }, { 0, 1 } };

    EXPECT_EQ((mesh.connections<2,0>().getDataByPos<0>()), expCon20);
    EXPECT_EQ((mesh.connections<3,0>().getDataByPos<0>()), expCon30);
    EXPECT_EQ((mesh.connections<0,3>().getDataByPos<0>()), expCon03);
    EXPECT_EQ((mesh.connections<1,3>().getDataByPos<0>()), expCon13);

    // Test neighbors
    std::vector<std::vector<size_t>> expNeighbors20 = { { 1, 2, 3, 4, 5, 8 }, { 0, 2, 3, 4, 5, 8 }, { 0, 1, 3, 5, 6, 7, 8 }, { 0, 1, 2, 4, 6, 7, 8 }, { 0, 1, 3, 5, 6, 7, 8 }, { 0, 1, 2, 4, 6, 7, 8 }, { 2, 3, 4, 5, 7, 8 }, { 2, 3, 4, 5, 6, 8 }, { 0, 1, 2, 3, 4, 5, 6, 7 } };
    std::vector<std::vector<size_t>> expNeighbors230 = { { 0, 1, 2, 4, 5, 6 }, { 1, 2, 3, 5, 6, 7 }, { 0, 1, 2, 4, 5, 6 }, { 0, 1, 2, 4, 5, 6 }, { 1, 2, 3, 5, 6, 7 }, { 1, 2, 3, 5, 6, 7 }, { 0, 1, 2, 4, 5, 6 }, { 1, 2, 3, 5, 6, 7 }, { 0, 1, 2, 3, 4, 5, 6, 7 } };
    std::vector<std::vector<size_t>> expNeighbors020 = { { 0, 1, 2, 4, 5, 6 }, { 1, 2, 3, 5, 6, 7 }, { 0, 1, 2, 4, 5, 6 }, { 0, 1, 2, 4, 5, 6 }, { 1, 2, 3, 5, 6, 7 }, { 1, 2, 3, 5, 6, 7 }, { 0, 1, 2, 4, 5, 6 }, { 1, 2, 3, 5, 6, 7 }, { 0, 1, 2, 3, 4, 5, 6, 7 } };
    std::vector<std::vector<size_t>> expNeighbors320 = { { 0, 1, 2, 4, 5, 6 }, { 1, 2, 3, 5, 6, 7 } };
    std::vector<std::vector<size_t>> expNeighbors02 = { { 1, 2, 4, 5, 6 }, { 0, 2, 3, 4, 5, 6, 7 }, { 0, 1, 3, 4, 5, 6, 7 }, { 1, 2, 5, 6, 7 }, { 0, 1, 2, 5, 6 }, { 0, 1, 2, 3, 4, 6, 7 }, { 0, 1, 2, 3, 4, 5, 7 }, { 1, 2, 3, 5, 6 } };

    EXPECT_EQ((mesh.neighborhood<2,0>().getDataByPos<0>()),   expNeighbors20);
    EXPECT_EQ((mesh.neighborhood<2,3,0>().getDataByPos<0>()), expNeighbors230);
    EXPECT_EQ((mesh.neighborhood<0,2,0>().getDataByPos<0>()), expNeighbors020);
    EXPECT_EQ((mesh.neighborhood<3,2,0>().getDataByPos<0>()), expNeighbors320);
    EXPECT_EQ((mesh.neighborhood<0,2>().getDataByPos<0>()),    expNeighbors02);

    // Test proper coloring
    EXPECT_TRUE((testProperColoring<1, 0>(mesh, mesh.coloring< 1, 0 >())));
    EXPECT_TRUE((testProperColoring<0, 1>(mesh, mesh.coloring< 0, 1 >())));
    EXPECT_TRUE((testProperColoring<0, 3>(mesh, mesh.coloring< 0, 3 >())));
    EXPECT_TRUE((testProperColoring<2, 0>(mesh, mesh.coloring< 2, 0 >())));

    EXPECT_TRUE((testProperColoring<1, 0>(mesh, mesh.coloring< 1, 0, METHOD_RANDOM >())));
    EXPECT_TRUE((testProperColoring<0, 1>(mesh, mesh.coloring< 0, 1, METHOD_RANDOM >())));
    EXPECT_TRUE((testProperColoring<0, 3>(mesh, mesh.coloring< 0, 3, METHOD_RANDOM >())));
    EXPECT_TRUE((testProperColoring<2, 0>(mesh, mesh.coloring< 2, 0, METHOD_RANDOM >())));

}

TEST( MeshRefineTest, 3DMeshTest ) {
    UnstructuredMesh<3, size_t, double, 6> mesh;
    twoPrisms(mesh);
    mesh.initializeCenters();

    auto col = mesh.coloring<3,1>();

    MeshDataContainer<CellData<3>, 3> cellData(mesh);
    for (auto& cell : mesh.getCells()){
        cellData[cell].color = col[cell];
        cellData[cell].center = cell.getCenter();
    }

    MeshDataContainer<MeshNativeType<3>::ElementType,3> types(mesh, MeshNativeType<3>::WEDGE);
    VTKMeshWriter<3, size_t, double> writer;
    std::ofstream out3D;
    out3D.open("mesh_refine_0.vtk");
    EXPECT_TRUE(bool(out3D));
    writer.writeHeader(out3D, "test data");
    writer.writeToStream(out3D, mesh, types);
    VTKMeshDataWriter<3>::writeToStream(out3D, cellData, writer);
    out3D.close();

    mesh.clear();

    MeshDataContainer<CellData<3>, 3> cellDataLoad(mesh);

    VTKMeshReader<3> reader;
    std::ifstream ifst("mesh_refine_0.vtk");
    EXPECT_TRUE(bool(ifst));
    reader.loadFromStream(ifst, mesh);
    VTKMeshDataReader<3, size_t>::readFromStream(ifst, cellDataLoad);
    ifst.close();

    mesh.initializeCenters();
    // compute centers
    auto centers = computeCenters<METHOD_DEFAULT>(mesh);
    std::vector<Vertex<3, double>> expectCenter = { {0.333333, 0.333333, 0.5}, {0.666667, 0.666667, 0.5} };
    EXPECT_TRUE(floatArrayCompare((centers.getDataByDim<3>()), expectCenter));

    // measures test
    auto measures = mesh.computeElementMeasures();

    std::vector<double> expectEdgeM = { 1, 1.41421, 1, 1, 1, 1, 1, 1.41421, 1, 1, 1, 1, 1, 1 };
    EXPECT_TRUE(floatArrayCompare(measures.getDataByDim<1>(), expectEdgeM));
    std::vector<double> expectFaceM = {  0.5, 1, 1.41421, 1, 0.5, 0.5, 1, 1, 0.5 };
    EXPECT_TRUE(floatArrayCompare(measures.getDataByDim<2>(), expectFaceM));
    std::vector<double> expectCellM = {  0.5, 0.5 };
    EXPECT_TRUE(floatArrayCompare(measures.getDataByDim<3>(), expectCellM));


    for (auto& cell : mesh.getCells()){
        EXPECT_EQ(cellDataLoad[cell].color, cellData[cell].color);
        EXPECT_TRUE(floatArrayCompare( cellDataLoad[cell].center, cellData[cell].center));
    }



    MeshDataContainer<MeshNativeType<3>::ElementType,3> types1(mesh, MeshNativeType<3>::POLYHEDRON);
    VTKMeshWriter<3, size_t, double> writer1;
    out3D.open("mesh_refine_1.vtk");
    EXPECT_TRUE(bool(out3D));
    writer.writeHeader(out3D, "test data");
    writer.writeToStream(out3D, mesh, types1);
    VTKMeshDataWriter<3>::writeToStream(out3D, cellData, writer1);
    out3D.close();

    UnstructuredMesh<3, size_t, double, 6> meshRefined;

    ifst.open("mesh_refine_1.vtk");
    EXPECT_TRUE(bool(ifst));
    reader.loadFromStream(ifst, meshRefined);
    MeshDataContainer<CellData<3>, 3> cellDataLoadRefined(mesh);
    VTKMeshDataReader<3, size_t>::readFromStream(ifst, cellDataLoad);
    ifst.close();

    for (auto& cell : meshRefined.getCells()){
        EXPECT_EQ(cellDataLoad[cell].color, cellData.getDataByPos<0>()[writer1.backwardCellIndexMapping[cell.getIndex()]].color);
        EXPECT_TRUE(floatArrayCompare( cellDataLoad[cell].center, cellData.getDataByPos<0>()[writer1.backwardCellIndexMapping[cell.getIndex()]].center));
    }
}


//#endif
#ifndef HAVE_GTEST
int main(){
    UnstructuredMesh2D_Functions_Test();
    UnstructuredMesh2DReadWrite();
    MeshRefineTest();
    UnstructuredMesh3D_Functions_Test();
}
#else
#include "UnitTests/main.h"
#endif

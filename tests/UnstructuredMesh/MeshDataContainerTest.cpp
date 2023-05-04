// Test of Traits class
#ifdef HAVE_GTEST
#include <gtest/gtest.h>
// #else
// #define TEST(_1,_2) void _1()
// #define EXPECT_TRUE(_1) (void)(_1 == true)
// #define EXPECT_FALSE(_1) (void)(_1 == false)
// #define EXPECT_EQ(_1,_2) (void)(_1 == _2)
// #endif
#include <list>
#include <map>
#include <array>
#include <string>
#include "GTMesh/UnstructuredMesh/UnstructuredMesh.h"
#include "MeshSetup.h"


TEST( MakeMeshDataContainerTest, basicTest )
{
    UnstructuredMesh<3, size_t, double, 6> mesh3;
    twoDeformedPrisms(mesh3);

    constexpr std::integer_sequence<unsigned int, 2> _1;

    MeshDataContainer<std::tuple<int, double, char, double>, 3,2,0> container(mesh3);


    MeshDataContainer<double, 3,3,1> cont(mesh3, 52.2);
    MeshDataContainer<double, 1,1> contAlloc;
    contAlloc.allocateData(cont, 42.2);

    EXPECT_EQ(cont.getDataByDim<1>(), std::vector<double>( 14, 52.2 ));
    EXPECT_EQ(cont.getDataByDim<1>().size(), contAlloc.getDataByDim<1>().size());
    EXPECT_EQ(cont.getDataByDim<1>(), cont[_1]);

    for(auto& c : container.getDataByDim<0>()) {
        c=42;
    }

    for(auto& v : mesh3.getVertices()){
        EXPECT_EQ(container.at(v), '*');
    }



    MeshDataContainer<std::tuple<int, double, char, int>, 3, 2, 0, 1> containerIni(mesh3,3, 42.15, 'a', 15);

    EXPECT_EQ(containerIni.getDataByPos<0>(), std::vector<int>(mesh3.getElements<3>().size(), 3));
    EXPECT_EQ(containerIni.getDataByPos<1>(), std::vector<double>(mesh3.getElements<2>().size(), 42.15));
    EXPECT_EQ(containerIni.getDataByPos<2>(), std::vector<char>(mesh3.getElements<0>().size(), 'a'));
    EXPECT_EQ(containerIni.getDataByPos<3>(), std::vector<int>(mesh3.getElements<1>().size(), 15));


    MeshDataContainer<std::tuple<int, double, char, int>, 3, 2, 0, 1> containerAssign;

    containerAssign = containerIni;

    contAlloc.getDataByDim<1>().clear();
    EXPECT_TRUE(contAlloc.getDataByDim<1>().empty());
    contAlloc.allocateData(containerAssign, 87);
    EXPECT_EQ(contAlloc.getDataByDim<1>(), std::vector<double>(mesh3.getElements<1>().size(), 87));
}


TEST(MakeMeshDataContainerTest, connections)
{
    UnstructuredMesh<3, size_t, double, 6> mesh3;
    twoDeformedPrisms(mesh3);

    MeshDataContainer<std::tuple<int, int>, 0, 3> meshData(mesh3, 0, 3);

    const auto conn03 = mesh3.connections<0, 3>();
    for (const auto &vert : mesh3.getElements<0>()) {
        for (const auto &elementIndex : conn03[vert]) {
            EXPECT_EQ(meshData[elementIndex], 3);
        }
    }

    const auto neighborhood = mesh3.neighborhood<0, 3>();
    for (const auto &vert : mesh3.getElements<0>()) {
        for (const auto &elementIndex : neighborhood[vert]) {
            EXPECT_EQ(meshData[elementIndex], 0);
        }
    }
}

#endif

#include "main.h"

// Test of Traits class
#ifdef HAVE_GTEST
#include <gtest/gtest.h>
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

    EXPECT_EQ(containerIni.getDataByPos<0>(), std::vector<int>(2, 3));
    EXPECT_EQ(containerIni.getDataByPos<1>(), std::vector<double>(2, 42.15));
    EXPECT_EQ(containerIni.getDataByPos<2>(), std::vector<char>(2, 'a'));
    EXPECT_EQ(containerIni.getDataByPos<3>(), std::vector<int>(2, 3));


    DBGMSG("assign test");
    MeshDataContainer<std::tuple<int, double, char, int>, 3, 2, 0, 1> containerAssign;

    containerAssign = containerIni;

    contAlloc.getDataByDim<1>().clear();
    EXPECT_TRUE(contAlloc.getDataByDim<1>().empty());
    contAlloc.allocateData(containerAssign, 87);
    EXPECT_EQ(contAlloc.getDataByDim<1>(), std::vector<double>(mesh3.getElements<1>().size(), 87));
}




#endif

#include "UnitTests/main.h"

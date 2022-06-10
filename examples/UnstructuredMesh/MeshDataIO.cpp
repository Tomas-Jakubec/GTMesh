#include <GTMesh/Debug/Debug.h>
#include <GTMesh/UnstructuredMesh/UnstructuredMesh.h>
#include <GTMesh/UnstructuredMesh/MeshDataContainer/MeshDataIO/VTKMeshDataReader.h>
#include <GTMesh/UnstructuredMesh/MeshDataContainer/MeshDataIO/VTKMeshDataWriter.h>
#include <GTMesh/Traits/Traits.h>
#include <fstream>

struct cellData{
    double invVol;
    size_t index;
};

MAKE_ATTRIBUTE_TRAIT(cellData, invVol, index);

int main () {
    UnstructuredMesh<3, size_t, double, 6> mesh;

    // load mesh from fpma file
    auto meshReaderFPMA = mesh.load("meshFile.fpma");

    MeshDataContainer<cellData, 3> meshData(mesh);

    mesh.initializeCenters();
    auto measures = mesh.computeElementMeasures();

    for (const auto& cell : mesh.getCells()){
        meshData[cell].invVol = 1.0 / measures[cell];
        meshData[cell].index = cell.getIndex();
    }

    std::ofstream oFile("meshFile-out.vtk");
    VTKMeshWriter<3, size_t, double> meshWriter;
    meshWriter.writeHeader(oFile, "mesh with data");
    meshWriter.writeToStream(oFile, mesh, meshReaderFPMA->getCellTypes());

    // append the data
    VTKMeshDataWriter<3>::writeToStream(oFile, meshData, meshWriter);

    // load the data back from meshFile-out
    UnstructuredMesh<3, size_t, double, 3> meshNew;
    auto meshReaderVTK = meshNew.load("meshFile-out.vtk");

    MeshDataContainer<cellData, 3> meshDataNew(meshNew);

    std::ifstream iFile("meshFile-out.vtk", std::ios::binary);

    VTKMeshDataReader<3, size_t>::readFromStream(iFile, meshDataNew);
    // Now the meshDataNew container contains the values stored in the meshFile-out.vtk
    DBGVAR(meshDataNew.getDataByDim<3>());
    return 0;
}

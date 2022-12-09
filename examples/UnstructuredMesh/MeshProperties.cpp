#include <iostream>
#include <fstream>
#include <unordered_set>
#include <unordered_map>
#include <GTMesh/Debug/Debug.h>
#include <GTMesh/UnstructuredMesh/UnstructuredMesh.h>
#include <GTMesh/UnstructuredMesh/MeshIO/MeshReader/FPMAMeshReader.h>

#include <GTMesh/UnstructuredMesh/MeshIO/MeshWriter/FPMAMeshWriter.h>
#include <GTMesh/UnstructuredMesh/MeshIO/MeshWriter/VTKMeshWriter.h>
using namespace std;

void meshSummary(const char* meshName){
    UnstructuredMesh<3, size_t, double> mesh;
    mesh.load(meshName);

    mesh.initializeCenters<METHOD_TESSELLATED>();
    mesh.setupBoundaryCells();

    mesh.setupBoundaryCellsCenters();

    auto measures = mesh.computeElementMeasures<METHOD_TESSELLATED>();

    auto normals = mesh.computeFaceNormals<METHOD_TESSELLATED>();
    auto cellsCenterDifference = computeCellsCenterDifference(mesh);

    // normalize the differences
    for (auto& diff : cellsCenterDifference.getDataByPos<0>()) diff /= diff.normEuclid();
    // the scalar product indicates how much parallel are the normal vectors and the connection lines
    // between cell centers
    auto scalarProduct = MeshDataContainer<double, 2>(mesh);
    for (const auto& face : mesh.getFaces()) {
        scalarProduct[face] = cellsCenterDifference[face] * normals[face];
    }

    double cellVolAvg = 0;
    for(auto cell : mesh.getCells()){
        cellVolAvg += measures[cell];
    }
    cellVolAvg /= mesh.getCells().size();

    double faceVolAvg = 0;
    for(auto face : mesh.getFaces()){
        faceVolAvg += measures[face];
    }
    faceVolAvg /= mesh.getFaces().size();

    // exports a json parsable using python
    DBGVAR_JSON(meshName, mesh.getCells().size(), mesh.getFaces().size(), cellVolAvg, faceVolAvg);

    DBGVAR_JSON(scalarProduct.getDataByPos<0>(), cellsCenterDifference.getDataByPos<0>());
}


int main(int argc, const char** argv)
{
    for(int i = 1; i < argc; i++){
        meshSummary(argv[i]);
    }
    return 0;
}

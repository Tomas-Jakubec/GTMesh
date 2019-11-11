#include <iostream>
//#define UNDEBUG
#include "../debug/Debug.h"
#include "UnstructuredMesh/UnstructuredMesh.h"
#include "UnstructuredMesh/MeshFunctions/MeshFunctions.h"
#include "UnstructuredMesh/MeshIO/MeshReader/VTKMeshReader.h"
#include "UnstructuredMesh/MeshIO/MeshWriter/VTKMeshWriter.h"
#include "UnstructuredMesh/MeshDataContainer/MeshDataIO/VTKMeshDataWriter.h"
#include "UnstructuredMesh/MeshDataContainer/MeshDataIO/VTKMeshDataReader.h"
#include "UnstructuredMesh/MeshIO/MeshReader/FPMAMeshReader.h"
#include "UnstructuredMesh/MeshIO/MeshWriter/FPMAMeshWriter.h"

#include "UnstructuredMesh/MeshDataContainer/MemberApproach.h"
#include <fstream>
#include <list>
using namespace std;



struct CellData {
    double invVol;

};


struct CompData {
    double T = 300;
};
MAKE_NAMED_ATTRIBUTE_TRAIT(CompData, "Temperature", T);

struct FaceData {
    double volOverDist;
};

void heatFlow(UnstructuredMesh<3,size_t, double, 12> mesh,
        MeshDataContainer<std::tuple<CellData, FaceData>, 3,2>& data,
        MeshDataContainer<CompData, 3>& compData,
        MeshDataContainer<double, 3>& outDeltas) {

    for (double& dT : outDeltas.getDataByPos<0>()){
        dT = 0;
    }

    for (auto& face : mesh.getFaces()) {
        double lT = 1000;
        bool lIn = false;
        if(!((face.getCellLeftIndex() & BOUNDARY_INDEX(size_t)) == BOUNDARY_INDEX(size_t))){
            lT = compData.getDataByDim<3>().at(face.getCellLeftIndex()).T;
            lIn = true;
        }

        double rT = 1000;
        bool rIn = false;
        if(!((face.getCellRightIndex() & BOUNDARY_INDEX(size_t)) == BOUNDARY_INDEX(size_t))){
            rT = compData.getDataByDim<3>().at(face.getCellRightIndex()).T;
            rIn = true;
        }


        double dT = data.at(face).volOverDist * (lT - rT);
        if (rIn)
            outDeltas.getDataByDim<3>().at(face.getCellRightIndex()) += dT;
        if (lIn)
            outDeltas.getDataByDim<3>().at(face.getCellLeftIndex()) -= dT;

    }
}

void explicitUpdate(UnstructuredMesh<3,size_t, double, 12> mesh,
                    MeshDataContainer<std::tuple<CellData, FaceData>, 3,2>& data,
                    MeshDataContainer<CompData, 3>& compData,
                    double tau,
                    double startTime,
                    double finalT){
    MeshDataContainer<double, 3> K1(mesh);
    double time = startTime;
    while (time < finalT) {

        heatFlow(mesh, data, compData, K1);

        for (size_t i = 0; i < mesh.getCells().size(); i++){
            auto& cellData = data.getDataByDim<3>().at(i);
            auto& _compData = compData.getDataByPos<0>().at(i);
            _compData.T += cellData.invVol * tau * K1.getDataByDim<3>().at(i);
        }
        cout << "time: " << time << "\r";
        time += tau;
    }
    DBGMSG("computation done");
}



void exportData(std::string filename,
                VTKMeshWriter<3, size_t, double>& writer,
                double time,
                UnstructuredMesh<3, size_t, double, 12>& mesh,
                MeshDataContainer<CompData, 3>& compData) {

    ofstream ofile(filename + "_" + to_string(time) + ".vtk");
    writer.writeHeader(ofile, "HC_test"s + to_string(time));
    writer.writeToStream(ofile, mesh, MeshDataContainer<MeshNativeType<3>::ElementType, 3>(mesh, MeshNativeType<3>::POLYHEDRON));

    VTKMeshDataWriter<3> dataWriter;
    dataWriter.writeToStream(ofile, compData, writer);

    ofile.close();
    DBGMSG("Data eported");

}


void testHeatConduction() {
    UnstructuredMesh<3, size_t, double, 12> mesh;
    FPMAMeshReader<3> reader;
    ifstream file("Poly_box.fpma");
    reader.loadFromStream(file, mesh);

    mesh.initializeCenters();

    mesh.setupBoundaryCells();



    mesh.setupBoundaryCellsCenters();

    auto measures = mesh.computeElementMeasures();

    auto dists = ComputeCellsDistance(mesh);

    MeshDataContainer<std::tuple<CellData, FaceData>, 3,2> meshData(mesh);

    MeshDataContainer<CompData, 3> compData(mesh);

    for (auto& face : mesh.getFaces()) {
        meshData.at(face).volOverDist = measures.at(face) / dists.at(face);
    }

    for (auto& cell : mesh.getCells()) {
        meshData.at(cell).invVol = 1.0 / measures.at(cell);
    }

    double startTime = 0.0, finalTime = 1e-2, tau = 1e-5;

    VTKMeshWriter<3, size_t, double> writer;

    exportData("Heatcond_test",
               writer,
               startTime,
               mesh,
               compData);


    explicitUpdate(mesh,
                   meshData,
                   compData,
                   tau,
                   startTime,
                   finalTime);

    exportData("Heatcond_test",
               writer,
               startTime,
               mesh,
               compData);


    startTime = finalTime;
    finalTime *= 2;
    explicitUpdate(mesh,
                   meshData,
                   compData,
                   tau,
                   startTime,
                   finalTime);

    exportData("Heatcond_test",
               writer,
               startTime,
               mesh,
               compData);



}



int main()
{
    testHeatConduction();
}

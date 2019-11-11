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
    double T = 300;
    double invVol;

};
struct FaceData {
    double volOverDist;
};

void heatFlow(UnstructuredMesh<3,size_t, double, 12> mesh,
        MeshDataContainer<std::tuple<CellData, FaceData>, 3,2>& data,
        MeshDataContainer<double, 3>& outDeltas) {

    for (double& dT : outDeltas.getDataByPos<0>()){
        dT = 0;
    }

    for (auto& face : mesh.getFaces()) {
        double lT = 1000;
        bool lIn = false;
        if(!((face.getCellLeftIndex() & BOUNDARY_INDEX(size_t)) == BOUNDARY_INDEX(size_t))){
            lT = data.getDataByDim<3>().at(face.getCellLeftIndex()).T;
            lIn = true;
        }

        double rT = 1000;
        bool rIn = false;
        if(!((face.getCellRightIndex() & BOUNDARY_INDEX(size_t)) == BOUNDARY_INDEX(size_t))){
            rT = data.getDataByDim<3>().at(face.getCellRightIndex()).T;
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
                    double tau,
                    double startTime,
                    double finalT){
    MeshDataContainer<double, 3> K1(mesh);
    double time = startTime;
    while (time < finalT) {

        heatFlow(mesh, data, K1);

        for (size_t i = 0; i < mesh.getCells().size(); i++){
            auto& cellData = data.getDataByDim<3>().at(i);
            cellData.T += cellData.invVol * tau * K1.getDataByDim<3>().at(i);
        }
        cout << "time: " << time << "\r";
        time += tau;
    }
    DBGMSG("computation done");
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

    MeshDataContainer<std::tuple<CellData, FaceData>, 3,2> compData(mesh);

    for (auto& face : mesh.getFaces()) {
        compData.at(face).volOverDist = measures.at(face) / dists.at(face);
    }

    for (auto& cell : mesh.getCells()) {
        compData.at(cell).invVol = 1.0 / measures.at(cell);
    }

    double startTime = 0.0, finalTime = 1e-2, tau = 1e-5;

    VTKMeshWriter<3, size_t, double> writer;

    exportData("Heatcond_test",
               writer,
               startTime,
               mesh,
               compData);


    explicitUpdate(mesh,
                   compData,
                   tau,
                   startTime,
                   finalTime);


    exportData("Heatcond_test",
               writer,
               finalTime,
               mesh,
               compData);

    startTime = finalTime;
    finalTime *= 2;
    explicitUpdate(mesh,
                   compData,
                   tau,
                   startTime,
                   finalTime);

    exportData("Heatcond_test",
               writer,
               finalTime,
               mesh,
               compData);

}



int main()
{
    testHeatConduction();
}

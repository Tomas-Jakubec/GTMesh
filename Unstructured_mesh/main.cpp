#include <iostream>
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

template <unsigned int ProblemDimension, typename Real>
class HeatCunduction {
public:
    UnstructuredMesh<ProblemDimension,size_t, double, 12> mesh;
    MeshDataContainer<std::tuple<CellData, FaceData>, 3,2> meshData;
    VTKMeshWriter<ProblemDimension, size_t, Real> writer;

    HeatCunduction(){
    }

    void loadMesh(const std::string& fileName) {
        FPMAMeshReader<3> reader;
        ifstream file(fileName);
        reader.loadFromStream(file, mesh);

        mesh.template initializeCenters<TESSELLATED>();

        mesh.setupBoundaryCells();



        mesh.setupBoundaryCellsCenters();


        auto measures = mesh.template computeElementMeasures<TESSELLATED>();

        auto dists = ComputeCellsDistance(mesh);


        meshData.alocateData(mesh);

        MeshDataContainer<CompData, 3> compData(mesh);

        for (auto& face : mesh.getFaces()) {
            meshData.at(face).volOverDist = measures.at(face) / dists.at(face);
        }

        for (auto& cell : mesh.getCells()) {
            meshData.at(cell).invVol = 1.0 / measures.at(cell);
        }
    }

    void calculateRHS(
            Real,//time is unused in this problem
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


            double dT = meshData.at(face).volOverDist * (lT - rT);
            if (rIn)
                outDeltas.getDataByDim<3>().at(face.getCellRightIndex()) += dT;
            if (lIn)
                outDeltas.getDataByDim<3>().at(face.getCellLeftIndex()) -= dT;

        }

        for (auto& cell : mesh.getCells()) {
            outDeltas.at(cell) *= meshData.at(cell).invVol;
        }
    }

    template<typename dataType>
    void exportData(double time,
                    MeshDataContainer<dataType, 3>& compData) {

        ofstream ofile("HeatCond"s + "_" + to_string(time) + ".vtk");
        writer.writeHeader(ofile, "HC_test"s + to_string(time));
        writer.writeToStream(ofile, mesh, MeshDataContainer<MeshNativeType<3>::ElementType, 3>(mesh, MeshNativeType<3>::POLYHEDRON));

        VTKMeshDataWriter<3> dataWriter;
        dataWriter.writeToStream(ofile, compData, writer);

        ofile.close();
        DBGMSG("Data eported");

    }
};


void RKMSolver(HeatCunduction<3,double>& problem,
                MeshDataContainer<CompData, 3>& compData,//x_ini
                double tau_ini,//tau_ini
                double startTime,
                double finalT,
                double delta){

    double tau = tau_ini;

    MeshDataContainer<CompData, 3> Ktemp(problem.mesh);

    MeshDataContainer<double, 3> K1(problem.mesh);
    MeshDataContainer<double, 3> K2(problem.mesh);
    MeshDataContainer<double, 3> K3(problem.mesh);
    MeshDataContainer<double, 3> K4(problem.mesh);
    MeshDataContainer<double, 3> K5(problem.mesh);
    double time = startTime;
    bool run = true;
    while (run) {
        if (time + tau > finalT) {
            tau = finalT - time;
            run = false;
        }

        problem.calculateRHS(time, compData, K1);

        for (size_t i = 0; i < Ktemp.getDataByPos<0>().size(); i++){
            Ktemp.getDataByPos<0>().at(i).T = compData.getDataByDim<3>().at(i).T + (tau * (1.0 / 3.0) * K1.getDataByDim<3>().at(i));
        }

        problem.calculateRHS(time, Ktemp, K2);

        for (size_t i = 0; i < Ktemp.getDataByPos<0>().size(); i++){
            Ktemp.getDataByPos<0>().at(i).T = compData.getDataByDim<3>().at(i).T + (tau * (1.0 / 6.0) * (K1.getDataByDim<3>().at(i) + K2.getDataByDim<3>().at(i)));
        }


        problem.calculateRHS(time, Ktemp, K3);

        for (size_t i = 0; i < Ktemp.getDataByPos<0>().size(); i++){
            Ktemp.getDataByPos<0>().at(i).T = compData.getDataByDim<3>().at(i).T + (tau * (0.125 * K1.getDataByDim<3>().at(i) + 0.375 * K3.getDataByDim<3>().at(i)));
        }

        problem.calculateRHS(time, Ktemp, K4);

        for (size_t i = 0; i < Ktemp.getDataByPos<0>().size(); i++){
            Ktemp.getDataByPos<0>().at(i).T = compData.getDataByDim<3>().at(i).T + (tau * ((0.5 * K1.getDataByDim<3>().at(i)) - (1.5 * K3.getDataByDim<3>().at(i)) + (2.0 * K4.getDataByPos<0>().at(i))));
        }

        problem.calculateRHS(time, Ktemp, K5);

        double error = 0.0;

        for (size_t i = 0; i < K4.getDataByPos<0>().size(); i++){
            double tmpE = abs(0.2 * K1.getDataByPos<0>().at(i) - 0.9 * K3.getDataByPos<0>().at(i) +
                    0.8 * K4.getDataByPos<0>().at(i) - 0.1 * K5.getDataByPos<0>().at(i));
            if (tmpE > error) {
                error = tmpE;
            }
        }

        error *= tau * (1.0 / 3.0);
        if (error < delta) {
            for (size_t i = 0; i < K4.getDataByPos<0>().size(); i++){
                auto& _compData = compData.getDataByPos<0>().at(i);
                _compData.T += tau * (1.0 / 6.0) * (((K1.getDataByDim<3>().at(i) + K5.getDataByDim<3>().at(i))) + (4.0 * K4.getDataByPos<0>().at(i)));
            }
            time += tau;
            tau *= 1.05;
            cout << "time: " << time << "\r";
        } else {
            tau *= std::pow(0.2, delta/error) * 0.8;
            run = true;
        }

    }
    DBGMSG("computation done");
}



struct Index {
    size_t index;
};

MAKE_ATTRIBUTE_TRAIT(Index, index);

void meshAdmisibility(const std::string& meshFileName) {
    UnstructuredMesh<3, size_t, double, 12> mesh;
    FPMAMeshReader<3> reader;
    ifstream file(meshFileName);
    reader.loadFromStream(file, mesh);

    mesh.initializeCenters<TESSELLATED>();

    mesh.setupBoundaryCells();



    mesh.setupBoundaryCellsCenters();


    auto faceNormals = mesh.computeFaceNormals<TESSELLATED>();

    MeshDataContainer<Vector<3, double>, 2> centerLines(mesh);

    double max = -1, min = 2, avg = 0;
    MeshDataContainer<double, 2> faceAdmisibility(mesh);

    auto faceVert = MeshConnections<2,0, ORDER_ORIGINAL>::connections(mesh);

    int cnt = 0;
    for( auto& face : mesh.getFaces()) {
        //DBGVAR(faceIndex);
        bool badFace = false;
        for (size_t i = 1; i < faceVert.at(face).size(); i++) {
            size_t index = faceVert.at(face).at(i);
            Vertex<3,double> centerToVert = mesh.getVertices().at(faceVert.at(face).at(0)) - mesh.getVertices().at(index);
            double scalarProduct = faceNormals.at(face) * (Vector<3,double>{centerToVert[0],centerToVert[1],centerToVert[2]});
            //HTMLDBGCOND(scalarProduct > 0.01, faceNormals.at(face), centerToVert, scalarProduct);
            if (abs(scalarProduct) > 1e-4) badFace = true; // the face is not planar
        }
        if (badFace) cnt++;
    }
    DBGVAR(cnt, mesh.getFaces().size());

    for (auto& face : mesh.getFaces()) {

        auto& cellLeft = (face.getCellLeftIndex() & BOUNDARY_INDEX(size_t)) == BOUNDARY_INDEX(size_t)?
                  mesh.getBoundaryCells().at(face.getCellLeftIndex()&EXTRACTING_INDEX(size_t)):
                  mesh.getCells().at(face.getCellLeftIndex());
        auto& cellRight = (face.getCellRightIndex() & BOUNDARY_INDEX(size_t)) == BOUNDARY_INDEX(size_t)?
                  mesh.getBoundaryCells().at(face.getCellRightIndex()&EXTRACTING_INDEX(size_t)):
                  mesh.getCells().at(face.getCellRightIndex());


        Vertex<3,double> line = cellRight.getCenter() - cellLeft.getCenter();
        centerLines.at(face) = {line[0],line[1],line[2]};
        centerLines.at(face) /= centerLines.at(face).normEukleid();

        faceAdmisibility.at(face) = faceNormals.at(face) * centerLines.at(face);
        if (faceAdmisibility.at(face) > max) {
            max = faceAdmisibility.at(face);
        }
        if (faceAdmisibility.at(face) < min) {
            min = faceAdmisibility.at(face);
        }
        avg += faceAdmisibility.at(face);
        //DBGVARCOND((face.getCellLeftIndex() == 691 || face.getCellLeftIndex() == 699) && (face.getCellRightIndex() == 691 || face.getCellRightIndex() == 699), face.getIndex());
        //DBGVARCOND_HTML(faceAdmisibility.at(face) < 0.0, faceNormals.at(face), centerLines.at(face), face.getIndex(), face.getCenter(), face.getCellLeftIndex(), cellLeft.getCenter(), face.getCellRightIndex(), cellRight.getCenter(), faceAdmisibility.at(face));
    }
    avg /= mesh.getFaces().size();
    DBGVAR(max, min, avg);

    VTKMeshWriter<3, size_t, double> writer;
    MeshDataContainer<Index, 3> indexes(mesh);
    for (size_t i = 0; i < mesh.getCells().size(); i++) {
        indexes.getDataByDim<3>().at(i).index = i;
    }
    //exportData("Admisibility", writer, 0, mesh, indexes);
}






void testHeatConduction1() {

    HeatCunduction<3,double> hcProblem;

    hcProblem.loadMesh("Poly_2_5_level2.fpma");

    double startTime = 0.0, finalTime = 0.5 , tau = 1e-4;

    VTKMeshWriter<3, size_t, double> writer;

    MeshDataContainer<CompData, 3> compData(hcProblem.mesh);

    hcProblem.exportData(startTime, compData);

    RKMSolver(hcProblem, compData,
                   tau,
                   startTime,
                   finalTime,
                   1e-2);

    hcProblem.exportData(finalTime, compData);



    startTime = finalTime;
    finalTime *= 2;
    RKMSolver(hcProblem, compData,
                   tau,
                   startTime,
                   finalTime,
                   1e-2);


    hcProblem.exportData(finalTime, compData);




}

int main()
{
    testHeatConduction1();
    //testHeatConduction();
    //meshAdmisibility("Poly_2_5_level1.fpma");

}

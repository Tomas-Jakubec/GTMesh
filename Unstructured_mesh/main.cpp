#include <iostream>
//#define UNDEBUG
#define CONSOLE_COLORED_OUTPUT
#include "../src/Debug/Debug.h"
#include "../src/UnstructuredMesh/UnstructuredMesh.h"
#include "../src/UnstructuredMesh/MeshFunctions/MeshFunctions.h"
#include "../src/UnstructuredMesh/MeshIO/MeshReader/VTKMeshReader.h"
#include "../src/UnstructuredMesh/MeshIO/MeshWriter/VTKMeshWriter.h"
#include "../src/UnstructuredMesh/MeshDataContainer/MeshDataIO/VTKMeshDataWriter.h"
#include "../src/UnstructuredMesh/MeshDataContainer/MeshDataIO/VTKMeshDataReader.h"
#include "../src/UnstructuredMesh/MeshIO/MeshReader/FPMAMeshReader.h"
#include "../src/UnstructuredMesh/MeshIO/MeshWriter/FPMAMeshWriter.h"

#include "../src/Traits/Traits.h"
#include "../src/Traits/TraitsAlgorithm/TraitsAlgorithm.h"
#include <fstream>
#include <list>
#include <chrono>
using namespace std;



struct CellData {
    double invVol;

};


struct CompData {
    double T;
    //CompData(const CompData&) = default;
    //CompData(CompData&&) = default;
    CompData(double _T = 0){T = _T;}

};
MAKE_NAMED_ATTRIBUTE_TRAIT(CompData, "Temperature", T);
MAKE_ATTRIBUTE_TRAIT_ARITHMETIC(CompData, T);

struct FaceData {
    double volOverDist;
};

template <unsigned int ProblemDimension, typename Real>
class HeatCunduction {
public:
    using MeshType = UnstructuredMesh<ProblemDimension,size_t, double, 12>;
    using ResultType = CompData;

public:
    MeshType mesh;
    MeshDataContainer<std::tuple<CellData, FaceData>, 3,2> meshData;
    VTKMeshWriter<ProblemDimension, size_t, Real> writer;

    HeatCunduction(){
    }

    void loadMesh(const std::string& fileName) {
        FPMAMeshReader<MeshType::meshDimension()> reader;
        ifstream file(fileName);
        reader.loadFromStream(file, mesh);

        mesh.template initializeCenters<TESSELLATED>();

        mesh.setupBoundaryCells();



        mesh.setupBoundaryCellsCenters();


        auto measures = mesh.template computeElementMeasures<TESSELLATED>();

        auto dists = ComputeCellsDistance(mesh);


        meshData.allocateData(mesh);

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
            MeshDataContainer<CompData, ProblemDimension>& compData,
            MeshDataContainer<CompData, ProblemDimension>& outDeltas) {

        for (CompData& dT : outDeltas.template getDataByPos<0>()){
            dT.T = 0;
        }

        for (auto& face : mesh.getFaces()) {
            CompData lT(1000);
            bool lIn = false;
            if(!((face.getCellLeftIndex() & BOUNDARY_INDEX(size_t)) == BOUNDARY_INDEX(size_t))){
                lT = compData.template getDataByDim<3>().at(face.getCellLeftIndex());
                lIn = true;
            }

            CompData rT(1000);
            bool rIn = false;
            if(!((face.getCellRightIndex() & BOUNDARY_INDEX(size_t)) == BOUNDARY_INDEX(size_t))){
                rT = compData.template getDataByDim<3>().at(face.getCellRightIndex());
                rIn = true;
            }


            CompData dT = meshData.at(face).volOverDist * (lT - rT);
            if (rIn)
                outDeltas.template getDataByDim<3>().at(face.getCellRightIndex()) += dT;
            if (lIn)
                outDeltas.template getDataByDim<3>().at(face.getCellLeftIndex()) -= dT;

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

struct watch{

    std::chrono::high_resolution_clock clock;
    std::chrono::time_point<std::chrono::high_resolution_clock> timeStart;
    long long duration;
    long long deviation;
    unsigned long lapCnt;

    watch(){
        reset();
    }

    void reset(){
        duration = 0;
        deviation = 0;
        lapCnt = 0;
    }

    void start(){
        timeStart = clock.now();
    }

    void lap(){
        auto lapTime = clock.now() - timeStart;
        duration += lapTime.count();
        deviation += lapTime.count() * lapTime.count();
        lapCnt++;
    }


    struct stats{
        double avg;
        double dev;
        unsigned long lapCnt;
    };

    stats getResult(){
        double avg = double(duration)/lapCnt;
        double dev = sqrt(((deviation) - lapCnt * pow(avg,2)) / (lapCnt - 1));
        return stats{avg, dev, lapCnt};
    }
};

MAKE_ATTRIBUTE_TRAIT_IO(watch::stats, avg, dev, lapCnt);

watch
wK1,
wK2,
wK3,
wK4,
wError;

template <typename Problem, typename = typename std::enable_if<HasDefaultArithmeticTraits<typename Problem::ResultType>::value>::type>
void RKMSolver(
        Problem& problem,
        MeshDataContainer<typename Problem::ResultType, Problem::MeshType::meshDimension()>& compData,//x_ini
        double tau_ini,//tau_ini
        double startTime,
        double finalT,
        double delta
        )
{
    constexpr unsigned int MeshDimension = Problem::MeshType::meshDimension();
    double tau = tau_ini;

    MeshDataContainer<typename Problem::ResultType, MeshDimension> Ktemp;
    Ktemp.allocateData(compData);

    MeshDataContainer<typename Problem::ResultType, MeshDimension> K1(compData);
    MeshDataContainer<typename Problem::ResultType, MeshDimension> K2(compData);
    MeshDataContainer<typename Problem::ResultType, MeshDimension> K3(compData);
    MeshDataContainer<typename Problem::ResultType, MeshDimension> K4(compData);
    MeshDataContainer<typename Problem::ResultType, MeshDimension> K5(compData);
    double time = startTime;
    bool run = true;
    while (run) {
        if (time + tau > finalT) {
            tau = finalT - time;
            run = false;
        }
        problem.calculateRHS(time, compData, K1);
wK1.start();
        for (size_t i = 0; i < Ktemp.template getDataByPos<0>().size(); i++){
            Ktemp.template getDataByPos<0>().at(i) = compData.template getDataByDim<MeshDimension>().at(i) + (tau * (1.0 / 3.0) * K1.template getDataByDim<MeshDimension>().at(i));
        }
wK1.lap();

        problem.calculateRHS(time, Ktemp, K2);
wK2.start();
        for (size_t i = 0; i < Ktemp.template getDataByPos<0>().size(); i++){
            Ktemp.template getDataByPos<0>().at(i) = compData.template getDataByDim<MeshDimension>().at(i) + (tau * (1.0 / 6.0) * (K1.template getDataByDim<MeshDimension>().at(i) + K2.template getDataByDim<MeshDimension>().at(i)));
        }
wK2.lap();

        problem.calculateRHS(time, Ktemp, K3);
wK3.start();
        for (size_t i = 0; i < Ktemp.template getDataByPos<0>().size(); i++){
            Ktemp.template getDataByPos<0>().at(i) = compData.template getDataByDim<MeshDimension>().at(i) + (tau * (0.125 * K1.template getDataByDim<MeshDimension>().at(i) + 0.375 * K3.template getDataByDim<MeshDimension>().at(i)));
        }
wK3.lap();
        problem.calculateRHS(time, Ktemp, K4);
wK4.start();
        for (size_t i = 0; i < Ktemp.template getDataByPos<0>().size(); i++){
            Ktemp.template getDataByPos<0>().at(i) = compData.template getDataByDim<MeshDimension>().at(i) + (tau * ((0.5 * K1.template getDataByDim<MeshDimension>().at(i)) - (1.5 * K3.template getDataByDim<MeshDimension>().at(i)) + (2.0 * K4.template getDataByPos<0>().at(i))));
        }

wK4.lap();
        problem.calculateRHS(time, Ktemp, K5);


        double error = 0.0;
wError.start();
        for (size_t i = 0; i < K4.template getDataByPos<0>().size(); i++){
            double tmpE = max(abs(0.2 * K1.template getDataByPos<0>().at(i) - 0.9 * K3.template getDataByPos<0>().at(i) +
                    0.8 * K4.template getDataByPos<0>().at(i) - 0.1 * K5.template getDataByPos<0>().at(i)));

            if (tmpE > error) {
                error = tmpE;
            }
        }
wError.lap();
//DBGVAR_CSV(error);

        error *= tau * (1.0 / 3.0);
        if (error < delta) {
            for (size_t i = 0; i < K4.template getDataByPos<0>().size(); i++){
                auto& _compData = compData.template getDataByPos<0>().at(i);
                _compData += tau * (1.0 / 6.0) * (((K1.template getDataByDim<MeshDimension>().at(i) + K5.template getDataByDim<MeshDimension>().at(i))) + (4.0 * K4.template getDataByPos<0>().at(i)));
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

    double startTime = 0.0, finalTime = 0.1 , tau = 1e-4;

    VTKMeshWriter<3, size_t, double> writer;

    MeshDataContainer<CompData, 3> compData(hcProblem.mesh, CompData(300));

    hcProblem.exportData(startTime, compData);

    RKMSolver(hcProblem, compData,
                   tau,
                   startTime,
                   finalTime,
                   1);
DBGVAR(wK1.getResult(),wK2.getResult(),wK3.getResult(),wK4.getResult(),wError.getResult());
wK1.reset();wK2.reset();wK3.reset();wK4.reset();wError.reset();

    hcProblem.exportData(finalTime, compData);



    startTime = finalTime;
    finalTime *= 2;
    RKMSolver(hcProblem, compData,
                   tau,
                   startTime,
                   finalTime,
                   1);


    hcProblem.exportData(finalTime, compData);

DBGVAR(wK1.getResult(),wK2.getResult(),wK3.getResult(),wK4.getResult(),wError.getResult());



}



int main()
{
    testHeatConduction1();
    //testHeatConduction();
    //meshAdmisibility("Poly_2_5_level1.fpma");

}

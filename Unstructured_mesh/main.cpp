#include <iostream>
//#define UNDEBUG
#define CONSOLE_COLORED_OUTPUT
#include "../src/Debug/Debug.h"
#include "multiphaseflow.h"

#include <fstream>
#include <list>
#include <chrono>

using namespace std;



// struct CellData {
//     double invVol;
//
// };
//
//
// struct CompData {
//     double T;
//     //CompData(const CompData&) = default;
//     //CompData(CompData&&) = default;
//     CompData(double _T = 0){T = _T;}
//
// };
// MAKE_NAMED_ATTRIBUTE_TRAIT(CompData, "Temperature", T);
// MAKE_ATTRIBUTE_TRAIT_ARITHMETIC(CompData, T);
//
// struct FaceData {
//     double volOverDist;
// };
//
// template <unsigned int ProblemDimension, typename Real>
// class HeatCunduction {
// public:
//     using MeshType = UnstructuredMesh<ProblemDimension,size_t, double, 12>;
//     using ResultType = CompData;
//
// public:
//     MeshType mesh;
//     MeshDataContainer<std::tuple<CellData, FaceData>, 3,2> meshData;
//     VTKMeshWriter<ProblemDimension, size_t, Real> writer;
//
//     HeatCunduction(){
//     }
//
//     void loadMesh(const std::string& fileName) {
//         FPMAMeshReader<MeshType::meshDimension()> reader;
//         ifstream file(fileName);
//         reader.loadFromStream(file, mesh);
//
//         mesh.template initializeCenters<TESSELLATED>();
//
//         mesh.setupBoundaryCells();
//
//
//
//         mesh.setupBoundaryCellsCenters();
//
//
//         auto measures = mesh.template computeElementMeasures<TESSELLATED>();
//
//         auto dists = ComputeCellsDistance(mesh);
//
//
//         meshData.allocateData(mesh);
//
//         MeshDataContainer<CompData, 3> compData(mesh);
//
//         for (auto& face : mesh.getFaces()) {
//             meshData.at(face).volOverDist = measures.at(face) / dists.at(face);
//         }
//
//         for (auto& cell : mesh.getCells()) {
//             meshData.at(cell).invVol = 1.0 / measures.at(cell);
//         }
//     }
//
//     void calculateRHS(
//             Real,//time is unused in this problem
//             MeshDataContainer<CompData, ProblemDimension>& compData,
//             MeshDataContainer<CompData, ProblemDimension>& outDeltas) {
//
//         for (CompData& dT : outDeltas.template getDataByPos<0>()){
//             dT.T = 0;
//         }
//
//         for (auto& face : mesh.getFaces()) {
//             CompData lT(1000);
//             bool lIn = false;
//             if(!((face.getCellLeftIndex() & BOUNDARY_INDEX(size_t)) == BOUNDARY_INDEX(size_t))){
//                 lT = compData.template getDataByDim<3>().at(face.getCellLeftIndex());
//                 lIn = true;
//             }
//
//             CompData rT(1000);
//             bool rIn = false;
//             if(!((face.getCellRightIndex() & BOUNDARY_INDEX(size_t)) == BOUNDARY_INDEX(size_t))){
//                 rT = compData.template getDataByDim<3>().at(face.getCellRightIndex());
//                 rIn = true;
//             }
//
//
//             CompData dT = meshData.at(face).volOverDist * (lT - rT);
//             if (rIn)
//                 outDeltas.template getDataByDim<3>().at(face.getCellRightIndex()) += dT;
//             if (lIn)
//                 outDeltas.template getDataByDim<3>().at(face.getCellLeftIndex()) -= dT;
//
//         }
//
//         for (auto& cell : mesh.getCells()) {
//             outDeltas.at(cell) *= meshData.at(cell).invVol;
//         }
//     }
//
//     template<typename dataType>
//     void exportData(double time,
//                     MeshDataContainer<dataType, 3>& compData) {
//
//         using std::string_literals;
//
//         ofstream ofile("HeatCond"s + "_" + to_string(time) + ".vtk");
//         writer.writeHeader(ofile, "HC_test"s + to_string(time));
//         writer.writeToStream(ofile, mesh, MeshDataContainer<MeshNativeType<3>::ElementType, 3>(mesh, MeshNativeType<3>::POLYHEDRON));
//
//         VTKMeshDataWriter<3> dataWriter;
//         dataWriter.writeToStream(ofile, compData, writer);
//
//         ofile.close();
//         DBGMSG("Data eported");
//
//     }
// };

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

        error *= tau * (1.0 / 3.0);
        if (error < delta) {
            for (size_t i = 0; i < K4.template getDataByPos<0>().size(); i++){
                auto& _compData = compData.template getDataByPos<0>().at(i);
                _compData += tau * (1.0 / 6.0) * (((K1.template getDataByDim<MeshDimension>().at(i) + K5.template getDataByDim<MeshDimension>().at(i))) + (4.0 * K4.template getDataByPos<0>().at(i)));
            }
            time += tau;
            tau *= 1.005;
            //cout << "time: " << time << "\r";
        } else {
            tau *= std::pow(0.2, delta/error) * 0.8;
            run = true;
        }

    }
    DBGMSG("computation done");
}





MAKE_ATTRIBUTE_TRAIT(CellData, invVolume);
MAKE_ATTRIBUTE_TRAIT(EdgeData, n,LengthOverDist, Length,LeftCellKoef, RightCellKoef);
MAKE_ATTRIBUTE_TRAIT(PointData, u_s, u_g, PointType, cellKoef);

void testHeatConduction1() {

    MultiphaseFlow mpf;

    mpf.setupMeshData("Boiler.vtk");

    FlowData::R_spec = 287;
    FlowData::T = 300;
    FlowData::rho_s = 1700;

    mpf.myu = 1e5;
    mpf.myu_s = 1.5;
    mpf.d_s = 0.06;
    mpf.outFlow.setPressure(1e5);
    //mpf.outFlow.rho_g = 1.3;
    //mpf.outFlow.eps_g = 1;
    mpf.outFlow.eps_s = 0;
    mpf.R_spec = 287;
    mpf.T = 300;
    mpf.artifitialDisspation = 1;
    mpf.phi_s = 1;
    mpf.rho_s = 1700;
    mpf.inFlow_eps_g = 1;
    mpf.inFlow_eps_s = 0;
    mpf.inFlow_u_g = {0.0,15};
    mpf.inFlow_u_s = {0, 0};


    FlowData ini;
//    ini.eps_g = 1;
    ini.eps_s = 0;
    ini.rho_g = 1.3;
    ini.setPressure(1e5);
    ini.p_g = {};
    //ini.setVelocityGas({0, 0});
    ini.p_s = {0, 0};


    MeshDataContainer<FlowData, 2> compData(mpf.mesh, ini);
    for (auto& cell : mpf.mesh.getCells()){
        if(cell.getCenter()[1] > 7 && cell.getCenter()[1] < 9){
            compData.at(cell).eps_s = 0.1;
        }
    }
    //MeshDataContainer<FlowData, 2> result(compData);

    //mpf.ActualizePointData(compData);

    //mpf.calculateRHS(0.0, compData, result);

    mpf.exportData(0.0, compData);
    for (double t = 0; t < 1; t += 0.05){
        RKMSolver(mpf, compData, 1e-4, t, t + 0.05, 1);
        mpf.exportData(t + 0.05, compData);
    }

/*
    HeatCunduction<3,double> hcProblem;
    hcProblem.loadMesh("Poly_2_5_level2.fpma");

    double startTime = 0.0, finalTime = 0.1 , tau = 1e-4;

    VTKMeshWriter<3, size_t, double> writer;

    MeshDataContainer<CompData, 3> compData(hcProblem.mesh, CompData(300));

//    hcProblem.exportData(startTime, compData);

    RKMSolver(hcProblem, compData,
                   tau,
                   startTime,
                   finalTime,
                   1);
DBGVAR(wK1.getResult(),wK2.getResult(),wK3.getResult(),wK4.getResult(),wError.getResult());
wK1.reset();wK2.reset();wK3.reset();wK4.reset();wError.reset();

//    hcProblem.exportData(finalTime, compData);



    startTime = finalTime;
    finalTime *= 2;
    RKMSolver(hcProblem, compData,
                   tau,
                   startTime,
                   finalTime,
                   1);


//    hcProblem.exportData(finalTime, compData);

DBGVAR(wK1.getResult(),wK2.getResult(),wK3.getResult(),wK4.getResult(),wError.getResult());
*/


}



int main()
{
    testHeatConduction1();
    //testHeatConduction();
    //meshAdmisibility("Poly_2_5_level1.fpma");

}

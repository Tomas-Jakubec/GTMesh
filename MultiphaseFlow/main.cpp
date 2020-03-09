#include <iostream>
//#define UNDEBUG
//#define CONSOLE_COLORED_OUTPUT
#include "../src/Debug/Debug.h"
#include "multiphaseflow.h"

#include <fstream>
#include <list>
#include <chrono>

using namespace std;



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
wK1;

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

    DBGMSG("RKM_start", run);

    while (run) {
#pragma omp single
        {
            wK1.start();
            if (time + tau > finalT) {
                tau = finalT - time;
                run = false;
            }
        }
#pragma omp parallel
        {
        problem.calculateRHS(time, compData, K1);

        #pragma omp for
        for (size_t i = 0; i < Ktemp.template getDataByPos<0>().size(); i++){
            Ktemp.template getDataByPos<0>().at(i) = compData.template getDataByDim<MeshDimension>().at(i) + (tau * (1.0 / 3.0) * K1.template getDataByDim<MeshDimension>().at(i));
        }


        problem.calculateRHS(time, Ktemp, K2);

        #pragma omp for
        for (size_t i = 0; i < Ktemp.template getDataByPos<0>().size(); i++){
            Ktemp.template getDataByPos<0>().at(i) = compData.template getDataByDim<MeshDimension>().at(i) + (tau * (1.0 / 6.0) * (K1.template getDataByDim<MeshDimension>().at(i) + K2.template getDataByDim<MeshDimension>().at(i)));
        }


        problem.calculateRHS(time, Ktemp, K3);


        #pragma omp for
        for (size_t i = 0; i < Ktemp.template getDataByPos<0>().size(); i++){
            Ktemp.template getDataByPos<0>().at(i) = compData.template getDataByDim<MeshDimension>().at(i) + (tau * (0.125 * K1.template getDataByDim<MeshDimension>().at(i) + 0.375 * K3.template getDataByDim<MeshDimension>().at(i)));
        }


        problem.calculateRHS(time, Ktemp, K4);

        #pragma omp for
        for (size_t i = 0; i < Ktemp.template getDataByPos<0>().size(); i++){
            Ktemp.template getDataByPos<0>().at(i) = compData.template getDataByDim<MeshDimension>().at(i) + (tau * ((0.5 * K1.template getDataByDim<MeshDimension>().at(i)) - (1.5 * K3.template getDataByDim<MeshDimension>().at(i)) + (2.0 * K4.template getDataByPos<0>().at(i))));
        }


        problem.calculateRHS(time, Ktemp, K5);
        }// end of parallel section

        double error = 0.0;

#pragma omp parallel
        #pragma omp for reduction(max : error)
        for (size_t i = 0; i < K4.template getDataByPos<0>().size(); i++){
            double tmpE = max(abs(0.2 * K1.template getDataByPos<0>().at(i) - 0.9 * K3.template getDataByPos<0>().at(i) +
                    0.8 * K4.template getDataByPos<0>().at(i) - 0.1 * K5.template getDataByPos<0>().at(i)));

            if (tmpE > error) {
                error = tmpE;
            }
        }



        error *= tau * (1.0 / 3.0);

        if (error < delta) {
            #pragma omp parallel
            #pragma omp for
            for (size_t i = 0; i < K4.template getDataByPos<0>().size(); i++){
                auto& _compData = compData.template getDataByPos<0>().at(i);
                _compData += tau * (1.0 / 6.0) * (((K1.template getDataByDim<MeshDimension>().at(i) + K5.template getDataByDim<MeshDimension>().at(i))) + (4.0 * K4.template getDataByPos<0>().at(i)));
            }

            time += tau;
            tau *= 1.005;
            cout << "time: " << time << "\r";
            wK1.lap();

        } else {

            wK1.lap();
            tau *= std::pow(0.2, delta/error) * 0.8;
            run = true;

        }

    }

    DBGMSG("compuatation done");
}


template <typename Problem, typename = typename std::enable_if<HasDefaultArithmeticTraits<typename Problem::ResultType>::value>::type>
void EulerSolver(
        Problem& problem,
        MeshDataContainer<typename Problem::ResultType, Problem::MeshType::meshDimension()>& compData,//x_ini
        double tau_ini,//tau_ini
        double startTime,
        double finalT
        )
{
    constexpr unsigned int MeshDimension = Problem::MeshType::meshDimension();
    double tau = tau_ini;


    MeshDataContainer<typename Problem::ResultType, MeshDimension> K1(compData);
    double time = startTime;
    bool run = true;
    while (run) {
        if (time + tau >= finalT) {
            tau = finalT - time;
            run = false;
        }
        problem.calculateRHS(time, compData, K1);


        for (size_t i = 0; i < compData.template getDataByPos<0>().size(); i++){

            compData.template getDataByPos<0>().at(i) += tau * K1.template getDataByDim<MeshDimension>().at(i);
        }

        time += tau;
    }
    DBGMSG("compuatation done");
}



void MultiphaseFlowCalculation() {
    constexpr unsigned int ProblemDim = 3;

    using MPFType = MultiphaseFlow<ProblemDim>;

    MPFType mpf;

    // setup constants
    MPFType::ResultType::R_spec = 287;
    MPFType::ResultType::T = 300;
    MPFType::ResultType::rho_s = 1700;


    mpf.artifitialDisspation = 0.01;
    mpf.R_spec = 287;
    mpf.myu = 1e-5;
    mpf.rho_s = 1700;
    mpf.myu_s = 0.5;//1.5;
    mpf.d_s = 0.002;//0.06;
    mpf.phi_s = 1;

    mpf.T = 300;

    mpf.setupMeshData("stack.fpma");



    mpf.outFlow.setPressure(1e5);
    mpf.outFlow.eps_s = 0;

    mpf.inFlow_eps_g = 1;
    mpf.inFlow_eps_s = 0;
    mpf.inFlow_u_g = {};
    mpf.inFlow_u_g[ProblemDim - 1] = 15;
    mpf.inFlow_u_s = {};


    FlowData<ProblemDim> ini;
    ini.eps_s = 0;
    ini.setPressure(1e5);
    ini.p_g = {};
    ini.p_s = {};

    MeshDataContainer<FlowData<ProblemDim>, ProblemDim> compData(mpf.mesh, ini);

/*
    for (auto& cell : mpf.mesh.getCells()){
        if(
                cell.getCenter()[1] > 1.0 && cell.getCenter()[1] < 7.0 &&
                cell.getCenter()[0] > 0 && cell.getCenter()[0] < 10
          ){
            compData.at(cell).eps_s = 0.2;
        }
    }
*/

    for (auto& cell : mpf.mesh.getCells()){
        if(
               cell.getCenter()[2] > -1.0 && cell.getCenter()[2] < 2.0
//           cell.getCenter()[1] > -0.5 && cell.getCenter()[1] < 0.5 &&
//           cell.getCenter()[0] > -0.5 && cell.getCenter()[0] < 0.5
        ) {
            compData.at(cell).eps_s = 0.5;
        }
    }


    mpf.exportData(0.0, compData);

    double exportStep = 1e-1;
    for (double t = 0; t < 100*exportStep; t += exportStep){


        RKMSolver(mpf, compData, 1e-3, t, t + exportStep, 1e-2);

        //EulerSolver(mpf, compData, 1e-5, t, t + exportStep);
        mpf.exportData((t + exportStep)*1e-2, compData);

        DBGVAR(t + exportStep, wK1.getResult());
    }




}

void atEnd(){
    DBGVAR(wK1.getResult());
}

int main()
{
    atexit(atEnd);
    MultiphaseFlowCalculation();

    //testHeatConduction();
    //meshAdmisibility("Poly_2_5_level1.fpma");

}

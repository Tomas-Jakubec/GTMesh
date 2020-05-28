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

    while (true) {
        #pragma omp single
        {
            wK1.start();
            if (time + tau > finalT) {
                tau = finalT - time;
                run = false;
            } else {
                run = true;
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

            if (!run) {
                break;
            }
            if (error == 0.0) continue;
            cout << "time: " << time << "\r";


        }

        wK1.lap();
        tau *= std::pow(delta/error, 0.2) * 0.8;
        if (tau == 0 ) {
            problem.exportData((time), compData, 1);
            throw std::runtime_error("zero tau at time: " + to_string(time));
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

enum BOUNDARY_SETUP{
    BOUNDARY_SETUP_BOILER2D,
    BOUNDARY_SETUP_BOILER3D,
    BOUNDARY_SETUP_CUBE,
    BOUNDARY_SETUP_STACK
};

template <BOUNDARY_SETUP>
struct BoundaryCondition {
    static constexpr unsigned int problemDim(){return 3;}
    static constexpr unsigned int g_axes(){return 2;}
    static std::string meshName() {return "";}
    template<unsigned int MeshDim>
    static double
    inFlowModulation(const Vertex<MeshDim, double>& x)
    {

        //double inFlowModulation = x[0];
        double inFlowModulation = sqrt(pow(x[0],2) + pow(x[2],2));

        //inFlowModulation = (- inFlowModulation * inFlowModulation + inFlowModulation - 0.0099) * 4.164931279;
        //inFlowModulation = -(inFlowModulation -1.9) * (inFlowModulation - 4.7) * 0.5102;//(-inFlowModulation * inFlowModulation + 6.6 * inFlowModulation - 8.93) * 0.5102;
        //inFlowModulation = - (inFlowModulation + 1.25) * (inFlowModulation - 1.25) * (1/1.5625);
        inFlowModulation = - (inFlowModulation + 0.05) * (inFlowModulation - 0.05) * 400;
        return inFlowModulation;
    }

    template<typename Cell>
    static Type
    TypeOfCell(const Cell &cell) {
        if (cell.getIndex() >= BOUNDARY_INDEX(size_t)){

            if (
                //    cell.getCenter()[1] <= 1e-5
                //    sqrt(pow(cell.getCenter()[0],2) + pow(cell.getCenter()[1],2)) < 1.25 &&
                //    cell.getCenter()[2] <= -1.249
                    sqrt(pow(cell.getCenter()[0],2) + pow(cell.getCenter()[2],2)) < 0.05 &&
                    cell.getCenter()[1] < 1e-5
                ) {
                return Type::INFLOW;
            }

            if (
                //    cell.getCenter()[1] >= 34.349
                    sqrt(pow(cell.getCenter()[0],2) + pow(cell.getCenter()[2],2)) < 0.075 &&
                    cell.getCenter()[1] > 2.099
                //    sqrt(pow(cell.getCenter()[0],2) + pow(cell.getCenter()[1],2)) < 1.0 &&
                //    cell.getCenter()[2] > 0
                //cell.GetCenter()[0] > 2 && (cell.GetCenter()[1] < 32.0 && cell.GetCenter()[1] > 30)
                //cell.GetCenter()[0] > 8 && cell.GetCenter()[1] >= 33.39
                //cell.getCenter()[1] >= 1.999
               ) {
                return Type::OUTFLOW;
            }

            return Type::WALL;

        } else {

            return Type(cell.getFlag());

        }


        throw(std::runtime_error ("cell type not recognized " + std::to_string(cell.getIndex())));
    }
};


template <> constexpr unsigned int BoundaryCondition<BOUNDARY_SETUP::BOUNDARY_SETUP_BOILER3D>::problemDim() {return 3;}
template <> constexpr unsigned int BoundaryCondition<BOUNDARY_SETUP::BOUNDARY_SETUP_BOILER3D>::g_axes() {return 1;}
template <> std::string BoundaryCondition<BOUNDARY_SETUP::BOUNDARY_SETUP_BOILER3D>::meshName() {return "boiler.vtk";}
template <>
template<unsigned int MeshDim>
double
BoundaryCondition<BOUNDARY_SETUP::BOUNDARY_SETUP_BOILER3D>::
inFlowModulation(const Vertex<MeshDim, double>& x)
{
    double inFlowModulation = sqrt(pow(x[0],2) + pow(x[2],2));

    inFlowModulation = - (inFlowModulation + 0.05) * (inFlowModulation - 0.05) * 400;
    return inFlowModulation;
}

template <>
template<typename Cell>
Type
BoundaryCondition<BOUNDARY_SETUP::BOUNDARY_SETUP_BOILER3D>::
TypeOfCell(const Cell &cell) {
    if (cell.getIndex() >= BOUNDARY_INDEX(size_t)){

        if (
                sqrt(pow(cell.getCenter()[0],2) + pow(cell.getCenter()[2],2)) < 0.05 &&
                cell.getCenter()[1] < 1e-5
            ) {
            return Type::INFLOW;
        }

        if (
                sqrt(pow(cell.getCenter()[0],2) + pow(cell.getCenter()[2],2)) < 0.075 &&
                cell.getCenter()[1] > 2.099
           ) {
            return Type::OUTFLOW;
        }

        return Type::WALL;

    } else {

        return Type(cell.getFlag());

    }
    throw(std::runtime_error ("cell type not recognized " + std::to_string(cell.getIndex())));
}


template <> constexpr unsigned int BoundaryCondition<BOUNDARY_SETUP::BOUNDARY_SETUP_CUBE>::problemDim() {return 3;}
template <> constexpr unsigned int BoundaryCondition<BOUNDARY_SETUP::BOUNDARY_SETUP_CUBE>::g_axes() {return 2;}
template <> std::string BoundaryCondition<BOUNDARY_SETUP::BOUNDARY_SETUP_CUBE>::meshName() {return "";}
template <>
template<unsigned int MeshDim>
double
BoundaryCondition<BOUNDARY_SETUP::BOUNDARY_SETUP_CUBE>::
inFlowModulation(const Vertex<MeshDim, double>& x)
{

    double inFlowModulation = sqrt(pow(x[0],2) + pow(x[1],2));

    inFlowModulation = - (inFlowModulation + 0.3) * (inFlowModulation - 0.3) * 1.0/(0.3 * 0.3);
    return inFlowModulation;
}


template <>
template<typename Cell>
Type
BoundaryCondition<BOUNDARY_SETUP::BOUNDARY_SETUP_CUBE>::
TypeOfCell(const Cell &cell) {
    if (cell.getIndex() >= BOUNDARY_INDEX(size_t)){

        if (
                sqrt(pow(cell.getCenter()[0],2) + pow(cell.getCenter()[1],2)) < 0.3 &&
                cell.getCenter()[2] <= -1.249
            ) {
            return Type::INFLOW;
        }

        if (
                sqrt(pow(cell.getCenter()[0],2) + pow(cell.getCenter()[1],2)) < 0.3 &&
                cell.getCenter()[2] > 0
           ) {
            return Type::OUTFLOW;
        }

        return Type::WALL;

    } else {

        return Type(cell.getFlag());

    }
    throw(std::runtime_error ("cell type not recognized " + std::to_string(cell.getIndex())));
}


template <> constexpr unsigned int BoundaryCondition<BOUNDARY_SETUP::BOUNDARY_SETUP_STACK>::problemDim() {return 3;}
template <> constexpr unsigned int BoundaryCondition<BOUNDARY_SETUP::BOUNDARY_SETUP_STACK>::g_axes() {return 2;}
template <> std::string BoundaryCondition<BOUNDARY_SETUP::BOUNDARY_SETUP_STACK>::meshName() {return "";}
template <>
template<unsigned int MeshDim>
double
BoundaryCondition<BOUNDARY_SETUP::BOUNDARY_SETUP_STACK>::
inFlowModulation(const Vertex<MeshDim, double>& x)
{

    double inFlowModulation = sqrt(pow(x[0],2) + pow(x[1],2));

    inFlowModulation = - (inFlowModulation + 1.0) * (inFlowModulation - 1.0);
    return inFlowModulation;
}


template <>
template<typename Cell>
Type
BoundaryCondition<BOUNDARY_SETUP::BOUNDARY_SETUP_STACK>::
TypeOfCell(const Cell &cell) {
    if (cell.getIndex() >= BOUNDARY_INDEX(size_t)){

        if (
                sqrt(pow(cell.getCenter()[0],2) + pow(cell.getCenter()[1],2)) < 1.0 &&
                cell.getCenter()[2] <= -1.249
            ) {
            return Type::INFLOW;
        }

        if (
                sqrt(pow(cell.getCenter()[0],2) + pow(cell.getCenter()[1],2)) < 1.0 &&
                cell.getCenter()[2] > 0
           ) {
            return Type::OUTFLOW;
        }

        return Type::WALL;

    } else {

        return Type(cell.getFlag());

    }
    throw(std::runtime_error ("cell type not recognized " + std::to_string(cell.getIndex())));
}



template <> constexpr unsigned int BoundaryCondition<BOUNDARY_SETUP::BOUNDARY_SETUP_BOILER2D>::problemDim() {return 2;}
template <> constexpr unsigned int BoundaryCondition<BOUNDARY_SETUP::BOUNDARY_SETUP_BOILER2D>::g_axes() {return 1;}
template <> std::string BoundaryCondition<BOUNDARY_SETUP::BOUNDARY_SETUP_BOILER2D>::meshName() {return "Boiler2D.vtk";}
template <>
template<unsigned int MeshDim>
double
BoundaryCondition<BOUNDARY_SETUP::BOUNDARY_SETUP_BOILER2D>::
inFlowModulation(const Vertex<MeshDim, double>& x)
{
    //return 1;
    double inFlowModulation = x[0];
    return -(inFlowModulation -1.9) * (inFlowModulation - 4.7) * 0.5102;
}

template <>
template<typename Cell>
Type
BoundaryCondition<BOUNDARY_SETUP::BOUNDARY_SETUP_BOILER2D>::
TypeOfCell(const Cell &cell) {
    if (cell.getIndex() >= BOUNDARY_INDEX(size_t)){
/*
        if (
                cell.getCenter()[1] <= 1e-5
            ) {
            return Type::INFLOW;
        }

        if (
                cell.getCenter()[1] >= 34.349
           ) {
            return Type::OUTFLOW;
        }
*/
        return Type::WALL;

    } else {

        return Type(cell.getFlag());

    }
    throw(std::runtime_error ("cell type not recognized " + std::to_string(cell.getIndex())));
}

void MultiphaseFlowCalculation(string name) {
    constexpr unsigned int Reserve = 12;
    using BC = BoundaryCondition<BOUNDARY_SETUP::BOUNDARY_SETUP_BOILER2D>;
    string meshName = "../" + BC::meshName();
    if (!name.empty()) {
        meshName = name;
    }
    using MPFType = MultiphaseFlow<BC::problemDim(), BC, Reserve>;

    MPFType mpf;

    // setup constants
    MPFType::ResultType::R_spec = 287;
    MPFType::ResultType::T = 300;
    MPFType::ResultType::rho_s = 1700;


    mpf.artifitialDisspationGas = 0.1;
    mpf.artifitialDisspationSolid = 0.1;
    mpf.R_spec = MPFType::ResultType::R_spec;
    mpf.myu = 1e-5;
    mpf.rho_s = MPFType::ResultType::rho_s;
    mpf.myu_s = 0.5;//1.5;
    mpf.d_s = 0.00078;
    mpf.phi_s = 1;

    mpf.T = MPFType::ResultType::T;

    mpf.setupMeshData(meshName);



    mpf.outFlow.eps_s = 0;
    mpf.outFlow.setPressure(1e5);

    mpf.inFlow_eps_g = 1;
    mpf.inFlow_eps_s = 0;
    mpf.inFlow_u_g = {};
    mpf.inFlow_u_g[BC::g_axes()] = 25;
    mpf.inFlow_u_s = {};


    FlowData<BC::problemDim()> ini;
    ini.eps_s = 0;
    ini.setPressure(1e5);
    //ini.setRho_g(0.3564);
    ini.p_g = {};
    ini.p_s = {};

    JSONLogger jsl("setup.json");
    jsl.writeVar("MultiphaseFlow setup", mpf);
    jsl.writeVar("mesh name", meshName);
    jsl.writeVar("MPFType::ResultType::R_spec", MPFType::ResultType::R_spec);
    jsl.writeVar("MPFType::ResultType::T", MPFType::ResultType::T);
    jsl.writeVar("MPFType::ResultType::rho_s", MPFType::ResultType::rho_s);


    MeshDataContainer<FlowData<BC::problemDim()>, BC::problemDim()> compData(mpf.mesh, ini);


    for (auto& cell : mpf.mesh.getCells()){
        if(
                //cell.getCenter()[1] > 1.5 && cell.getCenter()[1] < 1.6
                cell.getCenter()[1] > 2 && cell.getCenter()[1] < 14
                //&& cell.getCenter()[0] * cell.getCenter()[0] + cell.getCenter()[2] * cell.getCenter()[2] < 0.03*0.03
          ){
            compData.at(cell).eps_s = 0.2;
            compData.at(cell).setPressure(1e5);
        }
    }

/*
    for (auto& cell : mpf.mesh.getCells()){
        if(
                cell.getCenter()[1] < 1 && cell.getCenter()[1] > 0.05
             //  cell.getCenter()[2] > -1.0 && cell.getCenter()[2] < 2.0
//           cell.getCenter()[1] > -0.5 && cell.getCenter()[1] < 0.5 &&
//           cell.getCenter()[0] > -0.5 && cell.getCenter()[0] < 0.5
        ) {
            compData.at(cell).eps_s = 0.5;
            compData[cell].setPressure(1e5);
        }
    }
*/

    mpf.exportData(0.0, compData);
    double exportStep = 1e-1;
    for (double t = 0; t < 100 * exportStep; t += exportStep){


        RKMSolver(mpf, compData, 1e-7, t, t + exportStep, 1e-4);

        mpf.exportData((t + exportStep), compData, 1.0/exportStep);


        DBGVAR(t + exportStep, wK1.getResult());
    }




}

void atEnd(){
    DBGVAR(wK1.getResult());
}

int main(int argc, const char** argv)
{
    atexit(atEnd);
    MultiphaseFlowCalculation(argc > 1 ? argv[1] : "");

}

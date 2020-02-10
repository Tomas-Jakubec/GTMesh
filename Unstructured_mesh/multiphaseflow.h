#ifndef MULTIPHASEFLOW_H
#define MULTIPHASEFLOW_H


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

#include <list>
#include <valarray>
#include <fstream>

/**
 * @brief The FlowData struct
 * flow computation data
 */
struct FlowData {

// gasseous part
    double p;

    double rho_g;

    /**
     * @brief eps_g
     * volume fraction of gaseous part of the flow
     * holds \epsilon_{s} + \epsilon_{g} = 1
     */
    double eps_g;

    /**
     * @brief u_g
     *
     * velocity of gaseous phase
     */
    Vector<2, double> getVelocityGas(){
        return p_g / rho_g;
    }


    Vector<2, double> setVelocityGas(const Vector<2, double>& u_g){
        return p_g = u_g * rho_g;
    }

    /**
     * @brief p_g
     *
     * momentum of gaseous phase
     */
    Vector<2,double> p_g;

    double fluxRho_g; //flux of "mass" over cell boundary

    Vector<2,double> fluxP_g; //flux of momentum over cell boundary


// solid part
    /**
     * @brief u_s
     *
     * velocity of solid component
     */
    Vector<2, double> getVelocitySolid(){
        return p_s / rho_s;
    }


    Vector<2, double> setVelocitySolid(const Vector<2, double>& u_s){
        return p_s = u_s * rho_s;
    }

    static double rho_s;

    /**
     * @brief p_s
     *
     * momentum of solid phase
     */
    Vector<2,double> p_s;


    /**
     * @brief eps_g
     * volume fraction of solid part of the flow
     * holds \epsilon_{s} + \epsilon_{g} = 1
     */
    double eps_s;


    double fluxRho_s;
    Vector<2,double> fluxP_s;

};



MAKE_ATTRIBUTE_TRAIT_ARITHMETIC(FlowData, p, eps_g, rho_g, eps_s, p_g, p_s)


MAKE_CUSTOM_ATTRIBUTE_TRAIT_IO(FlowData,
                           "pressure", &FlowData::p,
                           "eps_g", &FlowData::eps_g,
                           "rho_g", &FlowData::rho_g,
                           "eps_s", &FlowData::eps_s
                           //"velocity_gas", std::make_pair(&FlowData::getVelocityGas, &FlowData::setVelocityGas),
                           //"velocity_solid", std::make_pair(&FlowData::getVelocitySolid, &FlowData::setVelocitySolid)
                            )




struct EdgeData {
    double LengthOverDist;
    double Length;
    Vector<2,double> n;

    // koeficients of konvex combination of cells
    // holds that LeftCellKoef + RightCellKoef = 1
    double LeftCellKoef;
    double RightCellKoef;


    /**
     * Next variables are supposed to store fluxes over the edge
     * in order to simply parallelize the computation
     */

    /**
     * @brief fluxP_g
     *
     * flux of momentum of gaseous phase over the edge
     */
    Vector<2,double> fluxP_g;


    /**
     * @brief fluxP_s
     *
     * flux of momentum of solid phase over the edge
     */
    Vector<2,double> fluxP_s;


    /**
     * @brief fluxRho_s
     *
     * flux of mass of solid phase over the edge
     */
    double fluxRho_s;

    /**
     * @brief fluxRho_g
     *
     * flux of mass of gaseous phase over the edge
     */
    double fluxRho_g; //flux of "mass" over cell boundary
};

// Enumerator Type for expessing other cell and point position
enum Type{
    INNER = 1,
    WALL = 2,
    INFLOW = 4,
    OUTFLOW = 8,
};

// Data mapped to point
struct PointData
{
    // velocity of gas in the particular point
    Vector<2,double> u_g;
    /**
     * @brief u_s
     * averaged velocity of solid component in vertices
     */
    Vector<2,double> u_s;

    // type of the point expessing its position
    Type PointType;

    // inverted number of neighboring cells
    double cellKoef;
};


struct CellData{
    double invVolume;
};

class MultiphaseFlow {
public:
    MultiphaseFlow() = default;
    ~MultiphaseFlow() = default;

public: // make some types public
    using MeshType = UnstructuredMesh<2,size_t, double>;
    using ResultType = FlowData;
public:
    // Structures containing state data

    MeshType mesh;

    MeshDataContainer<std::tuple<CellData,EdgeData,PointData>, 2,1,0> meshData;
private:

    VTKMeshWriter<2> writer;

public:

private:


// these atributes are public to be set by user
public:

    double myu, myu_s;
    double R_spec;
    /**
     * @brief T
     * gas temperature (remains constant)
     */
    double T;
    /**
     * @brief artifitialDisspation
     * Coefitient of artifitial diffusion added to calculation of convective flux at all inner edges in order to guarantee numerial stability.
     */
    double artifitialDisspation;
    /**
     * @brief d_s
     * Average diameter of the solid particle.
     */
    double d_s;
    /**
     * @brief phi_s
     * Sphericity of average solid particle.
     */
    double phi_s;

    /**
     * @brief rho_s
     * density od solid component (remains constant)
     * if the combustion is added, then the density of the combustible part may change in time
     */
    double rho_s;
    FlowData boundaryWall;
    FlowData inFlow;
    FlowData outFlow;


private:
    /**
      * Functions computing special factors in formula
      */

    /**
     * @brief G
     * @param eps_g
     * @return
     */
//    inline double G(double eps_g);

    /**
     * @brief Beta_s
     * @param fd
     * @return
     */
//    inline double Beta_s(const FlowData& fd);

    /**
     * @brief C_d
     * @param fd
     * @return
     */
//    inline double C_d(const FlowData& fd);

    /**
     * @brief Re_s
     * @param fd
     * @return
     */
//    inline double Re_s(const FlowData& fd);

    /**
     * @brief ComputeFluxGas_inner
     * @param leftData
     * @param rightData
     * @param edgeData
     * @param fcData
     */
//    inline void ComputeFluxGas_inner(FlowData& leftData, FlowData& rightData, EdgeData& edgeData, FluxComputeData &fcData);
//
//    inline void ComputeFluxGas_inflow(FlowData& innerCellData, const Cell& innerCell, EdgeData &edgeData, FluxComputeData &fcData);
//
//    inline void ComputeFluxGas_outflow(FlowData& innerCellData, EdgeData& edgeData, FluxComputeData &fcData);
//
//    inline void ComputeFluxGas_wall(FlowData& innerCellData, EdgeData &edgeData, FluxComputeData &fcData);

    /**
     * @brief ComputeFluxSolid_inner
     * @param leftData
     * @param rightData
     * @param edgeData
     * @param fcData
     */
//    inline void ComputeFluxSolid_inner(FlowData& leftData, FlowData& rightData,  EdgeData &edgeData, FluxComputeData& fcData);
//
//    inline void ComputeFluxSolid_inflow(FlowData& innerCellData, const Cell& innerCell,  EdgeData& edgeData, FluxComputeData &fcData);
//
//    inline void ComputeFluxSolid_outflow(FlowData& innerCellData,  EdgeData& edgeData, FluxComputeData &fcData);
//
//    inline void ComputeFluxSolid_wall(FlowData& innerCellData, EdgeData &edgeData, FluxComputeData &fcData);

    /**
     * @brief reg
     * regularizes number to aviod zero division.
     * Simply adds small non-zero constant (1e-7)
     * @param x
     * number to be regularized
     * @return
     * greatered x
     */
//    inline double reg(double x);

private:

    // compute flux over edges for all cells
    //class Flux;



    double FlowModulation(double x);


public:


//    void LoadInitialData(std::istream& ist);


public:

//    void ComputeFlux(FVMethod::FluxComputeData& fcData);
//
//    void ComputeSource(FVMethod::CellComputeData& ccData);
//
    void setupMeshData(const std::string& fileName);
//
//
//    void Compute(size_t NumberOfSteps);
//
//    void ComputeStep();
//
//    void DataWrite(std::ostream& ost);
//
//    void SetTimeStep(double TS);
//
//
//    // Computes velocities in the points on the start of the program
//    void InitializePointData();
//
//    void ActualizePointData();
//
//    // Calculates velocities during the computation step for each cell
//    void CalculateVertexData(Mesh::MeshCell &cell, const NumData<FlowData> &cellData);
//
//    Type TypeOfCell(const MeshType::Cell& cell);



    template<typename dataType>
    void exportData(double time,
                    MeshDataContainer<dataType, 3>& compData) {


        std::ofstream ofile(std::string("MultiphaseFlow") + "_" + std::to_string(time) + ".vtk");
        writer.writeHeader(ofile, std::string("MPF") + std::to_string(time));
        writer.writeToStream(ofile, mesh, MeshDataContainer<MeshNativeType<2>::ElementType, 2>(mesh, MeshNativeType<2>::POLYGON));

        VTKMeshDataWriter<2> dataWriter;
        dataWriter.writeToStream(ofile, compData, writer);

        ofile.close();
        DBGMSG("Data eported");

    }
};


#ifdef ___

/*
 * Definition of inlined functions
 *
*/
inline double MultiphaseFlow::G(double eps_g)
{
    return (std::pow(10, -8.76 * eps_g + 5.43));
}






inline double MultiphaseFlow::Beta_s(const FlowData &fd)
{
    if(fd.eps_s > 0.2) {
        double denominator = 1 / (fd.eps_g * d_s * phi_s);
        return (
                    150 * (std::pow(fd.eps_s, 2.0) * myu) * std::pow(denominator, 2.0) +
                    1.75 * (fd.u_g - fd.u_s).Norm() * fd.rho_g * fd.eps_s * denominator
                );
    } else {
        return (
                    (4.0/3.0) * C_d(fd) * (fd.u_g - fd.u_s).Norm() * fd.rho_g * fd.eps_s / (d_s * phi_s)
               );
    }
}






inline double MultiphaseFlow::C_d(const FlowData &fd)
{
    double re_s = reg(Re_s(fd));

    if (re_s <= 1000) {
        return (24/re_s)*(1.0 + 0.15 * std::pow(re_s, 0.687));
    } else {
        return 0.44;
    }
}






inline double MultiphaseFlow::Re_s(const FlowData &fd)
{
    // multiplying by inverted value of myu may make code faster
    return(
            (fd.u_g - fd.u_s).Norm() * d_s * fd.rho_g * fd.eps_g /  myu
          );
}




inline void MultiphaseFlow::ComputeFluxGas_inner(FlowData &leftData, FlowData &rightData, EdgeData &edgeData, FluxComputeData& fcData)
{
    /**
     * @brief eps_g_e
     * volume fraction aproximated at the edge
     */
    double eps_g_e = (leftData.eps_g * edgeData.LeftCellKoef + rightData.eps_g * edgeData.RightCellKoef);


    // computing derivatives of velocity
    Vector<2,double> du_g_dn = (rightData.u_g - leftData.u_g) * (edgeData.LengthOverDist);

    Vector<2,double> du_g_dt = (PointData.at(fcData.GetCurrEdge().GetPointBIndex()).u_g -
                      PointData.at(fcData.GetCurrEdge().GetPointAIndex()).u_g);

    //transform diferences to standard base
    Vector<2,double> du_g_dx = du_g_dn * (edgeData.n.x) - (du_g_dt * edgeData.n.y);
    Vector<2,double> du_g_dy = du_g_dt * (edgeData.n.x) + (du_g_dn * edgeData.n.y);


    double productOf_u_And_n = (((leftData.u_g * edgeData.LeftCellKoef) +
                                 (rightData.u_g * edgeData.RightCellKoef))
                                * edgeData.n);

    // flux of density
    double delta_rho = -(leftData.eps_g * leftData.rho_g * edgeData.LeftCellKoef +
                         rightData.eps_g * rightData.rho_g * edgeData.RightCellKoef) *
                       productOf_u_And_n * edgeData.Length +
                       (rightData.rho_g * rightData.eps_g - leftData.rho_g * leftData.eps_g) * edgeData.LengthOverDist * artifitialDisspation;


    // computing the flux of momentum
    Vector<2,double> fluxP_g = (-productOf_u_And_n) *
                  (leftData.p_g * edgeData.LeftCellKoef + rightData.p_g * edgeData.RightCellKoef);

    // adding the element of pressure gradient
    fluxP_g -= (leftData.p * edgeData.LeftCellKoef +
                 rightData.p * edgeData.RightCellKoef) *
               edgeData.n;


    fluxP_g *= edgeData.Length;



    // add artifitial dissipation
    fluxP_g += (edgeData.LengthOverDist * artifitialDisspation * ((rightData.p_g) - (leftData.p_g)));


    //viscose_x_x = (4.0/3) * myu * du_x_dx * edgeData.n.x * (leftData.eps_g * edgeData.LeftCellKoef + rightData.eps_g * edgeData.RightCellKoef);
    //viscose_y_y = (4.0/3) * myu * du_y_dy * edgeData.n.y * (leftData.eps_g * edgeData.LeftCellKoef + rightData.eps_g * edgeData.RightCellKoef);
    Vector<2,double> viscose_diag = ((4.0/3) * myu * eps_g_e) * Vector<2,double>(du_g_dx.x * edgeData.n.x, du_g_dy.y * edgeData.n.y);
    // Diffusion of the velocity
    // firstly i will consider diffusion only in the direction

    //viscose_x_y = myu * (du_y_dx + du_x_dy) * edgeDat (leftData.eps_g * edgeData.LeftCellKoef + rightData.eps_g * edgeData.RightCellKoef);
    //viscose_y_x = myu * (du_y_dx + du_x_dy) * edgeData.n.x * (leftData.eps_g * edgeData.LeftCellKoef + rightData.eps_g * edgeData.RightCellKoef);
    Vector<2,double> viscose_side_diag = (myu * eps_g_e *(du_g_dx.y + du_g_dy.x)) * Vector<2,double>(edgeData.n.y,edgeData.n.x);

    // adding of fluxes to cells
    edgeData.fluxP_g = fluxP_g + viscose_diag + viscose_side_diag;

    edgeData.fluxRho_g = delta_rho;



}

inline void MultiphaseFlow::ComputeFluxGas_inflow(FlowData &innerCellData, const Cell& innerCell, EdgeData &edgeData, FVMethod::FluxComputeData &fcData)
{

    Vector<2,double> modulatedU = inFlow.getVelocityGas();

    modulatedU.y *= FlowModulation(innerCell.GetCenter().GetX());

    double productOf_u_And_n = (modulatedU * edgeData.n);

    /*
    ** Preparing prtial derivatives wrt
    ** edge orientation
    */

    // computing derivatives of velocity
    Vector<2,double> du_g_dn(0.0, 0.0); // boundary cond
    Vector<2,double> du_g_dt = (PointData.at(fcData.GetCurrEdge().GetPointBIndex()).u_g -
                      PointData.at(fcData.GetCurrEdge().GetPointAIndex()).u_g);

    // transform derivatives to standard base
    Vector<2,double> du_g_dx = du_g_dn * (edgeData.n.x) - (du_g_dt * edgeData.n.y);
    Vector<2,double> du_g_dy = du_g_dt * (edgeData.n.x) + (du_g_dn * edgeData.n.y);

    // flux of density
    double delta_rho = -innerCellData.rho_g * inFlow.eps_g *
                        productOf_u_And_n * edgeData.Length;




    // computing the flux of momentum
    Vector<2,double> flux = - (innerCellData.rho_g) *
                    (productOf_u_And_n) * inFlow.eps_g *
                    (modulatedU);


    // adding the element of pressure gradient
    flux -= (innerCellData.p) * edgeData.n;


    flux *= edgeData.Length;



    // Diffusion of the velocity
    //viscose_x_x = (4.0/3) * myu * du_x_dx  * edgeData.n.x;
    //viscose_y_y = (4.0/3) * myu * du_y_dy * edgeData.n.y;
    Vector<2,double> viscose_diag = ((4.0/3.0) * myu * inFlow.eps_g) * Vector<2,double>(du_g_dx.x * edgeData.n.x, du_g_dy.y * edgeData.n.y);

    // Diffusion of the velocity

    //viscose_y_x = myu * (du_x_dy + du_y_dx) * edgeData.n.x;
    //viscose_x_y = myu * (du_x_dy + du_y_dx) * edgeData.n.y;
    Vector<2,double> viscose_side_diag = (myu * inFlow.eps_g * (du_g_dx.y + du_g_dy.x)) * Vector<2,double>(edgeData.n.y,edgeData.n.x);

    edgeData.fluxP_g = flux + viscose_side_diag + viscose_diag;

    edgeData.fluxRho_g = delta_rho;

}

inline void MultiphaseFlow::ComputeFluxGas_outflow(FlowData &innerCellData, EdgeData &edgeData, FVMethod::FluxComputeData &fcData)
{
    double productOf_u_And_n = ((innerCellData.u_g) * edgeData.n);
    double eps_g_e = productOf_u_And_n > 0 ? innerCellData.eps_g : outFlow.eps_g;

    /*
    ** Preparing prtial derivatives wrt
    ** edge orientation
    */

    // computing derivatives of velocity
    Vector<2,double> du_g_dn(0.0, 0.0); // boundary cond
    Vector<2,double> du_g_dt = (PointData.at(fcData.GetCurrEdge().GetPointBIndex()).u_g -
                      PointData.at(fcData.GetCurrEdge().GetPointAIndex()).u_g);

    // transform derivatives to standard base
    Vector<2,double> du_g_dx = du_g_dn * (edgeData.n.x) - (du_g_dt * edgeData.n.y);
    Vector<2,double> du_g_dy = du_g_dt * (edgeData.n.x) + (du_g_dn * edgeData.n.y);



    // flux of density
    double delta_rho = -outFlow.rho_g * eps_g_e *
                        productOf_u_And_n * edgeData.Length;


    // computing the flux of momentum
    Vector<2,double> flux = - (outFlow.rho_g *
                    productOf_u_And_n * eps_g_e) *
                    (innerCellData.u_g);


    // adding the element of pressure gradient
    flux -= (outFlow.p) * edgeData.n;

    flux *= edgeData.Length;



    // Diffusion of the velocity
    //viscose_x_x = (4.0/3) * myu * du_x_dx  * edgeData.n.x;
    //viscose_y_y = (4.0/3) * myu * du_y_dy * edgeData.n.y;
    Vector<2,double> viscose_diag = ((4.0/3.0) * myu * eps_g_e) * Vector<2,double>(du_g_dx.x * edgeData.n.x, du_g_dy.y * edgeData.n.y);

    // Diffusion of the velocity

    //viscose_y_x = myu * (du_x_dy + du_y_dx) * edgeData.n.x;
    //viscose_x_y = myu * (du_x_dy + du_y_dx) * edgeData.n.y;
    Vector<2,double> viscose_side_diag = (myu * eps_g_e * (du_g_dx.y + du_g_dy.x)) * Vector<2,double>(edgeData.n.y,edgeData.n.x);


    edgeData.fluxP_g = flux + viscose_side_diag + viscose_diag;

    edgeData.fluxRho_g = delta_rho;

}

inline void MultiphaseFlow::ComputeFluxGas_wall(FlowData &innerCellData, EdgeData &edgeData, FVMethod::FluxComputeData &fcData)
{
    // computing derivatives of velocity
    Vector<2,double> du_g_dn = edgeData.LengthOverDist * (fcData.GetCellLeft().GetCellFlag() == INNER ? -1 * innerCellData.u_g : innerCellData.u_g);  // boundary cond

    Vector<2,double> du_g_dt = Vector<2,double>(0.0, 0.0);

    // transform derivatives to standard base
    Vector<2,double> du_g_dx = du_g_dn * (edgeData.n.x) - (du_g_dt * edgeData.n.y);
    Vector<2,double> du_g_dy = du_g_dt * (edgeData.n.x) + (du_g_dn * edgeData.n.y);




    // the flux of momentum
    // is reduced only to pressure gradient
    Vector<2,double> flux = -(innerCellData.p) * edgeData.n;


    flux *= edgeData.Length;

    // Diffusion of the velocity
    Vector<2,double> viscose_diag = ((4.0/3.0) * myu * innerCellData.eps_g) * Vector<2,double>(du_g_dx.x * edgeData.n.x, du_g_dy.y * edgeData.n.y);
    Vector<2,double> viscose_side_diag = (myu * innerCellData.eps_g * (du_g_dx.y + du_g_dy.x)) * Vector<2,double>(edgeData.n.y,edgeData.n.x);



    edgeData.fluxP_g = flux + viscose_side_diag + viscose_diag;

}













inline void MultiphaseFlow::ComputeFluxSolid_inner(FlowData &leftData, FlowData &rightData, EdgeData &edgeData, FVMethod::FluxComputeData &fcData)
{
    /**
     * @brief eps_s_e
     * volume fraction aproximated at the edge
     */
    double eps_s_e = (leftData.eps_s * edgeData.LeftCellKoef + rightData.eps_s * edgeData.RightCellKoef);


    Vector<2,double> du_s_dn = (rightData.u_s - leftData.u_s) * edgeData.LengthOverDist;

    Vector<2,double> du_s_dt = (PointData.at(fcData.GetCurrEdge().GetPointBIndex()).u_s -
                      PointData.at(fcData.GetCurrEdge().GetPointAIndex()).u_s);

    //transform diferences to standard base
    Vector<2,double> du_s_dx = du_s_dn * (edgeData.n.x) - (du_s_dt * edgeData.n.y);

    Vector<2,double> du_s_dy = du_s_dt * (edgeData.n.x) + (du_s_dn * edgeData.n.y);

    double product_of_u_s_and_n = (((leftData.u_s * edgeData.LeftCellKoef) +
                                  (rightData.u_s * edgeData.RightCellKoef)) *
                                  edgeData.n);


    // flux of mass
    double fluxRho_s = (-rho_s * eps_s_e * edgeData.Length * product_of_u_s_and_n  +
                       (rightData.eps_s - leftData.eps_s) * rho_s * edgeData.LengthOverDist * artifitialDisspation) ;


    // computing the flux of momentum
    Vector<2,double> fluxP_s = -product_of_u_s_and_n *
                     (leftData.p_s * edgeData.LeftCellKoef + rightData.p_s * edgeData.RightCellKoef);


    // adding the element of pressure gradient
    fluxP_s -= G(leftData.eps_g * edgeData.LeftCellKoef + rightData.eps_g * edgeData.RightCellKoef) * eps_s_e *
              edgeData.n;

    fluxP_s *= edgeData.Length;

    // add artifitial dissipation
    // TODO this is the critical point it is necessary to look into this problem
    // maybe the artifitial disipation must be raise
    fluxP_s += (edgeData.LengthOverDist * artifitialDisspation ) * (rightData.p_s - leftData.p_s);




    //viscose_x_x = (4.0/3) * myu * du_x_dx * edgeData.n.x * (leftData.eps_g * edgeData.LeftCellKoef + rightData.eps_g * edgeData.RightCellKoef);
    //viscose_y_y = (4.0/3) * myu * du_y_dy * edgeData.n.y * (leftData.eps_g * edgeData.LeftCellKoef + rightData.eps_g * edgeData.RightCellKoef);
    //Vector<2,double> viscose_diag = ((4.0/3) * myu_s * eps_s_e) * Vector<2,double>(du_s_dx.x * edgeData.n.x, du_s_dy.y * edgeData.n.y);
    fluxP_s += ((4.0/3) * myu_s * eps_s_e) * Vector<2,double>(du_s_dx.x * edgeData.n.x, du_s_dy.y * edgeData.n.y);

    //Vector<2,double> viscose_diag = (4.0/3) * myu * du_s_dy * edgeData.n.y * (leftData.eps_g * edgeData.LeftCellKoef + rightData.eps_g * edgeData.RightCellKoef);


    //viscose_x_y = myu * (du_y_dx + du_x_dy) * edgeData.n.y * (leftData.eps_g * edgeData.LeftCellKoef + rightData.eps_g * edgeData.RightCellKoef);
    //viscose_y_x = myu * (du_y_dx + du_x_dy) * edgeData.n.x * (leftData.eps_g * edgeData.LeftCellKoef + rightData.eps_g * edgeData.RightCellKoef);
    //Vector<2,double> viscose_side_diag = (myu_s * eps_s_e *(du_s_dx.y + du_s_dy.x)) * Vector<2,double>(edgeData.n.y,edgeData.n.x);
    fluxP_s += (myu_s * eps_s_e *(du_s_dx.y + du_s_dy.x)) * Vector<2,double>(edgeData.n.y,edgeData.n.x);

    // storing result at the edge
    edgeData.fluxP_s = fluxP_s; // + viscose_diag + viscose_side_diag;

    edgeData.fluxRho_s = fluxRho_s;
}

inline void MultiphaseFlow::ComputeFluxSolid_inflow(FlowData &innerCellData, const Cell &innerCell, EdgeData &edgeData, FVMethod::FluxComputeData &fcData)
{
    (void)innerCellData;

    Vector<2,double> modulatedU = inFlow.u_s;

    modulatedU.y *= FlowModulation(innerCell.GetCenter().GetX());

    double productOf_u_And_n = (modulatedU * edgeData.n);

    /*
    ** Preparing prtial derivatives wrt
    ** edge orientation
    */

    // computing derivatives of velocity
    Vector<2,double> du_s_dn(0.0, 0.0); // boundary cond
    Vector<2,double> du_s_dt = (PointData.at(fcData.GetCurrEdge().GetPointBIndex()).u_s -
                      PointData.at(fcData.GetCurrEdge().GetPointAIndex()).u_s);

    // transform derivatives to standard base
    Vector<2,double> du_s_dx = du_s_dn * (edgeData.n.x) - (du_s_dt * edgeData.n.y);
    Vector<2,double> du_s_dy = du_s_dt * (edgeData.n.x) + (du_s_dn * edgeData.n.y);

    // flux of density
    double fluxRho_s = -rho_s * inFlow.eps_s *
                        productOf_u_And_n * edgeData.Length;




    // computing the flux of momentum
    Vector<2,double> fluxP_s = - (rho_s) *
                    (productOf_u_And_n) * inFlow.eps_s *
                    (modulatedU);

    // adding the element of pressure gradient
    fluxP_s -= -G(inFlow.eps_g) * inFlow.eps_s * edgeData.n;


    fluxP_s *= edgeData.Length;



    // Diffusion of the velocity
    //viscose_x_x = (4.0/3) * myu * du_x_dx  * edgeData.n.x;
    //viscose_y_y = (4.0/3) * myu * du_y_dy * edgeData.n.y;
    Vector<2,double> viscose_diag = ((4.0/3.0) * myu_s * inFlow.eps_s) * Vector<2,double>(du_s_dx.x * edgeData.n.x, du_s_dy.y * edgeData.n.y);

    // Diffusion of the velocity

    //viscose_y_x = myu * (du_x_dy + du_y_dx) * edgeData.n.x;
    //viscose_x_y = myu * (du_x_dy + du_y_dx) * edgeData.n.y;
    Vector<2,double> viscose_side_diag = (myu_s * inFlow.eps_s * (du_s_dx.y + du_s_dy.x)) * Vector<2,double>(edgeData.n.y,edgeData.n.x);

    edgeData.fluxP_s = fluxP_s + viscose_side_diag + viscose_diag;

    edgeData.fluxRho_s = fluxRho_s;

}

inline void MultiphaseFlow::ComputeFluxSolid_outflow(FlowData &innerCellData, EdgeData &edgeData, FVMethod::FluxComputeData &fcData)
{


    double productOf_u_And_n = ((innerCellData.u_s) * edgeData.n);
    double eps_s_e = productOf_u_And_n > 0 ? innerCellData.eps_s : outFlow.eps_s;

    /*
    ** Preparing prtial derivatives wrt
    ** edge orientation
    */
    // computing derivatives of velocity
    Vector<2,double> du_s_dn(0.0, 0.0); // boundary cond
    Vector<2,double> du_s_dt = (PointData.at(fcData.GetCurrEdge().GetPointBIndex()).u_s -
                      PointData.at(fcData.GetCurrEdge().GetPointAIndex()).u_s);

    // transform derivatives to standard base
    Vector<2,double> du_s_dx = du_s_dn * (edgeData.n.x) - (du_s_dt * edgeData.n.y);
    Vector<2,double> du_s_dy = du_s_dt * (edgeData.n.x) + (du_s_dn * edgeData.n.y);



    // flux of density
    double fluxRho_s = -rho_s * eps_s_e *
                        productOf_u_And_n * edgeData.Length;


    // computing the flux of momentum
    Vector<2,double> fluxP_s = - (rho_s) *
                    (productOf_u_And_n) * eps_s_e *
                    (innerCellData.u_s);


    // adding the element of pressure gradient
    fluxP_s -= -G(1.0 - eps_s_e) * eps_s_e * edgeData.n;


    fluxP_s *= edgeData.Length;



    // Diffusion of the velocity
    //viscose_x_x = (4.0/3) * myu * du_x_dx  * edgeData.n.x;
    //viscose_y_y = (4.0/3) * myu * du_y_dy * edgeData.n.y;
    Vector<2,double> viscose_diag = ((4.0/3.0) * myu_s * eps_s_e) * Vector<2,double>(du_s_dx.x * edgeData.n.x, du_s_dy.y * edgeData.n.y);

    // Diffusion of the velocity

    //viscose_y_x = myu * (du_x_dy + du_y_dx) * edgeData.n.x;
    //viscose_x_y = myu * (du_x_dy + du_y_dx) * edgeData.n.y;
    Vector<2,double> viscose_side_diag = (myu_s * eps_s_e * (du_s_dx.y + du_s_dy.x)) * Vector<2,double>(edgeData.n.y,edgeData.n.x);

    edgeData.fluxP_s = fluxP_s + viscose_side_diag + viscose_diag;

    edgeData.fluxRho_s = fluxRho_s;

}

inline void MultiphaseFlow::ComputeFluxSolid_wall(FlowData &innerCellData, EdgeData &edgeData, FVMethod::FluxComputeData &fcData)
{
    // computing derivatives of velocity
    Vector<2,double> du_s_dn = edgeData.LengthOverDist * (fcData.GetCellLeft().GetCellFlag() == INNER ? -1 * innerCellData.u_s : innerCellData.u_s);  // boundary cond

    Vector<2,double> du_s_dt = Vector<2,double>(0.0, 0.0);

    // transform derivatives to standard base
    Vector<2,double> du_s_dx = du_s_dn * (edgeData.n.x) - (du_s_dt * edgeData.n.y);
    Vector<2,double> du_s_dy = du_s_dt * (edgeData.n.x) + (du_s_dn * edgeData.n.y);




    // the flux of momentum
    // is reduced only to pressure gradient
    Vector<2,double> flux = -G(innerCellData.eps_g) * innerCellData.eps_s * edgeData.Length * edgeData.n;


    // Diffusion of the velocity
    Vector<2,double> viscose_diag = ((4.0/3.0) * myu_s * innerCellData.eps_s) * Vector<2,double>(du_s_dx.x * edgeData.n.x, du_s_dy.y * edgeData.n.y);
    Vector<2,double> viscose_side_diag = (myu_s * innerCellData.eps_s * (du_s_dx.y + du_s_dy.x)) * Vector<2,double>(edgeData.n.y,edgeData.n.x);

    edgeData.fluxP_s = flux + viscose_side_diag + viscose_diag;

}




inline double MultiphaseFlow::reg(double x)
{
    return x + 1.0e-7;
}

#endif


#endif // MULTIPHASEFLOW_H

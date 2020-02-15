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
 * @brief reg
 * regularizes number to aviod zero division.
 * Simply adds small non-zero constant (1e-7)
 * @param x
 * number to be regularized
 * @return
 * greatered x
 */
static double reg(double x)
{
    return x + 1.0e-7;
}


/**
 * @brief The FlowData struct
 * flow computation data
 */
struct FlowData {

// gasseous part
    // pressure is dependent on the R_spec and T and rho_g
    //double p;
    double getPressure() const {return rho_g * R_spec * T;}
    void setPressure(const double& pressure){rho_g = pressure / ( R_spec * T );}

    static double R_spec;
    static double T;

    double rho_g;

    /**
     * @brief eps_g
     * volume fraction of gaseous part of the flow
     * holds \epsilon_{s} + \epsilon_{g} = 1
     */
    //double eps_g;
    double getEps_g() const {return 1.0 - eps_s;}
    void setEps_g(const double& eps_g){eps_s = 1.0 - eps_g;}

    /**
     * @brief u_g
     *
     * velocity of gaseous phase
     */
    Vector<2, double> getVelocityGas() const {
        return p_g / reg(rho_g * getEps_g());
    }


    void setVelocityGas(const Vector<2, double>& u_g){
        p_g = u_g * rho_g * getEps_g();
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
    Vector<2, double> getVelocitySolid() const {
        return p_s / reg(rho_s * eps_s);
    }


    void setVelocitySolid(const Vector<2, double>& u_s){
        p_s = u_s * rho_s * eps_s;
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



MAKE_ATTRIBUTE_TRAIT_ARITHMETIC(FlowData, rho_g, eps_s, p_g, p_s)

// TODO only temporary
//MAKE_ATTRIBUTE_TRAIT_IO(FlowData, rho_g, eps_s, p_g, p_s)


MAKE_CUSTOM_ATTRIBUTE_TRAIT(FlowData,
                           "pressure", std::make_pair(&FlowData::getPressure, &FlowData::setPressure),
                           "eps_g", std::make_pair(&FlowData::getEps_g, &FlowData::setEps_g),
                           "rho_g", &FlowData::rho_g,
                           "eps_s", &FlowData::eps_s,
                           "velocity_gas", std::make_pair(&FlowData::getVelocityGas, &FlowData::setVelocityGas),
                           "velocity_solid", std::make_pair(&FlowData::getVelocitySolid, &FlowData::setVelocitySolid)
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

MAKE_ATTRIBUTE_TRAIT(EdgeData, fluxRho_g,fluxP_g,RightCellKoef,LeftCellKoef,n,Length,LengthOverDist)

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

MAKE_ATTRIBUTE_TRAIT(PointData, PointType, cellKoef, u_g, u_s)

struct CellData{
    double invVolume;
};

class MultiphaseFlow {
public:
    MultiphaseFlow() = default;
    ~MultiphaseFlow() = default;

public: // make some types public
    static constexpr unsigned int ProblemDimension = 2;
    using MeshType = UnstructuredMesh<ProblemDimension,size_t, double>;
    using ResultType = FlowData;
public:
    // Structures containing state data

    MeshType mesh;

    MeshDataContainer<std::tuple<CellData,EdgeData,PointData>, 2,1,0> meshData;
private:

    VTKMeshWriter<ProblemDimension> writer;

public:

    MeshDataContainer<std::vector<size_t>, 0> vertToCellCon;

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


    FlowData outFlow;

    Vector<ProblemDimension, double> inFlow_u_s, inFlow_u_g;
    double inFlow_eps_g, inFlow_eps_s;

private:
    /**
      * Functions computing special factors in formula
      */

    /**
     * @brief G
     * @param eps_g
     * @return
     */
    inline double G(double eps_g);

    /**
     * @brief Beta_s
     * @param fd
     * @return
     */
    inline double Beta_s(const FlowData& fd);

    /**
     * @brief C_d
     * @param fd
     * @return
     */
    inline double C_d(const FlowData& fd);

    /**
     * @brief Re_s
     * @param fd
     * @return
     */
    inline double Re_s(const FlowData& fd);

    /**
     * @brief ComputeFluxGas_inner
     * @param leftData
     * @param rightData
     * @param edgeData
     * @param fcData
     */
    inline void ComputeFluxGas_inner(const FlowData& leftData, const FlowData& rightData, EdgeData& edgeData, const MeshType::Face& fcData);

    inline void ComputeFluxGas_inflow(const FlowData& innerCellData, const MeshType::Cell& innerCell,  EdgeData& edgeData, const MeshType::Face& fcData);

    inline void ComputeFluxGas_outflow(const FlowData& innerCellData, EdgeData& edgeData, const MeshType::Face& fcData);

    inline void ComputeFluxGas_wall(const FlowData& innerCellData, EdgeData& edgeData, const MeshType::Face& fcData);

    /**
     * @brief ComputeFluxSolid_inner
     * @param leftData
     * @param rightData
     * @param edgeData
     * @param fcData
     */
    inline void ComputeFluxSolid_inner(const FlowData &leftData, const FlowData &rightData, EdgeData &edgeData, const MeshType::Face& fcData);

    inline void ComputeFluxSolid_inflow(const FlowData& innerCellData, const MeshType::Cell& innerCell,  EdgeData& edgeData, const MeshType::Face& fcData);

    inline void ComputeFluxSolid_outflow(const FlowData& innerCellData, EdgeData& edgeData, const MeshType::Face& fcData);

    inline void ComputeFluxSolid_wall(const FlowData& innerCellData, EdgeData& edgeData, const MeshType::Face& fcData);



private:

    // compute flux over edges for all cells
    //class Flux;



    double FlowModulation(const Vertex<MeshType::meshDimension(), double>& x);


public:


//    void LoadInitialData(std::istream& ist);


public:

    void ComputeFlux(const MeshType::Face& fcData,  const MeshDataContainer<ResultType, ProblemDimension>& compData);

    void ComputeSource(const MeshType::Cell& ccData,
                       MeshDataContainer<ResultType, ProblemDimension>& compData,
                       MeshDataContainer<MultiphaseFlow::ResultType, MultiphaseFlow::ProblemDimension> &outDeltas);

    void setupMeshData(const std::string& fileName);
//
//

    void calculateRHS(
             double time,//time is unused in this problem
             MeshDataContainer<ResultType, ProblemDimension>& compData,
             MeshDataContainer<ResultType, ProblemDimension>& outDeltas);
//
//    void DataWrite(std::ostream& ost);
//
//
//    // Computes velocities in the points on the start of the program
    void InitializePointData();
//
    void UpdateVertexData(const MeshDataContainer<FlowData, MeshType::meshDimension()> &data);
//
//    // Calculates velocities during the computation step for each cell
//    void CalculateVertexData(Mesh::MeshCell &cell, const NumData<FlowData> &cellData);
//
    Type TypeOfCell(const MeshType::Cell& cell);



    template<typename dataType>
    void exportData(double time,
                    MeshDataContainer<dataType, ProblemDimension>& compData) {


        std::ofstream ofile(std::string("MultiphaseFlow") + "_" + std::to_string(time) + ".vtk");
        writer.writeHeader(ofile, std::string("MPF") + std::to_string(time));
        writer.writeToStream(ofile, mesh, MeshDataContainer<MeshNativeType<ProblemDimension>::ElementType, ProblemDimension>(mesh, MeshNativeType<ProblemDimension>::POLYGON));

        VTKMeshDataWriter<ProblemDimension> dataWriter;
        dataWriter.writeToStream(ofile, compData, writer);

        ofile.close();

    }
};



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
        double denominator = 1 / (fd.getEps_g() * d_s * phi_s);
        return (
                    150 * (std::pow(fd.eps_s, 2.0) * myu) * std::pow(denominator, 2.0) +
                    1.75 * (fd.getVelocityGas() - fd.getVelocitySolid()).normEukleid() * fd.rho_g * fd.eps_s * denominator
                );
    } else {
        return (
                    (4.0/3.0) * C_d(fd) * (fd.getVelocityGas() - fd.getVelocitySolid()).normEukleid() * fd.rho_g * fd.eps_s / (d_s * phi_s)
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
            (fd.getVelocityGas() - fd.getVelocitySolid()).normEukleid() * d_s * fd.rho_g * fd.getEps_g() /  myu
          );
}




inline void MultiphaseFlow::ComputeFluxGas_inner(const FlowData &leftData, const FlowData &rightData, EdgeData &edgeData, const MeshType::Face& fcData)
{
    /**
     * @brief eps_g_e
     * volume fraction aproximated at the edge
     */
    double lEps_g = leftData.getEps_g();
    double rEps_g = rightData.getEps_g();

    double eps_g_e = (lEps_g * edgeData.LeftCellKoef + rEps_g * edgeData.RightCellKoef);

    Vector<ProblemDimension, double> ru_g = rightData.getVelocityGas();
    Vector<ProblemDimension, double> lu_g = leftData.getVelocityGas();

    // computing derivatives of velocity
    Vector<ProblemDimension, double> du_g_dn = (ru_g - lu_g) * (edgeData.LengthOverDist);

    Vector<ProblemDimension, double> du_g_dt = (
                meshData.getDataByDim<0>().at(fcData.getVertexBIndex()).u_g -
                meshData.getDataByDim<0>().at(fcData.getVertexAIndex()).u_g);

    //transform diferences to standard base
    Vector<ProblemDimension, double> du_g_dx = du_g_dn * (edgeData.n[0]) - (du_g_dt * edgeData.n[1]);
    Vector<ProblemDimension, double> du_g_dy = du_g_dt * (edgeData.n[0]) + (du_g_dn * edgeData.n[1]);


    double productOf_u_And_n = (((lu_g * edgeData.LeftCellKoef) +
                                 (ru_g * edgeData.RightCellKoef))
                                * edgeData.n);

    // flux of density
    double delta_rho = -(lEps_g * leftData.rho_g * edgeData.LeftCellKoef +
                         rEps_g * rightData.rho_g * edgeData.RightCellKoef) *
                       productOf_u_And_n * edgeData.Length +
                       (rightData.rho_g * rEps_g - leftData.rho_g * lEps_g) * edgeData.LengthOverDist * artifitialDisspation;


    // computing the flux of momentum
    Vector<ProblemDimension,double> fluxP_g = (-productOf_u_And_n) *
                  (leftData.p_g * edgeData.LeftCellKoef + rightData.p_g * edgeData.RightCellKoef);




    // adding the element of pressure gradient
    fluxP_g -= (leftData.getPressure() * edgeData.LeftCellKoef +
                rightData.getPressure() * edgeData.RightCellKoef) * edgeData.n;



    fluxP_g *= edgeData.Length;



    // add artifitial dissipation
    fluxP_g += (edgeData.LengthOverDist * artifitialDisspation * ((rightData.p_g) - (leftData.p_g)));


    //viscose_x_x = (4.0/3) * myu * du_x_dx * edgeData.n[0] * (leftData.eps_g * edgeData.LeftCellKoef + rightData.eps_g * edgeData.RightCellKoef);
    //viscose_y_y = (4.0/3) * myu * du_y_dy * edgeData.n[1] * (leftData.eps_g * edgeData.LeftCellKoef + rightData.eps_g * edgeData.RightCellKoef);
    Vector<ProblemDimension,double> viscose_diag = ((4.0/3) * myu * eps_g_e) * Vector<ProblemDimension,double>{du_g_dx[0] * edgeData.n[0], du_g_dy[1] * edgeData.n[1]};
    // Diffusion of the velocity
    // firstly i will consider diffusion only in the direction

    //viscose_x_y = myu * (du_y_dx + du_x_dy) * edgeDat (leftData.eps_g * edgeData.LeftCellKoef + rightData.eps_g * edgeData.RightCellKoef);
    //viscose_y_x = myu * (du_y_dx + du_x_dy) * edgeData.n[0] * (leftData.eps_g * edgeData.LeftCellKoef + rightData.eps_g * edgeData.RightCellKoef);
    Vector<ProblemDimension,double> viscose_side_diag = (myu * eps_g_e *(du_g_dx[1] + du_g_dy[0])) * Vector<ProblemDimension,double>{edgeData.n[1],edgeData.n[0]};

    // adding of fluxes to cells
    edgeData.fluxP_g = fluxP_g + viscose_diag + viscose_side_diag;

    edgeData.fluxRho_g = delta_rho;


}

inline void MultiphaseFlow::ComputeFluxGas_inflow(const FlowData& innerCellData, const MeshType::Cell& innerCell,  EdgeData& edgeData, const MeshType::Face& fcData)
{

    Vector<ProblemDimension,double> modulatedU = inFlow_u_g;

    modulatedU[1] *= FlowModulation(innerCell.getCenter());

    double productOf_u_And_n = (modulatedU * edgeData.n);

    /*
    ** Preparing partial derivatives wrt
    ** edge orientation
    */

    // computing derivatives of velocity
    Vector<ProblemDimension,double> du_g_dn({0.0, 0.0}); // boundary cond
    Vector<ProblemDimension,double> du_g_dt = (
                meshData.getDataByDim<0>().at(fcData.getVertexBIndex()).u_g -
                meshData.getDataByDim<0>().at(fcData.getVertexAIndex()).u_g);

    // transform derivatives to standard base
    Vector<ProblemDimension,double> du_g_dx = du_g_dn * (edgeData.n[0]) - (du_g_dt * edgeData.n[1]);
    Vector<ProblemDimension,double> du_g_dy = du_g_dt * (edgeData.n[0]) + (du_g_dn * edgeData.n[1]);

    // flux of density
    double delta_rho = -innerCellData.rho_g * inFlow_eps_g *
                        productOf_u_And_n * edgeData.Length;




    // computing the flux of momentum
    Vector<ProblemDimension,double> flux = - (innerCellData.rho_g) *
                    (productOf_u_And_n) * inFlow_eps_g *
                    (modulatedU);


    // adding the element of pressure gradient
    flux -= (innerCellData.getPressure()) * edgeData.n;


    flux *= edgeData.Length;



    // Diffusion of the velocity
    //viscose_x_x = (4.0/3) * myu * du_x_dx  * edgeData.n[0];
    //viscose_y_y = (4.0/3) * myu * du_y_dy * edgeData.n[1];
    Vector<ProblemDimension,double> viscose_diag = ((4.0/3.0) * myu * inFlow_eps_g) * Vector<ProblemDimension,double>{du_g_dx[0] * edgeData.n[0], du_g_dy[1] * edgeData.n[1]};

    // Diffusion of the velocity

    //viscose_y_x = myu * (du_x_dy + du_y_dx) * edgeData.n[0];
    //viscose_x_y = myu * (du_x_dy + du_y_dx) * edgeData.n[1];
    Vector<ProblemDimension,double> viscose_side_diag = (myu * inFlow_eps_g * (du_g_dx[1] + du_g_dy[0])) * Vector<ProblemDimension,double>{edgeData.n[1],edgeData.n[0]};

    edgeData.fluxP_g = flux + viscose_side_diag + viscose_diag;

    edgeData.fluxRho_g = delta_rho;

}

inline void MultiphaseFlow::ComputeFluxGas_outflow(const FlowData& innerCellData, EdgeData& edgeData, const MeshType::Face& fcData)
{
    Vector<ProblemDimension, double> in_u_g = innerCellData.getVelocityGas();
    double productOf_u_And_n = ((in_u_g) * edgeData.n);
    double eps_g_e = productOf_u_And_n > 0 ? innerCellData.getEps_g() : outFlow.getEps_g();

    /*
    ** Preparing partial derivatives wrt
    ** edge orientation
    */

    // computing derivatives of velocity
    Vector<ProblemDimension,double> du_g_dn({0.0, 0.0}); // boundary cond
    Vector<ProblemDimension,double> du_g_dt = (
                meshData.getDataByDim<0>().at(fcData.getVertexBIndex()).u_g -
                meshData.getDataByDim<0>().at(fcData.getVertexAIndex()).u_g);

    // transform derivatives to standard base
    Vector<ProblemDimension,double> du_g_dx = du_g_dn * (edgeData.n[0]) - (du_g_dt * edgeData.n[1]);
    Vector<ProblemDimension,double> du_g_dy = du_g_dt * (edgeData.n[0]) + (du_g_dn * edgeData.n[1]);



    // flux of density
    double delta_rho = -outFlow.rho_g * eps_g_e *
                        productOf_u_And_n * edgeData.Length;


    // computing the flux of momentum
    Vector<ProblemDimension,double> flux = - (outFlow.rho_g *
                    productOf_u_And_n * eps_g_e) *
                    (in_u_g);


    // adding the element of pressure gradient
    flux -= (outFlow.getPressure()) * edgeData.n;

    flux *= edgeData.Length;



    // Diffusion of the velocity
    //viscose_x_x = (4.0/3) * myu * du_x_dx  * edgeData.n[0];
    //viscose_y_y = (4.0/3) * myu * du_y_dy * edgeData.n[1];
    Vector<ProblemDimension,double> viscose_diag = ((4.0/3.0) * myu * eps_g_e) * Vector<ProblemDimension,double>({du_g_dx[0] * edgeData.n[0], du_g_dy[1] * edgeData.n[1]});

    // Diffusion of the velocity

    //viscose_y_x = myu * (du_x_dy + du_y_dx) * edgeData.n[0];
    //viscose_x_y = myu * (du_x_dy + du_y_dx) * edgeData.n[1];
    Vector<ProblemDimension,double> viscose_side_diag = (myu * eps_g_e * (du_g_dx[1] + du_g_dy[0])) * Vector<ProblemDimension,double>({edgeData.n[1],edgeData.n[0]});


    edgeData.fluxP_g = flux + viscose_side_diag + viscose_diag;

    edgeData.fluxRho_g = delta_rho;

}

inline void MultiphaseFlow::ComputeFluxGas_wall(const FlowData& innerCellData, EdgeData& edgeData, const MeshType::Face& fcData)
{
    // computing derivatives of velocity
    Vector<ProblemDimension,double> du_g_dn = edgeData.LengthOverDist * (fcData.getCellLeftIndex() < BOUNDARY_INDEX(size_t) ?
                                                              -1.0 * innerCellData.getVelocityGas() : innerCellData.getVelocityGas());  // boundary cond

    Vector<ProblemDimension,double> du_g_dt = Vector<ProblemDimension,double>({0.0, 0.0});

    // transform derivatives to standard base
    Vector<ProblemDimension,double> du_g_dx = du_g_dn * (edgeData.n[0]) - (du_g_dt * edgeData.n[1]);
    Vector<ProblemDimension,double> du_g_dy = du_g_dt * (edgeData.n[0]) + (du_g_dn * edgeData.n[1]);




    // the flux of momentum
    // is reduced only to pressure gradient
    Vector<ProblemDimension,double> flux = -(innerCellData.getPressure()) * edgeData.n;


    flux *= edgeData.Length;

    //prepare eps_g
    double iEps_g = innerCellData.getEps_g();

    // Diffusion of the velocity
    Vector<ProblemDimension,double> viscose_diag = ((4.0/3.0) * myu * iEps_g) * Vector<ProblemDimension,double>({du_g_dx[0] * edgeData.n[0], du_g_dy[1] * edgeData.n[1]});
    Vector<ProblemDimension,double> viscose_side_diag = (myu * iEps_g * (du_g_dx[1] + du_g_dy[0])) * Vector<ProblemDimension,double>({edgeData.n[1],edgeData.n[0]});



    edgeData.fluxP_g = flux + viscose_side_diag + viscose_diag;

}








inline void MultiphaseFlow::ComputeFluxSolid_inner(const FlowData &leftData, const FlowData &rightData, EdgeData &edgeData, const MeshType::Face& fcData)
{
    /**
     * @brief eps_s_e
     * volume fraction aproximated at the edge
     */
    double eps_s_e = (leftData.eps_s * edgeData.LeftCellKoef + rightData.eps_s * edgeData.RightCellKoef);


    Vector<ProblemDimension, double> ru_s = rightData.getVelocitySolid();
    Vector<ProblemDimension, double> lu_s = leftData.getVelocitySolid();

    Vector<ProblemDimension,double> du_s_dn = (ru_s - lu_s) * edgeData.LengthOverDist;

    Vector<ProblemDimension,double> du_s_dt = (
                meshData.getDataByDim<0>().at(fcData.getVertexBIndex()).u_s -
                meshData.getDataByDim<0>().at(fcData.getVertexAIndex()).u_s);

    //transform diferences to standard base
    Vector<ProblemDimension,double> du_s_dx = du_s_dn * (edgeData.n[0]) - (du_s_dt * edgeData.n[1]);

    Vector<ProblemDimension,double> du_s_dy = du_s_dt * (edgeData.n[0]) + (du_s_dn * edgeData.n[1]);

    double product_of_u_s_and_n = (((lu_s * edgeData.LeftCellKoef) +
                                  (ru_s * edgeData.RightCellKoef)) *
                                  edgeData.n);


    // flux of mass
    double fluxRho_s = (-rho_s * eps_s_e * edgeData.Length * product_of_u_s_and_n  +
                       (rightData.eps_s - leftData.eps_s) * rho_s * edgeData.LengthOverDist * artifitialDisspation) ;


    // computing the flux of momentum
    Vector<ProblemDimension,double> fluxP_s = -product_of_u_s_and_n *
                     (leftData.p_s * edgeData.LeftCellKoef + rightData.p_s * edgeData.RightCellKoef);


    // adding the element of "pressure" gradient
    fluxP_s -= G(leftData.getEps_g() * edgeData.LeftCellKoef + rightData.getEps_g() * edgeData.RightCellKoef) * eps_s_e *
              edgeData.n;

    fluxP_s *= edgeData.Length;

    // add artifitial dissipation
    // TODO this is the critical point it is necessary to look into this problem
    // maybe the artifitial disipation must be raise
    fluxP_s += (edgeData.LengthOverDist * artifitialDisspation ) * (rightData.p_s - leftData.p_s);




    //viscose_x_x = (4.0/3) * myu * du_x_dx * edgeData.n[0] * (leftData.eps_g * edgeData.LeftCellKoef + rightData.eps_g * edgeData.RightCellKoef);
    //viscose_y_y = (4.0/3) * myu * du_y_dy * edgeData.n[1] * (leftData.eps_g * edgeData.LeftCellKoef + rightData.eps_g * edgeData.RightCellKoef);
    //Vector<ProblemDimension,double> viscose_diag = ((4.0/3) * myu_s * eps_s_e) * Vector<ProblemDimension,double>(du_s_dx[0] * edgeData.n[0], du_s_dy[1] * edgeData.n[1]);
    fluxP_s += ((4.0/3) * myu_s * eps_s_e) * Vector<ProblemDimension,double>{du_s_dx[0] * edgeData.n[0], du_s_dy[1] * edgeData.n[1]};

    //Vector<ProblemDimension,double> viscose_diag = (4.0/3) * myu * du_s_dy * edgeData.n[1] * (leftData.eps_g * edgeData.LeftCellKoef + rightData.eps_g * edgeData.RightCellKoef);


    //viscose_x_y = myu * (du_y_dx + du_x_dy) * edgeData.n[1] * (leftData.eps_g * edgeData.LeftCellKoef + rightData.eps_g * edgeData.RightCellKoef);
    //viscose_y_x = myu * (du_y_dx + du_x_dy) * edgeData.n[0] * (leftData.eps_g * edgeData.LeftCellKoef + rightData.eps_g * edgeData.RightCellKoef);
    //Vector<ProblemDimension,double> viscose_side_diag = (myu_s * eps_s_e *(du_s_dx[1] + du_s_dy[0])) * Vector<ProblemDimension,double>(edgeData.n[1],edgeData.n[0]);
    fluxP_s += (myu_s * eps_s_e *(du_s_dx[1] + du_s_dy[0])) * Vector<ProblemDimension,double>{edgeData.n[1],edgeData.n[0]};

    // storing result at the edge
    edgeData.fluxP_s = fluxP_s; // + viscose_diag + viscose_side_diag;

    edgeData.fluxRho_s = fluxRho_s;
}


inline void MultiphaseFlow::ComputeFluxSolid_inflow(const FlowData& innerCellData, const MeshType::Cell& innerCell,  EdgeData& edgeData, const MeshType::Face& fcData)
{
    (void)innerCellData;

    Vector<ProblemDimension,double> modulatedU = inFlow_u_s;

    modulatedU[1] *= FlowModulation(innerCell.getCenter());

    double productOf_u_And_n = (modulatedU * edgeData.n);

    /*
    ** Preparing partial derivatives wrt
    ** edge orientation
    */

    // computing derivatives of velocity
    Vector<ProblemDimension,double> du_s_dn{0.0, 0.0}; // boundary cond
    Vector<ProblemDimension,double> du_s_dt = (
                meshData.getDataByDim<0>().at(fcData.getVertexBIndex()).u_s -
                meshData.getDataByDim<0>().at(fcData.getVertexAIndex()).u_s);

    // transform derivatives to standard base
    Vector<ProblemDimension,double> du_s_dx = du_s_dn * (edgeData.n[0]) - (du_s_dt * edgeData.n[1]);
    Vector<ProblemDimension,double> du_s_dy = du_s_dt * (edgeData.n[0]) + (du_s_dn * edgeData.n[1]);

    // flux of density
    double fluxRho_s = -rho_s * inFlow_eps_s *
                        productOf_u_And_n * edgeData.Length;




    // computing the flux of momentum
    Vector<ProblemDimension,double> fluxP_s = - (rho_s) *
                    (productOf_u_And_n) * inFlow_eps_s *
                    (modulatedU);

    // adding the element of pressure gradient
    fluxP_s -= -G(inFlow_eps_g) * inFlow_eps_s * edgeData.n;


    fluxP_s *= edgeData.Length;



    // Diffusion of the velocity
    //viscose_x_x = (4.0/3) * myu * du_x_dx  * edgeData.n[0];
    //viscose_y_y = (4.0/3) * myu * du_y_dy * edgeData.n[1];
    Vector<ProblemDimension,double> viscose_diag = ((4.0/3.0) * myu_s * inFlow_eps_s) * Vector<ProblemDimension,double>{du_s_dx[0] * edgeData.n[0], du_s_dy[1] * edgeData.n[1]};

    // Diffusion of the velocity

    //viscose_y_x = myu * (du_x_dy + du_y_dx) * edgeData.n[0];
    //viscose_x_y = myu * (du_x_dy + du_y_dx) * edgeData.n[1];
    Vector<ProblemDimension,double> viscose_side_diag = (myu_s * inFlow_eps_s * (du_s_dx[1] + du_s_dy[0])) * Vector<ProblemDimension,double>{edgeData.n[1],edgeData.n[0]};

    edgeData.fluxP_s = fluxP_s + viscose_side_diag + viscose_diag;

    edgeData.fluxRho_s = fluxRho_s;

}

inline void MultiphaseFlow::ComputeFluxSolid_outflow(const FlowData& innerCellData, EdgeData& edgeData, const MeshType::Face& fcData)
{

    Vector<ProblemDimension, double> in_u_s = innerCellData.getVelocitySolid();

    double productOf_u_And_n = ((in_u_s) * edgeData.n);
    double eps_s_e = productOf_u_And_n > 0 ? innerCellData.eps_s : outFlow.eps_s;

    /*
    ** Preparing partial derivatives wrt
    ** edge orientation
    */
    // computing derivatives of velocity
    Vector<ProblemDimension,double> du_s_dn({0.0, 0.0}); // boundary cond
    Vector<ProblemDimension,double> du_s_dt = (
                meshData.getDataByDim<0>().at(fcData.getVertexBIndex()).u_s -
                meshData.getDataByDim<0>().at(fcData.getVertexAIndex()).u_s);

    // transform derivatives to standard base
    Vector<ProblemDimension,double> du_s_dx = du_s_dn * (edgeData.n[0]) - (du_s_dt * edgeData.n[1]);
    Vector<ProblemDimension,double> du_s_dy = du_s_dt * (edgeData.n[0]) + (du_s_dn * edgeData.n[1]);



    // flux of density
    double fluxRho_s = -rho_s * eps_s_e *
                        productOf_u_And_n * edgeData.Length;


    // computing the flux of momentum
    Vector<ProblemDimension,double> fluxP_s = - (rho_s) *
                    (productOf_u_And_n) * eps_s_e *
                    (in_u_s);


    // adding the element of pressure gradient
    fluxP_s -= -G(1.0 - eps_s_e) * eps_s_e * edgeData.n;


    fluxP_s *= edgeData.Length;



    // Diffusion of the velocity
    //viscose_x_x = (4.0/3) * myu * du_x_dx  * edgeData.n[0];
    //viscose_y_y = (4.0/3) * myu * du_y_dy * edgeData.n[1];
    Vector<ProblemDimension,double> viscose_diag = ((4.0/3.0) * myu_s * eps_s_e) * Vector<ProblemDimension,double>({du_s_dx[0] * edgeData.n[0], du_s_dy[1] * edgeData.n[1]});

    // Diffusion of the velocity

    //viscose_y_x = myu * (du_x_dy + du_y_dx) * edgeData.n[0];
    //viscose_x_y = myu * (du_x_dy + du_y_dx) * edgeData.n[1];
    Vector<ProblemDimension,double> viscose_side_diag = (myu_s * eps_s_e * (du_s_dx[1] + du_s_dy[0])) * Vector<ProblemDimension,double>({edgeData.n[1],edgeData.n[0]});

    edgeData.fluxP_s = fluxP_s + viscose_side_diag + viscose_diag;

    edgeData.fluxRho_s = fluxRho_s;

}

inline void MultiphaseFlow::ComputeFluxSolid_wall(const FlowData& innerCellData, EdgeData& edgeData, const MeshType::Face& fcData)
{
    // computing derivatives of velocity
    Vector<ProblemDimension,double> du_s_dn = edgeData.LengthOverDist * (fcData.getCellLeftIndex() < BOUNDARY_INDEX(size_t) ?
                                                                             -1.0 * innerCellData.getVelocitySolid() : innerCellData.getVelocitySolid());  // boundary cond

    Vector<ProblemDimension,double> du_s_dt = {0.0, 0.0};

    // transform derivatives to standard base
    Vector<ProblemDimension,double> du_s_dx = du_s_dn * (edgeData.n[0]) - (du_s_dt * edgeData.n[1]);
    Vector<ProblemDimension,double> du_s_dy = du_s_dt * (edgeData.n[0]) + (du_s_dn * edgeData.n[1]);




    // the flux of momentum
    // is reduced only to pressure gradient
    Vector<ProblemDimension,double> flux = -G(innerCellData.getEps_g()) * innerCellData.eps_s * edgeData.Length * edgeData.n;


    // Diffusion of the velocity
    Vector<ProblemDimension,double> viscose_diag = ((4.0/3.0) * myu_s * innerCellData.eps_s) * Vector<ProblemDimension,double>({du_s_dx[0] * edgeData.n[0], du_s_dy[1] * edgeData.n[1]});
    Vector<ProblemDimension,double> viscose_side_diag = (myu_s * innerCellData.eps_s * (du_s_dx[1] + du_s_dy[0])) * Vector<ProblemDimension,double>({edgeData.n[1],edgeData.n[0]});

    edgeData.fluxP_s = flux + viscose_side_diag + viscose_diag;

}







#endif // MULTIPHASEFLOW_H

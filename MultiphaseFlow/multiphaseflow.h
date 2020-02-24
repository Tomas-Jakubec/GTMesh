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



/*
  Helper functions to work with temporary tensors
*/
namespace Impl {
template <unsigned int dim, typename Real, unsigned int Index = 0>
typename std::enable_if<(Index == dim - 1)>::type
tensorProduct(Vector<dim, Vector<dim, Real>>& res, const Vector<dim, Real>& v1, const Vector<dim, Real>& v2){
    res[Index] = v1[Index] * v2;
}


template <unsigned int dim, typename Real, unsigned int Index = 0>
typename std::enable_if<(Index < dim - 1)>::type
tensorProduct(Vector<dim, Vector<dim, Real>>& res, const Vector<dim, Real>& v1, const Vector<dim, Real>& v2){
    res[Index] = v1[Index] * v2;
    Impl::tensorProduct<dim, Real, Index + 1>(res, v1, v2);
}
}

/**
 * @brief Calculates tensor product of v1 and v2.
 * The result is a matrix of type Vector<dim, Vector<dim, Real>>
 * @param v1
 * @param v2
 * @return
 */
template <unsigned int dim, typename Real>
Vector<dim, Vector<dim, Real>> tensorProduct(const Vector<dim, Real> &v1, const Vector<dim, Real> &v2){
    Vector<dim, Vector<dim, Real>> res;
    Impl::tensorProduct(res, v1, v2);
    return res;
}

namespace Impl {
template <unsigned int dim, typename Real, unsigned int Index = 0>
typename std::enable_if<(Index == dim - 1)>::type
tensorVectorProduct(Vector<dim, Real>& res, const Vector<dim, Vector<dim, Real>>& t1, const Vector<dim, Real>& v2){
    res[Index] = t1[Index] * v2;
}


template <unsigned int dim, typename Real, unsigned int Index = 0>
typename std::enable_if<(Index < dim - 1)>::type
tensorVectorProduct(Vector<dim, Real>& res, const Vector<dim, Vector<dim, Real>>& t1, const Vector<dim, Real>& v2){
    res[Index] = t1[Index] * v2; // scalar product
    Impl::tensorVectorProduct<dim, Real, Index + 1>(res, t1, v2);
}
}

/**
 * @brief Calculates application of tensor t1 to v2.
 * The result is a vector of type Vector<dim, Real>
 * @param t1
 * @param v2
 * @return
 */
template <unsigned int dim, typename Real>
Vector<dim, Real> tensorVectorProduct(const Vector<dim, Vector<dim, Real>> &t1, const Vector<dim, Real> &v2){
    Vector<dim, Real> res;
    Impl::tensorVectorProduct(res, t1, v2);
    return res;
}

/**
 * @brief The FlowData struct
 * flow compuatation data
 */
template<unsigned int Dim>
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
    Vector<Dim, double> getVelocityGas() const {
        return p_g / reg(rho_g * getEps_g());
    }


    void setVelocityGas(const Vector<Dim, double>& u_g){
        p_g = u_g * rho_g * getEps_g();
    }

    /**
     * @brief p_g
     *
     * momentum of gaseous phase
     */
    Vector<Dim,double> p_g;



// solid part
    /**
     * @brief u_s
     *
     * velocity of solid component
     */
    Vector<Dim, double> getVelocitySolid() const {
        return p_s / reg(rho_s * eps_s);
    }


    void setVelocitySolid(const Vector<Dim, double>& u_s){
        p_s = u_s * rho_s * eps_s;
    }

    static double rho_s;

    /**
     * @brief p_s
     *
     * momentum of solid phase
     */
    Vector<Dim,double> p_s;


    /**
     * @brief eps_g
     * volume fraction of solid part of the flow
     * holds \epsilon_{s} + \epsilon_{g} = 1
     */
    double eps_s;



};

template<unsigned int Dim>
double FlowData<Dim>::rho_s = 0;
template<unsigned int Dim>
double FlowData<Dim>::R_spec = 0;
template<unsigned int Dim>
double FlowData<Dim>::T = 0;

MAKE_ATTRIBUTE_TRAIT_ARITHMETIC(FlowData<2>, rho_g, eps_s, p_g, p_s)

// TODO only temporary
//MAKE_ATTRIBUTE_TRAIT_IO(FlowData, rho_g, eps_s, p_g, p_s)


MAKE_CUSTOM_ATTRIBUTE_TRAIT(FlowData<2>,
                           "pressure", std::make_pair(&FlowData<2>::getPressure, &FlowData<2>::setPressure),
                           "eps_g", std::make_pair(&FlowData<2>::getEps_g, &FlowData<2>::setEps_g),
                           "rho_g", &FlowData<2>::rho_g,
                           "eps_s", &FlowData<2>::eps_s,
                           "velocity_gas", std::make_pair(&FlowData<2>::getVelocityGas, &FlowData<2>::setVelocityGas),
                           "velocity_solid", std::make_pair(&FlowData<2>::getVelocitySolid, &FlowData<2>::setVelocitySolid)
                            )

MAKE_ATTRIBUTE_TRAIT_ARITHMETIC(FlowData<3>, rho_g, eps_s, p_g, p_s)

// TODO only temporary
//MAKE_ATTRIBUTE_TRAIT_IO(FlowData, rho_g, eps_s, p_g, p_s)


MAKE_CUSTOM_ATTRIBUTE_TRAIT(FlowData<3>,
                           "pressure", std::make_pair(&FlowData<3>::getPressure, &FlowData<3>::setPressure),
                           "eps_g", std::make_pair(&FlowData<3>::getEps_g, &FlowData<3>::setEps_g),
                           "rho_g", &FlowData<3>::rho_g,
                           "eps_s", &FlowData<3>::eps_s,
                           "velocity_gas", std::make_pair(&FlowData<3>::getVelocityGas, &FlowData<3>::setVelocityGas),
                           "velocity_solid", std::make_pair(&FlowData<3>::getVelocitySolid, &FlowData<3>::setVelocitySolid)
                            )



template<unsigned int Dim>
struct FaceData {
    double LengthOverDist;
    double Length;
    Vector<Dim, double> n;

    // koeficients of konvex combination of cells
    // holds that LeftCellKoef + RightCellKoef = 1
    double LeftCellKoef;
    double RightCellKoef;


    /**
     * Next variables are supposed to store fluxes over the edge
     * in order to simply parallelize the compuatation
     */

    /**
     * @brief fluxP_g
     *
     * flux of momentum of gaseous phase over the edge
     */
    Vector<Dim, double> fluxP_g;


    /**
     * @brief fluxP_s
     *
     * flux of momentum of solid phase over the edge
     */
    Vector<Dim, double> fluxP_s;


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


    /**
     * @brief gradient of velocity of gas
     */
    Vector<Dim, Vector<Dim, double>> grad_u_g;

    /**
     * @brief gradient of velocity of gas
     */
    Vector<Dim, Vector<Dim, double>> grad_u_s;
};

MAKE_ATTRIBUTE_TRAIT(FaceData<2>, fluxRho_g,fluxP_g,RightCellKoef,LeftCellKoef,n,Length,LengthOverDist)
MAKE_ATTRIBUTE_TRAIT(FaceData<3>, fluxRho_g,fluxP_g,RightCellKoef,LeftCellKoef,n,Length,LengthOverDist)
// Enumerator Type for expessing other cell and point position
enum Type{
    INNER = 1,
    WALL = 2,
    INFLOW = 4,
    OUTFLOW = 8,
};




template<unsigned int Dim>
struct CellData{
    double invVolume;
    Vector<Dim, Vector<Dim, double>> grad_u_g;

    Vector<Dim, Vector<Dim, double>> grad_u_s;
};

MAKE_ATTRIBUTE_TRAIT(CellData<3>, invVolume);



class MultiphaseFlow {
public:
    MultiphaseFlow() = default;
    ~MultiphaseFlow() = default;

public: // make some types public
    static constexpr unsigned int ProblemDimension = 2;
    using MeshType = UnstructuredMesh<ProblemDimension,size_t, double>;
    using ResultType = FlowData<ProblemDimension>;

    using Cell = typename MeshType::Cell;
    using Face = typename MeshType::Face;
public:
    // Structures containing state data

    MeshType mesh;

    MeshDataContainer<std::tuple<CellData<ProblemDimension>,FaceData<ProblemDimension>>, ProblemDimension,ProblemDimension - 1> meshData;
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


    FlowData<ProblemDimension> outFlow;

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
    inline double Beta_s(const FlowData<ProblemDimension>& fd);

    /**
     * @brief C_d
     * @param fd
     * @return
     */
    inline double C_d(const FlowData<ProblemDimension>& fd);

    /**
     * @brief Re_s
     * @param fd
     * @return
     */
    inline double Re_s(const FlowData<ProblemDimension>& fd);

    /**
     * @brief ComputeFluxGas_inner
     * @param leftData
     * @param rightData
     * @param edgeData
     * @param fcData
     */
    inline void ComputeFluxGas_inner(const FlowData<ProblemDimension>& leftData, const FlowData<ProblemDimension>& rightData, FaceData<ProblemDimension>& edgeData, const Face& fcData);

    inline void ComputeFluxGas_inflow(const FlowData<ProblemDimension>& innerCellData, const Cell& innerCell,  FaceData<ProblemDimension>& edgeData, const Face& fcData);

    inline void ComputeFluxGas_outflow(const FlowData<ProblemDimension>& innerCellData, FaceData<ProblemDimension>& edgeData, const Face& fcData);

    inline void ComputeFluxGas_wall(const FlowData<ProblemDimension>& innerCellData, FaceData<ProblemDimension>& edgeData, const Face& fcData);

    inline void ComputeViscousFluxGas_inner(const FlowData<ProblemDimension>& leftData, const FlowData<ProblemDimension>& rightData, FaceData<ProblemDimension>& edgeData, const Face& fcData);

    inline void ComputeViscousFluxGas_inflow(const Cell& innerCell,  const Face& face);

    inline void ComputeViscousFluxGas_outflow(const Cell& innerCell,  const Face& face);

    inline void ComputeViscousFluxGas_wall(const Cell& innerCell,  const Face& face);

    /**
     * @brief ComputeFluxSolid_inner
     * @param leftData
     * @param rightData
     * @param edgeData
     * @param fcData
     */
    inline void ComputeFluxSolid_inner(const FlowData<ProblemDimension> &leftData, const FlowData<ProblemDimension> &rightData, FaceData<ProblemDimension> &edgeData, const Face& fcData);

    inline void ComputeFluxSolid_inflow(const FlowData<ProblemDimension>& innerCellData, const Cell& innerCell,  FaceData<ProblemDimension>& edgeData, const Face& fcData);

    inline void ComputeFluxSolid_outflow(const FlowData<ProblemDimension>& innerCellData, FaceData<ProblemDimension>& edgeData, const Face& fcData);

    inline void ComputeFluxSolid_wall(const FlowData<ProblemDimension>& innerCellData, FaceData<ProblemDimension>& edgeData, const Face& fcData);

    inline void ComputeViscousFluxSolid_inner(const FlowData<ProblemDimension>& leftData, const FlowData<ProblemDimension>& rightData, FaceData<ProblemDimension>& edgeData, const Face& fcData);

    inline void ComputeViscousFluxSolid_inflow(const Cell& innerCell,  const Face& face);

    inline void ComputeViscousFluxSolid_outflow(const Cell& innerCell,  const Face& face);

    inline void ComputeViscousFluxSolid_wall(const Cell& innerCell,  const Face& face);


private:

    // compute flux over edges for all cells
    //class Flux;



    double FlowModulation(const Vertex<MeshType::meshDimension(), double>& x);


public:


//    void LoadInitialData(std::istream& ist);


public:

    void ComputeFlux(const Face& fcData,  const MeshDataContainer<ResultType, ProblemDimension>& compData);


    void ComputeViscousFlux(const Face& fcData,  const MeshDataContainer<ResultType, ProblemDimension>& compData);


    /**
     * @brief Calculates gradient of velocity using Gauss-Green th. in a cell.
     * @param cell
     */
    inline void ComupteGradU(const Cell& cell);

    void ComputeSource(const Cell& ccData,
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
//
//    // Calculates velocities during the compuatation step for each cell
//    void CalculateVertexData(Mesh::MeshCell &cell, const NumData<FlowData> &cellData);
//
    Type TypeOfCell(const Cell& cell);



    template<typename dataType>
    void exportData(double time,
                    MeshDataContainer<dataType, ProblemDimension>& compData) {


        std::ofstream ofile(std::string("MultiphaseFlow") + "_" + std::to_string(time) + ".vtk");
        writer.writeHeader(ofile, std::string("MPF") + std::to_string(time));
        auto polyType = typename std::conditional<ProblemDimension == 2, MeshNativeType<2>::ElementType, MeshNativeType<3>::ElementType>::type(ProblemDimension == 2 ? MeshNativeType<2>::POLYGON : MeshNativeType<3>::HEXAHEDRON);
        writer.writeToStream(ofile, mesh, MeshDataContainer<MeshNativeType<ProblemDimension>::ElementType, ProblemDimension>(mesh,  polyType));

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






inline double MultiphaseFlow::Beta_s(const FlowData<ProblemDimension> &fd)
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






inline double MultiphaseFlow::C_d(const FlowData<ProblemDimension> &fd)
{
    double re_s = reg(Re_s(fd));

    if (re_s <= 1000) {
        return (24/re_s)*(1.0 + 0.15 * std::pow(re_s, 0.687));
    } else {
        return 0.44;
    }
}






inline double MultiphaseFlow::Re_s(const FlowData<ProblemDimension> &fd)
{
    // multiplying by inverted value of myu may make code faster
    return(
            (fd.getVelocityGas() - fd.getVelocitySolid()).normEukleid() * d_s * fd.rho_g * fd.getEps_g() /  myu
          );
}




inline void MultiphaseFlow::ComputeFluxGas_inner(const FlowData<ProblemDimension> &leftData, const FlowData<ProblemDimension> &rightData, FaceData<ProblemDimension> &edgeData, const Face&)
{

    auto edge_u_g = (leftData.getVelocityGas() * edgeData.LeftCellKoef) + (rightData.getVelocityGas() * edgeData.RightCellKoef);
    double product_of_u_and_n = (edge_u_g * edgeData.n);

    // prepare upwind values
    const FlowData<ProblemDimension>& faceVal = product_of_u_and_n > 0 ? leftData : rightData;

    // update the velocity at the face/edge
    //productOf_u_And_n = faceVal.getVelocityGas() * edgeData.n;

    // flux of density
    double delta_rho = - faceVal.getEps_g() * faceVal.rho_g * product_of_u_and_n * edgeData.Length +
                       (rightData.rho_g * rightData.getEps_g() - leftData.rho_g * leftData.getEps_g()) * edgeData.LengthOverDist * artifitialDisspation;


    // computing the flux of momentum
    Vector<ProblemDimension,double> fluxP_g = (-product_of_u_and_n) * faceVal.p_g;

    // adding the element of pressure gradient
    fluxP_g -= (leftData.getPressure() * edgeData.LeftCellKoef + rightData.getPressure() * edgeData.RightCellKoef) * edgeData.n;
    // multiply by the face/edge measure
    fluxP_g *= edgeData.Length;

    // add artifitial dissipation
    fluxP_g += (edgeData.LengthOverDist * artifitialDisspation * ((rightData.p_g) - (leftData.p_g)));


    // compuatation of grad(u)
    edgeData.grad_u_g = tensorProduct(edge_u_g, edgeData.n) * edgeData.Length;

    // adding of fluxes to cells
    edgeData.fluxP_g = fluxP_g;

    edgeData.fluxRho_g = delta_rho;


}

inline void MultiphaseFlow::ComputeFluxGas_inflow(const FlowData<ProblemDimension>& innerCellData, const Cell& ,  FaceData<ProblemDimension>& edgeData, const Face& face)
{

    Vector<ProblemDimension,double> modulatedU = inFlow_u_g;

    modulatedU *= FlowModulation(face.getCenter());

    double product_of_u_and_n = (modulatedU * edgeData.n);


    // flux of density
    double delta_rho = -innerCellData.rho_g * inFlow_eps_g *
                        product_of_u_and_n * edgeData.Length;




    // computing the flux of momentum
    Vector<ProblemDimension,double> flux = - (innerCellData.rho_g) *
                    (product_of_u_and_n) * inFlow_eps_g *
                    (modulatedU);


    // adding the element of pressure gradient
    flux -= (innerCellData.getPressure()) * edgeData.n;


    flux *= edgeData.Length;



    // compuatation of grad(u)
    edgeData.grad_u_g = tensorProduct(modulatedU, edgeData.n) * edgeData.Length;

    edgeData.fluxP_g = flux;// + viscose_side_diag + viscose_diag;

    edgeData.fluxRho_g = delta_rho;

}

inline void MultiphaseFlow::ComputeFluxGas_outflow(const FlowData<ProblemDimension>& innerCellData, FaceData<ProblemDimension>& edgeData, const Face&)
{
    Vector<ProblemDimension, double> in_u_g = innerCellData.getVelocityGas();
    double product_of_u_and_n = ((in_u_g) * edgeData.n);
    double eps_g_e = product_of_u_and_n > 0 ? innerCellData.getEps_g() : outFlow.getEps_g();



    // flux of density
    double delta_rho = -outFlow.rho_g * eps_g_e *
                        product_of_u_and_n * edgeData.Length;


    // computing the flux of momentum
    Vector<ProblemDimension,double> flux = - (outFlow.rho_g *
                    product_of_u_and_n * eps_g_e) *
                    (in_u_g);


    // adding the element of pressure gradient
    flux -= (outFlow.getPressure()) * edgeData.n;

    flux *= edgeData.Length;


    // compuatation of grad(u)
    edgeData.grad_u_g = tensorProduct(in_u_g, edgeData.n) * edgeData.Length;

    edgeData.fluxP_g = flux;// + viscose_side_diag + viscose_diag;

    edgeData.fluxRho_g = delta_rho;

}

inline void MultiphaseFlow::ComputeFluxGas_wall(const FlowData<ProblemDimension>& innerCellData, FaceData<ProblemDimension>& edgeData, const Face&)
{





    // the flux of momentum
    // is reduced only to pressure gradient
    Vector<ProblemDimension,double> flux = -(innerCellData.getPressure()) * edgeData.n;


    flux *= edgeData.Length;

    //prepare eps_g

    // compuatation of grad(u)
    edgeData.grad_u_g = {};

    edgeData.fluxP_g = flux;// + viscose_side_diag + viscose_diag;
    edgeData.fluxRho_g = 0;
}







void MultiphaseFlow::ComputeViscousFluxGas_inner(const FlowData<ProblemDimension> &leftData, const FlowData<ProblemDimension> &rightData, FaceData<ProblemDimension> &edgeData, const Face &fcData)
{

    double lEps_g = leftData.getEps_g();
    double rEps_g = rightData.getEps_g();

    double eps_g_e = (lEps_g * edgeData.LeftCellKoef + rEps_g * edgeData.RightCellKoef);

    auto edge_grad_u_g = (meshData.getDataByDim<ProblemDimension>()[fcData.getCellLeftIndex()].grad_u_g * edgeData.LeftCellKoef + meshData.getDataByDim<ProblemDimension>()[fcData.getCellRightIndex()].grad_u_g * edgeData.RightCellKoef);
    decltype (edge_grad_u_g) viscousT;

    for (unsigned int i = 0; i < edge_grad_u_g.size(); i++) {
        for (unsigned int j = 0; j < edge_grad_u_g.size(); j++) {
            // grad_u + transposed grad_u
            viscousT[i][j] = (edge_grad_u_g[i][j] + edge_grad_u_g[i][j]) * eps_g_e * myu;
        }
    }

    for (unsigned int i = 0; i < viscousT.size(); i++) {
        viscousT[i][i] *= 2.0 / 3.0;
    }

    edgeData.fluxP_g += tensorVectorProduct(viscousT, edgeData.n);

}

void MultiphaseFlow::ComputeViscousFluxGas_inflow(const Cell &innerCell, const Face& face)
{

    double eps_g_e = inFlow_eps_g;

    auto edge_grad_u_g = (meshData[innerCell].grad_u_g);
    decltype (edge_grad_u_g) viscousT;

    for (unsigned int i = 0; i < edge_grad_u_g.size(); i++) {
        for (unsigned int j = 0; j < edge_grad_u_g.size(); j++) {
            // grad_u + transposed grad_u
            viscousT[i][j] = (edge_grad_u_g[i][j] + edge_grad_u_g[i][j]) * eps_g_e * myu;
        }
    }

    for (unsigned int i = 0; i < viscousT.size(); i++) {
        viscousT[i][i] *= 2.0 / 3.0;
    }

    meshData[face].fluxP_g += tensorVectorProduct(viscousT, meshData[face].n);

}

void MultiphaseFlow::ComputeViscousFluxGas_outflow(const Cell& innerCell,  const Face& face)
{
    ComputeViscousFluxGas_inflow(innerCell, face);
}

void MultiphaseFlow::ComputeViscousFluxGas_wall(const Cell& innerCell,  const Face& face)
{
    ComputeViscousFluxGas_inflow(innerCell, face);
}








inline void MultiphaseFlow::ComputeFluxSolid_inner(const FlowData<ProblemDimension> &leftData, const FlowData<ProblemDimension> &rightData, FaceData<ProblemDimension> &edgeData, const Face&)
{

    auto edge_u_s = (leftData.getVelocitySolid() * edgeData.LeftCellKoef) + (rightData.getVelocitySolid() * edgeData.RightCellKoef);

    double product_of_u_s_and_n = (edge_u_s * edgeData.n);

    // prepare upwind values at the face/edge
    const FlowData<ProblemDimension>& faceVal = product_of_u_s_and_n > 0 ? leftData : rightData;

    // update the velocity at the face/edge
    // product_of_u_s_and_n = faceVal.getVelocityGas() * edgeData.n;


    // flux of mass
    double fluxRho_s = (-rho_s * faceVal.eps_s * edgeData.Length * product_of_u_s_and_n  +
                       (rightData.eps_s - leftData.eps_s) * rho_s * edgeData.LengthOverDist * artifitialDisspation) ;


    // computing the flux of momentum
    Vector<ProblemDimension,double> fluxP_s = -product_of_u_s_and_n * faceVal.p_s;


    // adding the element of "pressure" gradient
    fluxP_s -= G(leftData.getEps_g() * edgeData.LeftCellKoef + rightData.getEps_g() * edgeData.RightCellKoef) *
               (leftData.eps_s * edgeData.LeftCellKoef + rightData.eps_s * edgeData.RightCellKoef)* edgeData.n;

    fluxP_s *= edgeData.Length;

    // add artifitial dissipation
    fluxP_s += (edgeData.LengthOverDist * artifitialDisspation ) * (rightData.p_s - leftData.p_s);


    edgeData.grad_u_s = tensorProduct(edge_u_s, edgeData.n) * edgeData.Length;


    // storing result at the edge
    edgeData.fluxP_s = fluxP_s; // + viscose_diag + viscose_side_diag;

    edgeData.fluxRho_s = fluxRho_s;
}


inline void MultiphaseFlow::ComputeFluxSolid_inflow(const FlowData<ProblemDimension>& innerCellData, const Cell&,  FaceData<ProblemDimension>& edgeData, const Face& face)
{
    (void)innerCellData;

    Vector<ProblemDimension,double> modulatedU = inFlow_u_s;

    modulatedU *= FlowModulation(face.getCenter());

    double productOf_u_And_n = (modulatedU * edgeData.n);


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


    edgeData.grad_u_s = tensorProduct(modulatedU, edgeData.n) * edgeData.Length;


    edgeData.fluxP_s = fluxP_s;// + viscose_side_diag + viscose_diag;

    edgeData.fluxRho_s = fluxRho_s;

}

inline void MultiphaseFlow::ComputeFluxSolid_outflow(const FlowData<ProblemDimension>& innerCellData, FaceData<ProblemDimension>& edgeData, const Face&)
{

    Vector<ProblemDimension, double> in_u_s = innerCellData.getVelocitySolid();

    double productOf_u_And_n = ((in_u_s) * edgeData.n);
    double eps_s_e = productOf_u_And_n > 0 ? innerCellData.eps_s : outFlow.eps_s;




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

    edgeData.grad_u_s = tensorProduct(in_u_s, edgeData.n) * edgeData.Length;


    edgeData.fluxP_s = fluxP_s;// + viscose_side_diag + viscose_diag;

    edgeData.fluxRho_s = fluxRho_s;

}

inline void MultiphaseFlow::ComputeFluxSolid_wall(const FlowData<ProblemDimension>& innerCellData, FaceData<ProblemDimension>& edgeData, const Face&)
{



    // the flux of momentum
    // is reduced only to pressure gradient
    Vector<ProblemDimension,double> flux = -G(innerCellData.getEps_g()) * innerCellData.eps_s * edgeData.Length * edgeData.n;


    // the velocity is 0
     edgeData.grad_u_s = {};// tensorProduct({}, edgeData.n) * edgeData.Length;

    edgeData.fluxP_s = flux;// + viscose_side_diag + viscose_diag;

}



void MultiphaseFlow::ComputeViscousFluxSolid_inner(const FlowData<ProblemDimension> &leftData, const FlowData<ProblemDimension> &rightData, FaceData<ProblemDimension> &edgeData, const Face &fcData)
{

    double lEps_s = leftData.eps_s;
    double rEps_s = rightData.eps_s;

    double eps_s_e = (lEps_s * edgeData.LeftCellKoef + rEps_s * edgeData.RightCellKoef);

    auto edge_grad_u_s = (meshData.getDataByDim<ProblemDimension>()[fcData.getCellLeftIndex()].grad_u_s * edgeData.LeftCellKoef + meshData.getDataByDim<ProblemDimension>()[fcData.getCellRightIndex()].grad_u_s * edgeData.RightCellKoef);
    decltype (edge_grad_u_s) viscousT;

    for (unsigned int i = 0; i < edge_grad_u_s.size(); i++) {
        for (unsigned int j = 0; j < edge_grad_u_s.size(); j++) {
            // grad_u + transposed grad_u
            viscousT[i][j] = (edge_grad_u_s[i][j] + edge_grad_u_s[i][j]) * eps_s_e * myu;
        }
    }

    for (unsigned int i = 0; i < viscousT.size(); i++) {
        viscousT[i][i] *= 2.0 / 3.0;
    }

    edgeData.fluxP_g += tensorVectorProduct(viscousT, edgeData.n);

}

void MultiphaseFlow::ComputeViscousFluxSolid_inflow(const Cell &innerCell, const Face& face)
{

    double eps_s_e = inFlow_eps_s;

    auto edge_grad_u_s = (meshData[innerCell].grad_u_s);
    decltype (edge_grad_u_s) viscousT;

    for (unsigned int i = 0; i < edge_grad_u_s.size(); i++) {
        for (unsigned int j = 0; j < edge_grad_u_s.size(); j++) {
            // grad_u + transposed grad_u
            viscousT[i][j] = (edge_grad_u_s[i][j] + edge_grad_u_s[i][j]) * eps_s_e * myu;
        }
    }

    for (unsigned int i = 0; i < viscousT.size(); i++) {
        viscousT[i][i] *= 2.0 / 3.0;
    }

    meshData[face].fluxP_s += tensorVectorProduct(viscousT, meshData[face].n);

}

void MultiphaseFlow::ComputeViscousFluxSolid_outflow(const Cell& innerCell,  const Face& face)
{
    ComputeViscousFluxSolid_inflow(innerCell, face);
}

void MultiphaseFlow::ComputeViscousFluxSolid_wall(const Cell& innerCell,  const Face& face)
{
    ComputeViscousFluxSolid_inflow(innerCell, face);
}


void MultiphaseFlow::ComupteGradU(const Cell &cell)
{
    meshData.at(cell).grad_u_g = {};
    meshData.at(cell).grad_u_s = {};

    MeshApply<ProblemDimension, ProblemDimension - 1>::apply(
                cell.getIndex(),
                mesh,
                [&](size_t cellIndex, size_t faceIndex){
            const FaceData<ProblemDimension>& eData = meshData.getDataByDim<ProblemDimension - 1>().at(faceIndex);

            if (cellIndex == mesh.getFaces().at(faceIndex).getCellLeftIndex()){
                meshData.at(cell).grad_u_g += eData.grad_u_g;
                meshData.at(cell).grad_u_s += eData.grad_u_s;
            } else {
                meshData.at(cell).grad_u_g -= eData.grad_u_g;
                meshData.at(cell).grad_u_s -= eData.grad_u_s;
            }
        }
                );
    meshData.at(cell).grad_u_g *= meshData.at(cell).invVolume;
    meshData.at(cell).grad_u_s *= meshData.at(cell).invVolume;
}







#endif // MULTIPHASEFLOW_H

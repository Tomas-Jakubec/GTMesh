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


    static double R_spec;
    static double T;
    static double rho_s;


    // ** Primary quantities **
    /**
     * @brief Density of the gaseous part
     */
    double rho_g_x_eps_g;
    /**
     * @brief p_g
     *
     * momentum of gaseous phase
     */
    Vector<Dim,double> p_g;

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

    // ** Dependent quantities **
    /**
     * @brief eps_g
     * volume fraction of gaseous part of the flow
     * holds \epsilon_{s} + \epsilon_{g} = 1
     */
    //double eps_g;
    double getEps_g() const {return 1.0 - eps_s;}
    void setEps_g(const double& eps_g){eps_s = 1.0 - eps_g;}

    /**
     * @brief
     * density of the gaseous phase
     */
    double getRho_g() const {return rho_g_x_eps_g / reg(getEps_g());}
    void setRho_g(const double& rho_g){rho_g_x_eps_g = getEps_g() * rho_g;}


    /**
     * @brief u_g
     *
     * velocity of gaseous phase
     */
    Vector<Dim, double> getVelocityGas() const {
        return p_g / reg(rho_g_x_eps_g);
    }
    void setVelocityGas(const Vector<Dim, double>& u_g){
        p_g = u_g * rho_g_x_eps_g;
    }


    // gasseous part
    // pressure is dependent on the R_spec and T and rho_g
    //double p;
    double getPressure() const { return getRho_g() * R_spec * T; }
    void setPressure(const double& pressure){ setRho_g( pressure / ( R_spec * T ) ); }


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



};

template<unsigned int Dim>
double FlowData<Dim>::rho_s = 0;
template<unsigned int Dim>
double FlowData<Dim>::R_spec = 0;
template<unsigned int Dim>
double FlowData<Dim>::T = 0;

MAKE_ATTRIBUTE_TRAIT_ARITHMETIC(FlowData<2>, rho_g_x_eps_g, eps_s, p_g, p_s);

// TODO only temporary
//MAKE_ATTRIBUTE_TRAIT_IO(FlowData, rho_g, eps_s, p_g, p_s)


MAKE_CUSTOM_TRAIT( FlowData<2>,
                   "eps_g", std::make_pair(&FlowData<2>::getEps_g, &FlowData<2>::setEps_g),
                   "pressure", std::make_pair(&FlowData<2>::getPressure, &FlowData<2>::setPressure),
                   "rho_g", std::make_pair(&FlowData<2>::getRho_g, &FlowData<2>::setRho_g),
                   "eps_s", &FlowData<2>::eps_s,
                   "velocity_gas", std::make_pair(&FlowData<2>::getVelocityGas, &FlowData<2>::setVelocityGas),
                   "velocity_solid", std::make_pair(&FlowData<2>::getVelocitySolid, &FlowData<2>::setVelocitySolid) );

MAKE_ATTRIBUTE_TRAIT_ARITHMETIC(FlowData<3>, rho_g_x_eps_g, eps_s, p_g, p_s);

// TODO only temporary
//MAKE_ATTRIBUTE_TRAIT_IO(FlowData, rho_g, eps_s, p_g, p_s)


MAKE_CUSTOM_TRAIT( FlowData<3>,
                   "eps_g", std::make_pair(&FlowData<3>::getEps_g, &FlowData<3>::setEps_g),
                   "pressure", std::make_pair(&FlowData<3>::getPressure, &FlowData<3>::setPressure),
                   "rho_g", std::make_pair(&FlowData<3>::getRho_g, &FlowData<3>::setRho_g),
                   "eps_s", &FlowData<3>::eps_s,
                   "velocity_gas", std::make_pair(&FlowData<3>::getVelocityGas, &FlowData<3>::setVelocityGas),
                   "velocity_solid", std::make_pair(&FlowData<3>::getVelocitySolid, &FlowData<3>::setVelocitySolid) );



template<unsigned int Dim>
struct FaceData {
    double MeasureOverDist;
    double Measure;
    Vector<Dim, double> n;

    // koeficients of konvex combination of cells
    // holds that LeftCellKoef + RightCellKoef = 1
    double LeftCellKoef;
    double RightCellKoef;

    // // koeficients of konvex combination for normals on tessellated mesh
    // // holds that LeftCellNormal + RightCellNormal = n
    // Vector<Dim, double> LeftCellNormal;
    // Vector<Dim, double> RightCellNormal;




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

MAKE_ATTRIBUTE_TRAIT(FaceData<2>, fluxRho_g,fluxP_g,RightCellKoef,LeftCellKoef,n,Measure,MeasureOverDist);
MAKE_ATTRIBUTE_TRAIT(FaceData<3>, fluxRho_g,fluxP_g,RightCellKoef,LeftCellKoef,n,Measure,MeasureOverDist);
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





/**
 * @brief Calculates the problem of two-phase flow.
 */
template< unsigned int Dimension,
          typename BoundaryCond,
          unsigned int ... Reserve >
class MultiphaseFlow {
public:
    MultiphaseFlow() = default;
    ~MultiphaseFlow() = default;

public: // make some types public
    static constexpr unsigned int ProblemDimension = Dimension;
    using MeshType = UnstructuredMesh<ProblemDimension,size_t, double, Reserve...>;
    using ResultType = FlowData<ProblemDimension>;

    using Cell = typename MeshType::Cell;
    using Face = typename MeshType::Face;
    using Vert = typename MeshType::Vertex;
public:
    // Structures containing state data

    MeshType mesh;

    MeshDataContainer<std::tuple<CellData<ProblemDimension>,FaceData<ProblemDimension>>, ProblemDimension,ProblemDimension - 1> meshData;
private:

    VTKMeshWriter<ProblemDimension> writer;
    std::unique_ptr<MeshReader<ProblemDimension>> reader;

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

    template<unsigned int _Dimension = Dimension>
    typename std::enable_if<_Dimension == 3>::type
    setupMeshData(const std::string& fileName);

    template<unsigned int _Dimension = Dimension>
    typename std::enable_if<_Dimension == 2>::type
    setupMeshData(const std::string& fileName);
//
//

    void calculateRHS(
             double time,//time is unused in this problem
             MeshDataContainer<ResultType, ProblemDimension>& compData,
             MeshDataContainer<ResultType, ProblemDimension>& outDeltas);


    Type TypeOfCell(const Cell& cell);


    void exportData( double time,
                     MeshDataContainer<ResultType, ProblemDimension>& compData,
                     double timeModifier = 1e2);

};

MAKE_ATTRIBUTE_TEMPLATE_TRAIT(PASS(MultiphaseFlow<Dim, BC, Res...>),
                              PASS(unsigned int Dim, typename BC, unsigned int... Res),
                              myu, myu_s, R_spec, T, artifitialDisspation, d_s, phi_s, rho_s,
                              outFlow, inFlow_u_s, inFlow_u_g, inFlow_eps_s, inFlow_eps_g);


/*
 * Definition of inlined functions
 *
*/
template<unsigned int Dimension, typename BoundaryCond, unsigned int ... Reserve>
inline double MultiphaseFlow< Dimension, BoundaryCond, Reserve... >::G(double eps_g)
{
    return (std::pow(10, -8.76 * eps_g + 5.43));
}





template<unsigned int Dimension, typename BoundaryCond, unsigned int ... Reserve>
inline double MultiphaseFlow< Dimension, BoundaryCond, Reserve... >::Beta_s(const FlowData<ProblemDimension> &fd)
{
    if(fd.eps_s > 0.2) {
        double denominator = 1 / (fd.getEps_g() * d_s * phi_s);
        return (
                    150 * (std::pow(fd.eps_s, 2.0) * myu) * std::pow(denominator, 2.0) +
                    1.75 * (fd.getVelocityGas() - fd.getVelocitySolid()).normEuclid() * fd.getRho_g() * fd.eps_s * denominator
                );
    } else {
        return (
                    (4.0/3.0) * C_d(fd) * (fd.getVelocityGas() - fd.getVelocitySolid()).normEuclid() * fd.getRho_g() * fd.eps_s / (d_s * phi_s)
               );
    }
}





template<unsigned int Dimension, typename BoundaryCond, unsigned int ... Reserve>
inline double MultiphaseFlow< Dimension, BoundaryCond, Reserve... >::C_d(const FlowData<ProblemDimension> &fd)
{
    double re_s = reg(Re_s(fd));

    if (re_s <= 1000) {
        return (24/re_s)*(1.0 + 0.15 * std::pow(re_s, 0.687));
    } else {
        return 0.44;
    }
}





template<unsigned int Dimension, typename BoundaryCond, unsigned int ... Reserve>
inline double MultiphaseFlow< Dimension, BoundaryCond, Reserve... >::Re_s(const FlowData<ProblemDimension> &fd)
{
    // multiplying by inverted value of myu may make code faster
    return(
            (fd.getVelocityGas() - fd.getVelocitySolid()).normEuclid() * d_s * fd.rho_g_x_eps_g /  myu
          );
}



template<unsigned int Dimension, typename BoundaryCond, unsigned int ... Reserve>
inline void MultiphaseFlow< Dimension, BoundaryCond, Reserve... >::ComputeFluxGas_inner(const FlowData<ProblemDimension> &leftData, const FlowData<ProblemDimension> &rightData, FaceData<ProblemDimension> &edgeData, const Face&)
{

    auto edge_u_g = (leftData.getVelocityGas() * edgeData.LeftCellKoef) + (rightData.getVelocityGas() * edgeData.RightCellKoef);
    double product_of_u_and_n = (edge_u_g * edgeData.n);

    // prepare upwind values
    const FlowData<ProblemDimension>& faceVal = product_of_u_and_n > 0 ? leftData : rightData;

    // update the velocity at the face/edge
    //product_of_u_and_n = faceVal.getVelocityGas() * edgeData.n;

    // flux of density
    double delta_rho = - faceVal.rho_g_x_eps_g * product_of_u_and_n * edgeData.Measure +
                       (rightData.rho_g_x_eps_g - leftData.rho_g_x_eps_g) * edgeData.MeasureOverDist * artifitialDisspation;


    // computing the flux of momentum
    Vector<ProblemDimension,double> fluxP_g = (-product_of_u_and_n) * faceVal.p_g;

    // adding the element of pressure gradient
    fluxP_g -= (leftData.getPressure() * edgeData.LeftCellKoef + rightData.getPressure() * edgeData.RightCellKoef) *
               (leftData.getEps_g() * edgeData.LeftCellKoef + rightData.getEps_g() * edgeData.RightCellKoef) * edgeData.n;

    // multiply by the face/edge measure
    fluxP_g *= edgeData.Measure;

    // add artifitial dissipation
    fluxP_g += (edgeData.MeasureOverDist * artifitialDisspation * ((rightData.p_g) - (leftData.p_g)));


    // compuatation of grad(u)
    edgeData.grad_u_g = tensorProduct(edge_u_g, edgeData.n) * edgeData.Measure;

    // adding of fluxes to cells
    edgeData.fluxP_g = fluxP_g;

    edgeData.fluxRho_g = delta_rho;


}


template<unsigned int Dimension, typename BoundaryCond, unsigned int ... Reserve>
inline void MultiphaseFlow< Dimension, BoundaryCond, Reserve... >::ComputeFluxGas_inflow(const FlowData<ProblemDimension>& innerCellData, const Cell& ,  FaceData<ProblemDimension>& edgeData, const Face& face)
{

    Vector<ProblemDimension,double> modulatedU = inFlow_u_g;

    modulatedU *= FlowModulation(face.getCenter());

    double product_of_u_and_n = (modulatedU * edgeData.n);


    // flux of density
    double delta_rho = -innerCellData.getRho_g() * inFlow_eps_g *
                        product_of_u_and_n * edgeData.Measure;




    // computing the flux of momentum
    Vector<ProblemDimension,double> flux = - (innerCellData.getRho_g()) *
                    (product_of_u_and_n) * inFlow_eps_g *
                    (modulatedU);


    // adding the element of pressure gradient
    flux -= (innerCellData.getPressure()) * edgeData.n;


    flux *= edgeData.Measure;



    // compuatation of grad(u)
    edgeData.grad_u_g = tensorProduct(modulatedU, edgeData.n) * edgeData.Measure;

    edgeData.fluxP_g = flux;// + viscose_side_diag + viscose_diag;

    edgeData.fluxRho_g = delta_rho;

}



template<unsigned int Dimension, typename BoundaryCond, unsigned int ... Reserve>
inline void MultiphaseFlow< Dimension, BoundaryCond, Reserve... >::ComputeFluxGas_outflow(const FlowData<ProblemDimension>& innerCellData, FaceData<ProblemDimension>& edgeData, const Face&)
{
    Vector<ProblemDimension, double> in_u_g = innerCellData.getVelocityGas();
    double product_of_u_and_n = ((in_u_g) * edgeData.n);
    double eps_g_e = product_of_u_and_n > 0 ? innerCellData.getEps_g() : outFlow.getEps_g();



    // flux of density
    double delta_rho = -outFlow.getRho_g() * eps_g_e *
                        product_of_u_and_n * edgeData.Measure;


    // computing the flux of momentum
    Vector<ProblemDimension,double> flux = - (outFlow.getRho_g() *
                    product_of_u_and_n * eps_g_e) *
                    (in_u_g);


    // adding the element of pressure gradient
    flux -= (outFlow.getPressure()) * edgeData.n;

    flux *= edgeData.Measure;


    // compuatation of grad(u)
    edgeData.grad_u_g = tensorProduct(in_u_g, edgeData.n) * edgeData.Measure;

    edgeData.fluxP_g = flux;// + viscose_side_diag + viscose_diag;

    edgeData.fluxRho_g = delta_rho;

}


template<unsigned int Dimension, typename BoundaryCond, unsigned int ... Reserve>
inline void MultiphaseFlow< Dimension, BoundaryCond, Reserve... >::ComputeFluxGas_wall(const FlowData<ProblemDimension>& innerCellData, FaceData<ProblemDimension>& edgeData, const Face&)
{





    // the flux of momentum
    // is reduced only to pressure gradient
    Vector<ProblemDimension,double> flux = -(innerCellData.getPressure()) * edgeData.n;


    flux *= edgeData.Measure;

    //prepare eps_g

    // compuatation of grad(u)
    edgeData.grad_u_g = {};

    edgeData.fluxP_g = flux;// + viscose_side_diag + viscose_diag;
    edgeData.fluxRho_g = 0;
}






template<unsigned int Dimension, typename BoundaryCond, unsigned int ... Reserve>
void MultiphaseFlow< Dimension, BoundaryCond, Reserve... >::ComputeViscousFluxGas_inner(const FlowData<ProblemDimension> &leftData, const FlowData<ProblemDimension> &rightData, FaceData<ProblemDimension> &edgeData, const Face &fcData)
{

    double lEps_g = leftData.getEps_g();
    double rEps_g = rightData.getEps_g();

    double eps_g_e = (lEps_g * edgeData.LeftCellKoef + rEps_g * edgeData.RightCellKoef);

    auto edge_grad_u_g = (meshData.template getDataByDim<ProblemDimension>()[fcData.getCellLeftIndex()].grad_u_g * edgeData.LeftCellKoef + meshData. template getDataByDim<ProblemDimension>()[fcData.getCellRightIndex()].grad_u_g * edgeData.RightCellKoef);
    decltype (edge_grad_u_g) viscousT;

    for (unsigned int i = 0; i < edge_grad_u_g.size(); i++) {
        for (unsigned int j = 0; j < edge_grad_u_g.size(); j++) {
            // grad_u + transposed grad_u
            viscousT[i][j] = (edge_grad_u_g[i][j] + edge_grad_u_g[i][j]) * eps_g_e * myu;
        }
    }

    double div_u = 0;
    for (unsigned int i = 0; i < viscousT.size(); i++) {
        div_u += edge_grad_u_g[i][i];
    }


    for (unsigned int i = 0; i < viscousT.size(); i++) {
        viscousT[i][i] -= (2.0 / 3.0) * eps_g_e * myu * div_u;
    }


    edgeData.fluxP_g += tensorVectorProduct(viscousT, edgeData.n);

}


template<unsigned int Dimension, typename BoundaryCond, unsigned int ... Reserve>
void MultiphaseFlow< Dimension, BoundaryCond, Reserve... >::ComputeViscousFluxGas_inflow(const Cell &innerCell, const Face& face)
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

    double div_u = 0;
    for (unsigned int i = 0; i < viscousT.size(); i++) {
        div_u += edge_grad_u_g[i][i];
    }


    for (unsigned int i = 0; i < viscousT.size(); i++) {
        viscousT[i][i] -= (2.0 / 3.0) * eps_g_e * myu * div_u;
    }


    meshData[face].fluxP_g += tensorVectorProduct(viscousT, meshData[face].n);

}

template<unsigned int Dimension, typename BoundaryCond, unsigned int ... Reserve>
void MultiphaseFlow< Dimension, BoundaryCond, Reserve... >::ComputeViscousFluxGas_outflow(const Cell& innerCell,  const Face& face)
{
    ComputeViscousFluxGas_inflow(innerCell, face);
}

template<unsigned int Dimension, typename BoundaryCond, unsigned int ... Reserve>
void MultiphaseFlow< Dimension, BoundaryCond, Reserve... >::ComputeViscousFluxGas_wall(const Cell& innerCell,  const Face& face)
{
    ComputeViscousFluxGas_inflow(innerCell, face);
}







template<unsigned int Dimension, typename BoundaryCond, unsigned int ... Reserve>
inline void MultiphaseFlow< Dimension, BoundaryCond, Reserve... >::ComputeFluxSolid_inner(const FlowData<ProblemDimension> &leftData, const FlowData<ProblemDimension> &rightData, FaceData<ProblemDimension> &edgeData, const Face&)
{

    auto edge_u_s = (leftData.getVelocitySolid() * edgeData.LeftCellKoef) + (rightData.getVelocitySolid() * edgeData.RightCellKoef);

    double product_of_u_s_and_n = (edge_u_s * edgeData.n);

    // prepare upwind values at the face/edge
    const FlowData<ProblemDimension>& faceVal = product_of_u_s_and_n > 0 ? leftData : rightData;

    // flux of mass
    double fluxRho_s = (-rho_s * faceVal.eps_s * edgeData.Measure * product_of_u_s_and_n  +
                       (rightData.eps_s - leftData.eps_s) * rho_s * edgeData.MeasureOverDist * artifitialDisspation) ;


    // computing the flux of momentum
    Vector<ProblemDimension,double> fluxP_s = -product_of_u_s_and_n * faceVal.p_s;


    // adding the element of "pressure" gradient
    // The gaseous pressure does not bring better results.
    //fluxP_s -= (leftData.getPressure() * edgeData.LeftCellKoef + rightData.getPressure() * edgeData.RightCellKoef) *
    //               (leftData.eps_s * edgeData.LeftCellKoef + rightData.eps_s * edgeData.RightCellKoef) * edgeData.n;

    fluxP_s -= G(leftData.getEps_g() * edgeData.LeftCellKoef + rightData.getEps_g() * edgeData.RightCellKoef) *
               ((leftData.eps_s * edgeData.LeftCellKoef) + (rightData.eps_s * edgeData.RightCellKoef))* edgeData.n;


    fluxP_s *= edgeData.Measure;

    // add artifitial dissipation
    fluxP_s += (edgeData.MeasureOverDist * artifitialDisspation ) * (rightData.p_s - leftData.p_s);


    edgeData.grad_u_s = tensorProduct(edge_u_s, edgeData.n) * edgeData.Measure;


    // storing result at the edge
    edgeData.fluxP_s = fluxP_s;

    edgeData.fluxRho_s = fluxRho_s;
}

template<unsigned int Dimension, typename BoundaryCond, unsigned int ... Reserve>
inline void MultiphaseFlow< Dimension, BoundaryCond, Reserve... >::ComputeFluxSolid_inflow(const FlowData<ProblemDimension>& innerCellData, const Cell&,  FaceData<ProblemDimension>& edgeData, const Face& face)
{
    (void)innerCellData;

    Vector<ProblemDimension,double> modulatedU = inFlow_u_s;

    modulatedU *= FlowModulation(face.getCenter());

    double productOf_u_And_n = (modulatedU * edgeData.n);


    // flux of density
    double fluxRho_s = -rho_s * inFlow_eps_s *
                        productOf_u_And_n * edgeData.Measure;




    // computing the flux of momentum
    Vector<ProblemDimension,double> fluxP_s = - (rho_s) *
                    (productOf_u_And_n) * inFlow_eps_s *
                    (modulatedU);

    // adding the element of pressure gradient
    fluxP_s -= -G(inFlow_eps_g) * inFlow_eps_s * edgeData.n;


    fluxP_s *= edgeData.Measure;


    edgeData.grad_u_s = tensorProduct(modulatedU, edgeData.n) * edgeData.Measure;


    edgeData.fluxP_s = fluxP_s;// + viscose_side_diag + viscose_diag;

    edgeData.fluxRho_s = fluxRho_s;

}


template<unsigned int Dimension, typename BoundaryCond, unsigned int ... Reserve>
inline void MultiphaseFlow< Dimension, BoundaryCond, Reserve... >::ComputeFluxSolid_outflow(const FlowData<ProblemDimension>& innerCellData, FaceData<ProblemDimension>& edgeData, const Face&)
{

    Vector<ProblemDimension, double> in_u_s = innerCellData.getVelocitySolid();

    double productOf_u_And_n = ((in_u_s) * edgeData.n);
    double eps_s_e = productOf_u_And_n > 0 ? innerCellData.eps_s : outFlow.eps_s;




    // flux of density
    double fluxRho_s = -rho_s * eps_s_e *
                        productOf_u_And_n * edgeData.Measure;


    // computing the flux of momentum
    Vector<ProblemDimension,double> fluxP_s = - (rho_s) *
                    (productOf_u_And_n) * eps_s_e *
                    (in_u_s);


    // adding the element of pressure gradient
    fluxP_s -= -G(1.0 - eps_s_e) * eps_s_e * edgeData.n;


    fluxP_s *= edgeData.Measure;

    edgeData.grad_u_s = tensorProduct(in_u_s, edgeData.n) * edgeData.Measure;


    edgeData.fluxP_s = fluxP_s;// + viscose_side_diag + viscose_diag;

    edgeData.fluxRho_s = fluxRho_s;

}


template<unsigned int Dimension, typename BoundaryCond, unsigned int ... Reserve>
inline void MultiphaseFlow< Dimension, BoundaryCond, Reserve... >::ComputeFluxSolid_wall(const FlowData<ProblemDimension>& innerCellData, FaceData<ProblemDimension>& edgeData, const Face&)
{



    // the flux of momentum
    // is reduced only to pressure gradient
    Vector<ProblemDimension,double> flux = -G(innerCellData.getEps_g()) * innerCellData.eps_s * edgeData.Measure * edgeData.n;


    // the velocity is 0
     edgeData.grad_u_s = {};// tensorProduct({}, edgeData.n) * edgeData.Length;

    edgeData.fluxP_s = flux;// + viscose_side_diag + viscose_diag;

}


template<unsigned int Dimension, typename BoundaryCond, unsigned int ... Reserve>
void MultiphaseFlow< Dimension, BoundaryCond, Reserve... >::ComputeViscousFluxSolid_inner(const FlowData<ProblemDimension> &leftData, const FlowData<ProblemDimension> &rightData, FaceData<ProblemDimension> &edgeData, const Face &fcData)
{

    double lEps_s = leftData.eps_s;
    double rEps_s = rightData.eps_s;

    double eps_s_e = (lEps_s * edgeData.LeftCellKoef + rEps_s * edgeData.RightCellKoef);

    auto edge_grad_u_s = (meshData.template getDataByDim<ProblemDimension>()[fcData.getCellLeftIndex()].grad_u_s * edgeData.LeftCellKoef + meshData.template getDataByDim<ProblemDimension>()[fcData.getCellRightIndex()].grad_u_s * edgeData.RightCellKoef);
    decltype (edge_grad_u_s) viscousT;

    for (unsigned int i = 0; i < edge_grad_u_s.size(); i++) {
        for (unsigned int j = 0; j < edge_grad_u_s.size(); j++) {
            // grad_u + transposed grad_u
            viscousT[i][j] = (edge_grad_u_s[i][j] + edge_grad_u_s[i][j]) * eps_s_e * myu_s;
        }
    }

    double div_u = 0;
    for (unsigned int i = 0; i < viscousT.size(); i++) {
        div_u += edge_grad_u_s[i][i];
    }


    for (unsigned int i = 0; i < viscousT.size(); i++) {
        viscousT[i][i] -= (2.0 / 3.0) * eps_s_e * myu_s * div_u;
    }

    edgeData.fluxP_g += tensorVectorProduct(viscousT, edgeData.n);

}


template<unsigned int Dimension, typename BoundaryCond, unsigned int ... Reserve>
void MultiphaseFlow< Dimension, BoundaryCond, Reserve... >::ComputeViscousFluxSolid_inflow(const Cell &innerCell, const Face& face)
{

    double eps_s_e = inFlow_eps_s;

    auto edge_grad_u_s = (meshData[innerCell].grad_u_s);
    decltype (edge_grad_u_s) viscousT;

    for (unsigned int i = 0; i < edge_grad_u_s.size(); i++) {
        for (unsigned int j = 0; j < edge_grad_u_s.size(); j++) {
            // grad_u + transposed grad_u
            viscousT[i][j] = (edge_grad_u_s[i][j] + edge_grad_u_s[i][j]) * eps_s_e * myu_s;
        }
    }

    double div_u = 0;
    for (unsigned int i = 0; i < viscousT.size(); i++) {
        div_u += edge_grad_u_s[i][i];
    }


    for (unsigned int i = 0; i < viscousT.size(); i++) {
        viscousT[i][i] -= (2.0 / 3.0) * eps_s_e * myu_s * div_u;
    }

    meshData[face].fluxP_s += tensorVectorProduct(viscousT, meshData[face].n);

}


template<unsigned int Dimension, typename BoundaryCond, unsigned int ... Reserve>
void MultiphaseFlow< Dimension, BoundaryCond, Reserve... >::ComputeViscousFluxSolid_outflow(const Cell& innerCell,  const Face& face)
{
    ComputeViscousFluxSolid_inflow(innerCell, face);
}


template<unsigned int Dimension, typename BoundaryCond, unsigned int ... Reserve>
void MultiphaseFlow< Dimension, BoundaryCond, Reserve... >::ComputeViscousFluxSolid_wall(const Cell& innerCell,  const Face& face)
{
    ComputeViscousFluxSolid_inflow(innerCell, face);
}

template<unsigned int Dimension, typename BoundaryCond, unsigned int ... Reserve>
void MultiphaseFlow< Dimension, BoundaryCond, Reserve... >::ComupteGradU(const Cell &cell)
{
    meshData.at(cell).grad_u_g = {};
    meshData.at(cell).grad_u_s = {};

    MeshApply<ProblemDimension, ProblemDimension - 1>::apply(
                cell.getIndex(),
                mesh,
                [&](size_t cellIndex, size_t faceIndex){
            const FaceData<ProblemDimension>& eData = meshData.template getDataByDim<ProblemDimension - 1>().at(faceIndex);

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



// From multiphaseflow.cpp
// in order to make dimension as template parameter


template<unsigned int Dimension, typename BoundaryCond, unsigned int ... Reserve>
void MultiphaseFlow< Dimension, BoundaryCond, Reserve... >::calculateRHS(double, MeshDataContainer<MultiphaseFlow< Dimension, BoundaryCond, Reserve... >::ResultType, MultiphaseFlow< Dimension, BoundaryCond, Reserve... >::ProblemDimension> &compData, MeshDataContainer<MultiphaseFlow< Dimension, BoundaryCond, Reserve... >::ResultType, MultiphaseFlow< Dimension, BoundaryCond, Reserve... >::ProblemDimension> &outDeltas)
{
    // in single compute step
    // we have to:

//#pragma omp parallel
{

    // first: correct negative volume fractions
    #pragma omp for
    for (size_t cellIndex = 0; cellIndex < mesh.getCells().size(); cellIndex++){
        auto& cell = mesh.getCells()[cellIndex];
        FlowData<ProblemDimension>& cellData = compData.at(cell);
        if (cellData.eps_s <= 0.0) {
             cellData.eps_s = 0.0;
             cellData.p_s = {};
         }
    }

    //UpdateVertexData(compData);

    // second: compute fluxes over edges

    #pragma omp for
    for (size_t faceIndex = 0; faceIndex < mesh.getFaces().size(); faceIndex++){
        auto& face = mesh.getFaces()[faceIndex];

        ComputeFlux(face, compData);

    }

    // third: compute approximations of gradients

    #pragma omp for
    for (size_t cellIndex = 0; cellIndex < mesh.getCells().size(); cellIndex++){
        auto& cell = mesh.getCells()[cellIndex];

        ComupteGradU(cell);
    }

    // fourth: compute source terms
    #pragma omp for
    for (size_t faceIndex = 0; faceIndex < mesh.getFaces().size(); faceIndex++){
        auto& face = mesh.getFaces()[faceIndex];

        ComputeViscousFlux(face, compData);

    }


    // fifth: compute source terms

    #pragma omp for
    for (size_t cellIndex = 0; cellIndex < mesh.getCells().size(); cellIndex++){
        auto& cell = mesh.getCells()[cellIndex];

        ComputeSource(cell, compData, outDeltas);
    }

}

}





template<unsigned int Dimension, typename BoundaryCond, unsigned int ... Reserve>
void MultiphaseFlow< Dimension, BoundaryCond, Reserve... >::ComputeFlux(const Face &fcData, const MeshDataContainer<ResultType, ProblemDimension>& compData)
{
    // first compute all variables interpolated
    // in the center of the edge
    if (fcData.getCellRightIndex() < BOUNDARY_INDEX(size_t) && fcData.getCellLeftIndex() < BOUNDARY_INDEX(size_t)) {

        const FlowData<ProblemDimension> &leftData = compData.template getDataByDim<ProblemDimension>().at(fcData.getCellLeftIndex());

        const FlowData<ProblemDimension> &rightData = compData.template getDataByDim<ProblemDimension>().at(fcData.getCellRightIndex());

        FaceData<ProblemDimension> &currFaceData = meshData.at(fcData);

        // Computation of fluxes of density and momentum of gaseous part
        ComputeFluxGas_inner(leftData, rightData, currFaceData, fcData);

        // Compute flux of solid phase
        ComputeFluxSolid_inner(leftData, rightData, currFaceData, fcData);

    } else {
        // Applying boundary conditions
        const Cell* innerCell = nullptr;
        const Cell* outerCell = nullptr;
        if (fcData.getCellRightIndex() < BOUNDARY_INDEX(size_t)){
            innerCell = &mesh.getCells().at(fcData.getCellRightIndex());
            outerCell = &mesh.getBoundaryCells().at(fcData.getCellLeftIndex() & EXTRACTING_INDEX(size_t));
        } else {
            innerCell = &mesh.getCells().at(fcData.getCellLeftIndex());
            outerCell = &mesh.getBoundaryCells().at(fcData.getCellRightIndex() & EXTRACTING_INDEX(size_t));
        }

        const FlowData<ProblemDimension> *innerCellData = &compData.at(*innerCell);

        FaceData<ProblemDimension> *currFaceData = &meshData.at(fcData);

        switch(outerCell->getFlag()){
        case INFLOW : {


            ComputeFluxGas_inflow(*innerCellData, *innerCell, *currFaceData, fcData);
            //ComputeFluxSolid_inflow(*innerCellData, *innerCell, *currFaceData, fcData);
            ComputeFluxSolid_wall(*innerCellData, *currFaceData, fcData);
        } break;

        case OUTFLOW :{

            ComputeFluxGas_outflow(*innerCellData, *currFaceData, fcData);
            ComputeFluxSolid_outflow(*innerCellData, *currFaceData, fcData);
            //ComputeFluxSolid_wall(*innerCellData, *currFaceData, fcData);
        } break;

        case WALL : {


                ComputeFluxGas_wall(*innerCellData, *currFaceData, fcData);
                ComputeFluxSolid_wall(*innerCellData, *currFaceData, fcData);

            } break;

        default : {
            throw std::runtime_error("unknown cell flag: " + std::to_string(outerCell->getFlag()) +
                                     " in cell number: " + std::to_string(outerCell->getIndex()));
            }
        }
    }

}


template<unsigned int Dimension, typename BoundaryCond, unsigned int ... Reserve>
void MultiphaseFlow< Dimension, BoundaryCond, Reserve... >::ComputeViscousFlux(const Face &fcData, const MeshDataContainer<ResultType, ProblemDimension> &compData)
{
    // first compute all variables interpolated
    // in the center of the edge
    if (fcData.getCellRightIndex() < BOUNDARY_INDEX(size_t) && fcData.getCellLeftIndex() < BOUNDARY_INDEX(size_t)) {

        const FlowData<ProblemDimension> &leftData = compData.template getDataByDim<ProblemDimension>().at(fcData.getCellLeftIndex());

        const FlowData<ProblemDimension> &rightData = compData.template getDataByDim<ProblemDimension>().at(fcData.getCellRightIndex());

        FaceData<ProblemDimension> &currFaceData = meshData.at(fcData);

        // Computation of fluxes of density and momentum of gaseous part
        ComputeViscousFluxGas_inner(leftData, rightData, currFaceData, fcData);

        // Compute flux of solid phase
        ComputeViscousFluxSolid_inner(leftData, rightData, currFaceData, fcData);

    } else {
        // Applying boundary conditions
        const Cell* innerCell = nullptr;
        const Cell* outerCell = nullptr;
        if (fcData.getCellRightIndex() < BOUNDARY_INDEX(size_t)){
            innerCell = &mesh.getCells().at(fcData.getCellRightIndex());
            outerCell = &mesh.getBoundaryCells().at(fcData.getCellLeftIndex() & EXTRACTING_INDEX(size_t));
        } else {
            innerCell = &mesh.getCells().at(fcData.getCellLeftIndex());
            outerCell = &mesh.getBoundaryCells().at(fcData.getCellRightIndex() & EXTRACTING_INDEX(size_t));
        }


        switch(outerCell->getFlag()){
        case INFLOW : {


            ComputeViscousFluxGas_inflow(*innerCell, fcData);
            ComputeViscousFluxSolid_inflow(*innerCell, fcData);
            //ComputeFluxSolid_wall(*innerCellData, *currFaceData, fcData);
        } break;

        case OUTFLOW :{

            ComputeViscousFluxGas_outflow(*innerCell, fcData);
            ComputeViscousFluxSolid_outflow(*innerCell, fcData);
            //ComputeFluxSolid_wall(*innerCellData, *currFaceData, fcData);
        } break;

        case WALL : {


                ComputeViscousFluxGas_wall(*innerCell, fcData);
                ComputeViscousFluxSolid_wall(*innerCell, fcData);

            } break;

        default : {
            throw std::runtime_error("unknown cell flag: " + std::to_string(outerCell->getFlag()) +
                                     " in cell number: " + std::to_string(outerCell->getIndex()));
            }
        }
    }

}






template<unsigned int Dimension, typename BoundaryCond, unsigned int ... Reserve>
void MultiphaseFlow< Dimension, BoundaryCond, Reserve... >::ComputeSource(const Cell& cell,
                                   MeshDataContainer<ResultType, ProblemDimension>& compData,
                                   MeshDataContainer<ResultType, ProblemDimension> &result)
{
    FlowData<ProblemDimension>& resData = result[cell];
    FlowData<ProblemDimension>& cellData = compData[cell];


    resData.rho_g_x_eps_g = 0;
    resData.eps_s = 0; // Firstly, the flux of rho_s is calculated. Finally, the eps_s is computed as flux_rho_s / rho_s
    resData.p_g = {};
    resData.p_s = {};

    MeshApply<ProblemDimension, ProblemDimension - 1>::apply(
                cell.getIndex(),
                mesh,
                // aplication of sum lambda to all cell faces
                [&](size_t cellIndex, size_t faceIndex){
            const FaceData<ProblemDimension>& eData = meshData.template getDataByDim<ProblemDimension - 1>().at(faceIndex);

            if (cellIndex == mesh.getFaces().at(faceIndex).getCellLeftIndex()){
                resData.rho_g_x_eps_g += eData.fluxRho_g;
                resData.eps_s         += eData.fluxRho_s;
                resData.p_g           += eData.fluxP_g;
                resData.p_s           += eData.fluxP_s;
            } else {
                resData.rho_g_x_eps_g -= eData.fluxRho_g;
                resData.eps_s         -= eData.fluxRho_s;
                resData.p_g           -= eData.fluxP_g;
                resData.p_s           -= eData.fluxP_s;
            }
        }
    );

    resData.rho_g_x_eps_g *= meshData[cell].invVolume;
    resData.eps_s         *= meshData[cell].invVolume;
    resData.p_g           *= meshData[cell].invVolume;
    resData.p_s           *= meshData[cell].invVolume;



    Vector<ProblemDimension, double> drag = Beta_s(cellData) * (cellData.getVelocitySolid() - cellData.getVelocityGas());

    Vector<ProblemDimension,double> g_acceleration = {};

    g_acceleration[ProblemDimension - 2] = -9.81;

    resData.p_g += (cellData.getRho_g() * g_acceleration + drag);

    resData.p_s += ((rho_s - cellData.getRho_g()) * cellData.eps_s * g_acceleration - drag);


    resData.eps_s /= rho_s;

}







template<unsigned int Dimension, typename BoundaryCond, unsigned int ... Reserve>
double MultiphaseFlow< Dimension, BoundaryCond, Reserve... >::FlowModulation(const Vertex<MeshType::meshDimension(), double>& x)
{
    return BoundaryCond::inFlowModulation(x);
}





/**
 * @brief MultiphaseFlow::SetData
 * @param initialValue
 */
template<unsigned int Dimension, typename BoundaryCond, unsigned int ... Reserve>
template<unsigned int _Dimension>
typename std::enable_if<_Dimension == 3>::type
MultiphaseFlow< Dimension, BoundaryCond, Reserve... >::setupMeshData(const std::string& fileName){

    reader = mesh.load(fileName);

    DBGVAR(mesh.getCells().size(), mesh.getVertices().size());


    mesh.template initializeCenters<METHOD_TESSELLATED>();

    mesh.setupBoundaryCells();



    mesh.setupBoundaryCellsCenters();


    auto measures = mesh.template computeElementMeasures<METHOD_TESSELLATED>();


    auto dists = ComputeCellsDistance(mesh);

    //vertToCellCon = MeshConnections<0, MeshType::meshDimension()>::connections(mesh);

    // Calculation of mesh properties
    auto faceNormals = mesh.template computeFaceNormals<METHOD_TESSELLATED>();
    meshData.allocateData(mesh);

    // setting cells volumes
    for (const auto& cell : mesh.getCells()) {
        meshData.at(cell).invVolume = 1.0 / measures.at(cell);
    }



    DBGMSG("Calculating edges properties");
    // Calculating of edge properties (length, length over cells distance, normal vector)
    for(const auto& face : mesh.getFaces()) {



        meshData.at(face).Measure = measures.at(face);

        meshData.at(face).MeasureOverDist = meshData.at(face).Measure / dists.at(face);

        meshData.at(face).n = faceNormals.at(face);

        Vertex<ProblemDimension, double>& lv = (face.getCellLeftIndex() >= BOUNDARY_INDEX(size_t)) ?
                    mesh.getBoundaryCells().at(EXTRACTING_INDEX(size_t) & face.getCellLeftIndex()).getCenter() :
                    mesh.getCells().at(face.getCellLeftIndex()).getCenter();

        Vertex<ProblemDimension, double>& rv = (face.getCellRightIndex() >= BOUNDARY_INDEX(size_t)) ?
                    mesh.getBoundaryCells().at(EXTRACTING_INDEX(size_t) & face.getCellRightIndex()).getCenter() :
                    mesh.getCells().at(face.getCellRightIndex()).getCenter();

        const Vertex<ProblemDimension, double>& faceCenter = face.getCenter();

        double leftCellKoef = 0;

        MeshApply<2,1>::apply(face.getIndex(), mesh, [&](size_t , size_t edgeIndex){
            Vertex<3,double>& vertA = mesh.getVertices().at(mesh.getEdges().at(edgeIndex).getVertexAIndex());
            Vertex<3,double>& vertB = mesh.getVertices().at(mesh.getEdges().at(edgeIndex).getVertexBIndex());

            std::array<Vertex<3,double>, 3> pyramidVec = {faceCenter - vertA, faceCenter - vertB, faceCenter};
            std::array<double, 3> norms;

            gramSchmidt<3, 3, size_t, double>(pyramidVec, norms);


            auto tessFaceCenter = (vertA + vertB + faceCenter)*(1.0/3.0);

            auto lck = (tessFaceCenter - lv).normEuclid();

            auto rck = (tessFaceCenter - rv).normEuclid();

            lck /= lck + rck;

            leftCellKoef += lck * (0.5 * norms[0] * norms[1]);

        });

        leftCellKoef /= meshData[face].Measure;

        meshData.at(face).LeftCellKoef = leftCellKoef;

        meshData.at(face).RightCellKoef = 1 - leftCellKoef;



    }

    for (size_t i = 0; i < mesh.getBoundaryCells().size(); i++) {

        Type cellType = TypeOfCell(mesh.getBoundaryCells().at(i));

        mesh.getBoundaryCells().at(i).getFlag() = cellType;
    }

}


/**
 * @brief MultiphaseFlow::SetData
 * @param initialValue
 */
template<unsigned int Dimension, typename BoundaryCond, unsigned int ... Reserve>
template<unsigned int _Dimension>
typename std::enable_if<_Dimension == 2>::type
MultiphaseFlow< Dimension, BoundaryCond, Reserve... >::setupMeshData(const std::string& fileName){

    reader = mesh.load(fileName);

    DBGVAR(mesh.getCells().size(), mesh.getVertices().size());

    mesh.initializeCenters();

    mesh.setupBoundaryCells();



    mesh.setupBoundaryCellsCenters();


    auto measures = mesh.template computeElementMeasures<>();


    auto dists = ComputeCellsDistance(mesh);

    //vertToCellCon = MeshConnections<0, MeshType::meshDimension()>::connections(mesh);

    // Calculation of mesh properties
    auto faceNormals = mesh.computeFaceNormals();
    meshData.allocateData(mesh);

    // setting cells volumes
    for (const auto& cell : mesh.getCells()) {
        meshData.at(cell).invVolume = 1.0 / measures.at(cell);
    }



    DBGMSG("Calculating edges properties");
    // Calculating of edge properties (length, length over cells distance, normal vector)
    for(const auto& face : mesh.getFaces()) {



        meshData.at(face).Measure = measures.at(face);

        meshData.at(face).MeasureOverDist = meshData.at(face).Measure / dists.at(face);

        meshData.at(face).n = faceNormals.at(face);

        Vertex<ProblemDimension, double>& lv = (face.getCellLeftIndex() >= BOUNDARY_INDEX(size_t)) ?
                    mesh.getBoundaryCells().at(EXTRACTING_INDEX(size_t) & face.getCellLeftIndex()).getCenter() :
                    mesh.getCells().at(face.getCellLeftIndex()).getCenter();

        Vertex<ProblemDimension, double>& rv = (face.getCellRightIndex() >= BOUNDARY_INDEX(size_t)) ?
                    mesh.getBoundaryCells().at(EXTRACTING_INDEX(size_t) & face.getCellRightIndex()).getCenter() :
                    mesh.getCells().at(face.getCellRightIndex()).getCenter();

        meshData.at(face).LeftCellKoef = (face.getCenter() - lv).normEuclid() / dists.at(face);

        meshData.at(face).RightCellKoef = (face.getCenter() - rv).normEuclid() / dists.at(face);

        double sum = meshData.at(face).LeftCellKoef + meshData.at(face).RightCellKoef;

        meshData.at(face).LeftCellKoef /= sum;

        meshData.at(face).RightCellKoef /= sum;
    }

    for (size_t i = 0; i < mesh.getBoundaryCells().size(); i++) {

        Type cellType = TypeOfCell(mesh.getBoundaryCells().at(i));

        mesh.getBoundaryCells().at(i).getFlag() = cellType;
    }

}







template<unsigned int Dimension, typename BoundaryCond, unsigned int ... Reserve>
Type MultiphaseFlow< Dimension, BoundaryCond, Reserve... >::TypeOfCell(const typename MeshType::Cell &cell) {
    return BoundaryCond::TypeOfCell(cell);
}

template <unsigned int Dimension, typename BoundaryCond, unsigned int ... Reserve>
void MultiphaseFlow< Dimension, BoundaryCond, Reserve... >::exportData( double time,
                                                          MeshDataContainer<ResultType, ProblemDimension>& compData,
                                                          double timeModifier)  {

    char timeStr[20];
    sprintf(timeStr, "%04ld", lround(time*timeModifier));

    std::ofstream ofile(std::string("MultiphaseFlow") + "_" + timeStr + ".vtk");
    writer.writeHeader(ofile, std::string("MPF ") + std::to_string(time));
    writer.writeToStream(ofile, mesh, reader->getCellTypes());

    VTKMeshDataWriter<ProblemDimension> dataWriter;
    dataWriter.writeToStream(ofile, compData, writer);

    ofile.close();

}



#endif // MULTIPHASEFLOW_H

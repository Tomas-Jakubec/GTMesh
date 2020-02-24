#include "multiphaseflow.h"
#include <stdexcept>








void MultiphaseFlow::calculateRHS(double, MeshDataContainer<MultiphaseFlow::ResultType, MultiphaseFlow::ProblemDimension> &compData, MeshDataContainer<MultiphaseFlow::ResultType, MultiphaseFlow::ProblemDimension> &outDeltas)
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






void MultiphaseFlow::ComputeFlux(const MeshType::Face &fcData, const MeshDataContainer<ResultType, ProblemDimension>& compData)
{
    // first compute all variables interpolated
    // in the center of the edge
    if (fcData.getCellRightIndex() < BOUNDARY_INDEX(size_t) && fcData.getCellLeftIndex() < BOUNDARY_INDEX(size_t)) {

        const FlowData<ProblemDimension> &leftData = compData.getDataByDim<ProblemDimension>().at(fcData.getCellLeftIndex());

        const FlowData<ProblemDimension> &rightData = compData.getDataByDim<ProblemDimension>().at(fcData.getCellRightIndex());

        FaceData<ProblemDimension> &currFaceData = meshData.at(fcData);

        // Computation of fluxes of density and momentum of gaseous part
        ComputeFluxGas_inner(leftData, rightData, currFaceData, fcData);

        // Compute flux of solid phase
        ComputeFluxSolid_inner(leftData, rightData, currFaceData, fcData);

    } else {
        // Applying boundary conditions
        const MeshType::Cell* innerCell = nullptr;
        const MeshType::Cell* outerCell = nullptr;
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
            ComputeFluxSolid_inflow(*innerCellData, *innerCell, *currFaceData, fcData);
            //ComputeFluxSolid_wall(*innerCellData, *currFaceData, fcData);
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

void MultiphaseFlow::ComputeViscousFlux(const MeshType::Face &fcData, const MeshDataContainer<ResultType, ProblemDimension> &compData)
{
    // first compute all variables interpolated
    // in the center of the edge
    if (fcData.getCellRightIndex() < BOUNDARY_INDEX(size_t) && fcData.getCellLeftIndex() < BOUNDARY_INDEX(size_t)) {

        const FlowData<ProblemDimension> &leftData = compData.getDataByDim<ProblemDimension>().at(fcData.getCellLeftIndex());

        const FlowData<ProblemDimension> &rightData = compData.getDataByDim<ProblemDimension>().at(fcData.getCellRightIndex());

        FaceData<ProblemDimension> &currFaceData = meshData.at(fcData);

        // Computation of fluxes of density and momentum of gaseous part
        ComputeViscousFluxGas_inner(leftData, rightData, currFaceData, fcData);

        // Compute flux of solid phase
        ComputeViscousFluxSolid_inner(leftData, rightData, currFaceData, fcData);

    } else {
        // Applying boundary conditions
        const MeshType::Cell* innerCell = nullptr;
        const MeshType::Cell* outerCell = nullptr;
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







void MultiphaseFlow::ComputeSource(const MeshType::Cell& cell,
                                   MeshDataContainer<ResultType, ProblemDimension>& compData,
                                   MeshDataContainer<MultiphaseFlow::ResultType, MultiphaseFlow::ProblemDimension> &result)
{
    FlowData<ProblemDimension>& resData = result[cell];
    FlowData<ProblemDimension>& cellData = compData[cell];


    resData.rho_g = 0;
    resData.eps_s = 0; // Firstly, the flux of rho_s is calculated. Finally, the eps_s is computed as flux_rho_s / rho_s
    resData.p_g = {};
    resData.p_s = {};

    MeshApply<ProblemDimension, ProblemDimension - 1>::apply(
                cell.getIndex(),
                mesh,
                // aplication of sum lambda to all cell faces
                [&](size_t cellIndex, size_t faceIndex){
            const FaceData<ProblemDimension>& eData = meshData.getDataByDim<ProblemDimension - 1>().at(faceIndex);

            if (cellIndex == mesh.getFaces().at(faceIndex).getCellLeftIndex()){
                resData.rho_g += eData.fluxRho_g;
                resData.eps_s += eData.fluxRho_s;
                resData.p_g += eData.fluxP_g;
                resData.p_s += eData.fluxP_s;
            } else {
                resData.rho_g -= eData.fluxRho_g;
                resData.eps_s -= eData.fluxRho_s;
                resData.p_g -= eData.fluxP_g;
                resData.p_s -= eData.fluxP_s;
            }
        }
    );

    resData.rho_g *= meshData[cell].invVolume;
    resData.eps_s *= meshData[cell].invVolume;
    resData.p_g *= meshData[cell].invVolume;
    resData.p_s *= meshData[cell].invVolume;



    Vector<ProblemDimension, double> drag = Beta_s(cellData) * (cellData.getVelocitySolid() - cellData.getVelocityGas());

    Vector<ProblemDimension,double> g_acceleration = {};

    g_acceleration[g_acceleration.size() - 1] = -9.81;

    resData.p_g += (cellData.rho_g * g_acceleration + drag);

    resData.p_s += ((rho_s - cellData.rho_g) * cellData.eps_s * g_acceleration - drag);

    // now prepare the result
    FlowData<ProblemDimension>& resultData = result.at(cell);

    resultData.eps_s /= rho_s;

    resultData.rho_g /= reg(cellData.getEps_g());

}








double MultiphaseFlow::FlowModulation(const Vertex<MeshType::meshDimension(), double>& x)
{
    double inFlowModulation = x[0];//sqrt(pow(x[0],2) + pow(x[1],2));

    //inFlowModulation = (- inFlowModulation * inFlowModulation + inFlowModulation - 0.0099) * 4.164931279;
    inFlowModulation = -(inFlowModulation -1.9) * (inFlowModulation - 4.7) * 0.5102;//(-inFlowModulation * inFlowModulation + 6.6 * inFlowModulation - 8.93) * 0.5102;
    //inFlowModulation = (inFlowModulation - 0.3) * (inFlowModulation - 0.3) * (1/0.09);
    return inFlowModulation;
}





/**
 * @brief MultiphaseFlow::SetData
 * @param initialValue
 */
void MultiphaseFlow::setupMeshData(const std::string& fileName){
    VTKMeshReader<MeshType::meshDimension()> reader;
    std::ifstream file(fileName, std::ios::binary);

    reader.loadFromStream(file, mesh);

    DBGVAR(mesh.getCells().size(), mesh.getVertices().size());


    mesh.template initializeCenters<>();

    mesh.setupBoundaryCells();



    mesh.setupBoundaryCellsCenters();


    auto measures = mesh.template computeElementMeasures<TESSELLATED>();


    auto dists = ComputeCellsDistance(mesh);

    vertToCellCon = MeshConnections<0, MeshType::meshDimension()>::connections(mesh);

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



        meshData.at(face).Length = measures.at(face);

        meshData.at(face).LengthOverDist = meshData.at(face).Length / dists.at(face);

        meshData.at(face).n = faceNormals.at(face);

        Vertex<ProblemDimension, double>& lv = (face.getCellLeftIndex() >= BOUNDARY_INDEX(size_t)) ?
                    mesh.getBoundaryCells().at(EXTRACTING_INDEX(size_t) & face.getCellLeftIndex()).getCenter() :
                    mesh.getCells().at(face.getCellLeftIndex()).getCenter();

        Vertex<ProblemDimension, double>& rv = (face.getCellRightIndex() >= BOUNDARY_INDEX(size_t)) ?
                    mesh.getBoundaryCells().at(EXTRACTING_INDEX(size_t) & face.getCellRightIndex()).getCenter() :
                    mesh.getCells().at(face.getCellRightIndex()).getCenter();

        meshData.at(face).LeftCellKoef = (face.getCenter() - lv).normEukleid() / dists.at(face);

        meshData.at(face).RightCellKoef = (face.getCenter() - rv).normEukleid() / dists.at(face);

        double sum = meshData.at(face).LeftCellKoef + meshData.at(face).RightCellKoef;

        meshData.at(face).LeftCellKoef /= sum;

        meshData.at(face).RightCellKoef /= sum;
    }

    for (size_t i = 0; i < mesh.getBoundaryCells().size(); i++) {

        Type cellType = TypeOfCell(mesh.getBoundaryCells().at(i));

        mesh.getBoundaryCells().at(i).getFlag() = cellType;
    }

}









Type MultiphaseFlow::TypeOfCell(const MeshType::Cell &cell) {
    if (cell.getIndex() >= BOUNDARY_INDEX(size_t)){

        if (
                cell.getCenter()[1] <= 1e-5
            //    sqrt(pow(cell.getCenter()[0],2) + pow(cell.getCenter()[1],2)) < 0.3 &&
            //    cell.getCenter()[2] < 0
            ) {
            return Type::INFLOW;
        }

        if (
                cell.getCenter()[1] >= 34.349
            //    sqrt(pow(cell.getCenter()[0],2) + pow(cell.getCenter()[1],2)) < 0.3 &&
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




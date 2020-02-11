#include "multiphaseflow.h"
double FlowData::rho_s = 0;
double FlowData::R_spec = 0;
double FlowData::T = 0;
#include <stdexcept>




void MultiphaseFlow::calculateRHS(double, MeshDataContainer<MultiphaseFlow::ResultType, MultiphaseFlow::ProblemDimension> &compData, MeshDataContainer<MultiphaseFlow::ResultType, MultiphaseFlow::ProblemDimension> &outDeltas)
{
    // in single compute step
    // we have to:

//#pragma omp parallel
{
    for (auto& cell : mesh.getCells()){
        FlowData& cellData = compData.at(cell);
        if (cellData.eps_s <= 0.0) {
             cellData.eps_s = 0.0;
             cellData.p_s = {};
         }
    }
    // first: compute vertexes values
    ActualizePointData(compData);
    // second: compute fluxes over edges
    for (auto& face : mesh.getFaces())
        DBGTRY(ComputeFlux(face, compData);)
    //FluxComputationParallel(*this);
//#pragma omp barrier
    //FVMethod::CellSourceComputationParallel(*this);
    for (auto& cell : mesh.getCells()){
        DBGTRY(ComputeSource(cell, compData, outDeltas);)
    }
//#pragma omp barrier
    // third: compute values of new time step
    //FVMethod::TimeStepParallel(*this);
}

}


/*
void MultiphaseFlow::ComputeStep()
{
    // in single compute step
    // we have to:

//#pragma omp parallel
{
    // first: compute vertexes values
    ActualizePointData();
    // second: compute fluxes over edges
    FluxComputationParallel(*this);
//#pragma omp barrier
    FVMethod::CellSourceComputationParallel(*this);
//#pragma omp barrier
    // third: compute values of new time step
    FVMethod::TimeStepParallel(*this);
}

}
*/



void MultiphaseFlow::ComputeFlux(const MeshType::Face &fcData, const MeshDataContainer<ResultType, ProblemDimension>& compData)
{
    // first compute all variables interpolated
    // in the center of the edge
    if (fcData.getCellRightIndex() < BOUNDARY_INDEX(size_t) && fcData.getCellLeftIndex() < BOUNDARY_INDEX(size_t)) {

        const FlowData &leftData = compData.getDataByDim<ProblemDimension>().at(fcData.getCellLeftIndex());

        const FlowData &rightData = compData.getDataByDim<ProblemDimension>().at(fcData.getCellRightIndex());

        EdgeData &currEdgeData = meshData.at(fcData);

        // Computation of fluxes of density and momentum of gaseous part
        ComputeFluxGas_inner(leftData, rightData, currEdgeData, fcData);


        // Compute flux of solid phase
        ComputeFluxSolid_inner(leftData, rightData, currEdgeData, fcData);

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

        const FlowData *innerCellData = &compData.at(*innerCell);

        EdgeData *currEdgeData = &meshData.at(fcData);

        switch(outerCell->getFlag()){
        case INFLOW : {


            ComputeFluxGas_inflow(*innerCellData, *innerCell, *currEdgeData, fcData);
            ComputeFluxSolid_inflow(*innerCellData, *innerCell, *currEdgeData, fcData);
            //ComputeFluxSolid_wall(*innerCellData, *currEdgeData, fcData);
        } break;

        case OUTFLOW :{

            ComputeFluxGas_outflow(*innerCellData, *currEdgeData, fcData);
            ComputeFluxSolid_outflow(*innerCellData, *currEdgeData, fcData);
            //ComputeFluxSolid_wall(*innerCellData, *currEdgeData, fcData);
        } break;

        case WALL : {


                ComputeFluxGas_wall(*innerCellData, *currEdgeData, fcData);
                ComputeFluxSolid_wall(*innerCellData, *currEdgeData, fcData);

            } break;

        default : {
            throw std::runtime_error("unknown cell flag: " + std::to_string(outerCell->getFlag()) +
                                     " in cell number: " + std::to_string(outerCell->getIndex()));
            }
        }
    }

}







void MultiphaseFlow::ComputeSource(const MeshType::Cell& ccData,
                                   MeshDataContainer<ResultType, ProblemDimension>& compData,
                                   MeshDataContainer<MultiphaseFlow::ResultType, MultiphaseFlow::ProblemDimension> &outDeltas)
{
    FlowData& cellData = compData.at(ccData);



    cellData.fluxRho_g = 0;
    cellData.fluxRho_s = 0;
    cellData.fluxP_g = {};
    cellData.fluxP_s = {};

    MeshApply<ProblemDimension, ProblemDimension - 1>::apply(
                ccData.getIndex(),
                mesh,
                // aplication of sum lambda to all cell faces
                [&](size_t cellIndex, size_t faceIndex){
            EdgeData& eData = meshData.getDataByDim<ProblemDimension - 1>().at(faceIndex);
            if (cellIndex == mesh.getFaces().at(faceIndex).getCellLeftIndex()){
                cellData.fluxRho_g += eData.fluxRho_g;
                cellData.fluxRho_s += eData.fluxRho_s;
                cellData.fluxP_g += eData.fluxP_g;
                cellData.fluxP_s += eData.fluxP_s;
            } else {
                cellData.fluxRho_g -= eData.fluxRho_g;
                cellData.fluxRho_s -= eData.fluxRho_s;
                cellData.fluxP_g -= eData.fluxP_g;
                cellData.fluxP_s -= eData.fluxP_s;
            }
        }
    );


    cellData.fluxP_g *= meshData.at(ccData).invVolume;
    cellData.fluxRho_g *= meshData.at(ccData).invVolume;
    cellData.fluxP_s *= meshData.at(ccData).invVolume;
    cellData.fluxRho_s *= meshData.at(ccData).invVolume;


    Vector<ProblemDimension, double> drag = Beta_s(cellData) * (cellData.getVelocitySolid() - cellData.getVelocityGas());

    Vector<ProblemDimension,double> g_acceleration = {0, -9.81};

    cellData.fluxP_g += (cellData.rho_g * g_acceleration + drag);

    cellData.fluxP_s += ((rho_s - cellData.rho_g) * cellData.eps_s * g_acceleration - drag);

    // now prepare the result
    FlowData& resultData = outDeltas.at(ccData);

    resultData.p_s = cellData.fluxP_s;
    resultData.eps_s = cellData.fluxRho_s / rho_s;


    resultData.p_g = cellData.fluxP_g;

    resultData.rho_g = cellData.fluxRho_g / reg(cellData.getEps_g());




    //cellData.T = cellData.T;
//These two can be considered as functions of other variables
    //resultData.eps_g = 1.0 - resultData.eps_s;
    //resultData.p = resultData.rho_g * R_spec * T;

}






// void MultiphaseFlow::ComputeTimeStep(FVMethod::CellComputeData &ccData)
// {
//
//     if (ccData.GetCurCell().GetCellFlag() == Type::WALL)
//         return;
//
//
//     FlowData& cellData = PrevState.GetDataAt(ccData.GetCurCellIndex());
//
//
//     // solid component
//     cellData.p_s += (cellData.fluxP_s * TimeStep);
//
//     cellData.eps_s += cellData.fluxRho_s * (TimeStep / rho_s);
//
//
//
// /*
//     cellData.u_s = cellData.u_s + cellData.fluxP_s *
//                                          (TimeStep  / reg(rho_s * cellData.eps_s)) ;
// */
//     if (cellData.eps_s <= 0.0) {
//         cellData.eps_s = 0.0;
//         cellData.p_s = Vector(0.0, 0.0);
//     }
//
//
//
//     cellData.u_s = cellData.p_s * (1 / reg(rho_s * cellData.eps_s));
// /* //I don't thing this cause problem
//     if (cellData.eps_s > 0.7){
//         cellData.eps_s = 0.7; // sand packaging limit
//     }*/
//     cellData.p_g += (cellData.fluxP_g * TimeStep);
//
//     cellData.rho_g += cellData.fluxRho_g * (TimeStep / reg(cellData.eps_g));
//
//     cellData.eps_g = 1.0 - cellData.eps_s;
//
//     cellData.u_g = cellData.p_g * (1  / reg(cellData.rho_g * cellData.eps_g));
//
//     //cellData.T = cellData.T;
//
//     cellData.p = cellData.rho_g * R_spec * T;
//
// }



double MultiphaseFlow::FlowModulation(const Vertex<MeshType::meshDimension(), double>& x)
{
    double inFlowModulation = x[0];

    //inFlowModulation = (- inFlowModulation * inFlowModulation + inFlowModulation - 0.0099) * 4.164931279;
    inFlowModulation = -(inFlowModulation -1.9) * (inFlowModulation - 4.7) * 0.5102;//(-inFlowModulation * inFlowModulation + 6.6 * inFlowModulation - 8.93) * 0.5102;
    return inFlowModulation;
}





/**
 * @brief MultiphaseFlow::SetData
 * @param initialValue
 */
void MultiphaseFlow::setupMeshData(const std::string& fileName){
    VTKMeshReader<MeshType::meshDimension()> reader;
    std::ifstream file(fileName);

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

        Vertex<2, double>& lv = (face.getCellLeftIndex() >= BOUNDARY_INDEX(size_t)) ?
                    mesh.getBoundaryCells().at(EXTRACTING_INDEX(size_t) & face.getCellLeftIndex()).getCenter() :
                    mesh.getCells().at(face.getCellLeftIndex()).getCenter();

        Vertex<2, double>& rv = (face.getCellRightIndex() >= BOUNDARY_INDEX(size_t)) ?
                    mesh.getBoundaryCells().at(EXTRACTING_INDEX(size_t) & face.getCellRightIndex()).getCenter() :
                    mesh.getCells().at(face.getCellRightIndex()).getCenter();

        meshData.at(face).LeftCellKoef = (face.getCenter() - lv).normEukleid() / dists.at(face);

        meshData.at(face).RightCellKoef = (face.getCenter() - rv).normEukleid() / dists.at(face);

    }


    InitializePointData();
}










#define BOUNDARY_NEIGHBOR 13


/**
 * @brief MultiphaseFlow::InitializePointData
 * Initializes point data this means point type
 * and first point velocity value
 */
void MultiphaseFlow::InitializePointData()
{
    PointData ini;
    ini.PointType = Type::INNER;
    ini.cellKoef = 0;

    for (PointData& pd : meshData.getDataByDim<0>()) {
        pd = ini;
    }

    /**
    ****
    **** Initialization of point types
    **** according the type of the point
    **** will decided which condition use
    ****
    **** if it is inner point then velocity
    **** in the point is average of values
    **** over cells neigboring with the point
    ****
    **** if it is outer point with
    **** - wall condition
    **** then its velocity is equal to zero
    **** - input condition
    **** then its velocity is equal to
    **** inflow velocity
    **** - output condition
    **** then the velocity is equal to average
    **** of the two neigboring cells velocities
    ****
    **** outer has higher priority than inner
    **** wall condition has the highest
    **** priority
    ****
    */
    DBGCHECK;
    for (size_t i = 0; i < mesh.getBoundaryCells().size(); i++) {

        Type cellType = TypeOfCell(mesh.getBoundaryCells().at(i));

        mesh.getBoundaryCells().at(i).getFlag() = cellType;

        size_t tmpEdgeIndex = mesh.getBoundaryCells().at(i).getBoundaryElementIndex();

        size_t tmpPointAIndex = mesh.getEdges().at(tmpEdgeIndex).getVertexAIndex();

        size_t tmpPointBIndex = mesh.getEdges().at(tmpEdgeIndex).getVertexBIndex();

        meshData.getDataByDim<0>().at(tmpPointAIndex).PointType = Type(meshData.getDataByDim<0>().at(tmpPointAIndex).PointType | cellType);

        meshData.getDataByDim<0>().at(tmpPointAIndex).cellKoef++;

        meshData.getDataByDim<0>().at(tmpPointBIndex).PointType = Type(meshData.getDataByDim<0>().at(tmpPointBIndex).PointType | cellType);

        meshData.getDataByDim<0>().at(tmpPointBIndex).cellKoef++;

        // mark the cells next to boundary
        mesh.getCells().at(mesh.getEdges().at(tmpEdgeIndex).getOtherCellIndex(i | BOUNDARY_INDEX(size_t))).getFlag() = BOUNDARY_NEIGHBOR;
    }
    // Now all points are marked with its type

    for (const auto& vert : mesh.getVertices()){
        meshData.at(vert).cellKoef += vertToCellCon.at(vert).size();
        meshData.at(vert).cellKoef = 1.0 / meshData.at(vert).cellKoef;
    }


}


void MultiphaseFlow::ActualizePointData(const MeshDataContainer<FlowData, MeshType::meshDimension()>& data) {

    // Initialize vector with zeros
    for(size_t i = 0; i < meshData.getDataByDim<0>().size(); i++) {
        meshData.getDataByDim<0>().at(i).u_g = {};
        meshData.getDataByDim<0>().at(i).u_s = {};
    }

    // can run in parallel without limitations
    for (const auto& vert : mesh.getVertices()) {
        for (const auto& cellIndex : vertToCellCon.at(vert)) {

            const MeshType::Cell& cell = mesh.getCells().at(cellIndex);

            if ((meshData.at(vert).PointType & Type::WALL) == Type::WALL){

                meshData.at(vert).u_g = {};
                meshData.at(vert).u_s = {};

            } else {

                if ((meshData.at(vert).PointType & Type::OUTFLOW) == Type::OUTFLOW){

                    if (cell.getFlag() == BOUNDARY_NEIGHBOR){
                        // if the cell over the edge is boundary then
                        // the boundary condition set in this cell is
                        // Neuman => the velocity remains same

                        meshData.at(vert).u_g += data.at(cell).getVelocityGas() * 2 * meshData.at(vert).cellKoef;

                        meshData.at(vert).u_s += data.at(cell).getVelocitySolid() * 2 * meshData.at(vert).cellKoef;

                    } else {


                        meshData.at(vert).u_g += data.at(cell).getVelocityGas() * meshData.at(vert).cellKoef;

                        meshData.at(vert).u_s += data.at(cell).getVelocitySolid() * meshData.at(vert).cellKoef;

                    }


                } else {

                    if ((meshData.at(vert).PointType & Type::INFLOW) == Type::INFLOW){

                        meshData.at(vert).u_g = inFlow_u_g * FlowModulation(vert);

                        meshData.at(vert).u_s = inFlow_u_s * FlowModulation(vert);

                    } else {

                        if ((meshData.at(vert).PointType & Type::INNER) == Type::INNER){

                            meshData.at(vert).u_g += data.at(cell).getVelocityGas() * meshData.at(vert).cellKoef;

                            meshData.at(vert).u_s += data.at(cell).getVelocitySolid() * meshData.at(vert).cellKoef;


                        }
                    }
                }
            }

        }
    }



}



/*

void MultiphaseFlow::CalculateVertexData(Mesh::MeshCell &cell, const NumData<FlowData>& cellData) {

    do{

        size_t tmpPointIndex;


        if (cell.GetCellIsLeft()){
            tmpPointIndex = cell.GetCurrentCellEdge().GetPointAIndex();
        } else {
            tmpPointIndex = cell.GetCurrentCellEdge().GetPointBIndex();
        }

        PointDatum* pDat = &PointData.at(tmpPointIndex);

        if ((pDat->PointType & Type::WALL) == Type::WALL){

            pDat->u_g = Vector(0,0);
            pDat->u_s = Vector(0,0);

        } else {

            if ((pDat->PointType & Type::OUTPUT) == Type::OUTPUT){

                if (cell.GetOtherCell().GetCellType() == TYPE_CELL_BOUNDARY){
                    // if the cell over the edge is boundary then
                    // the boundary condition set in this cell is
                    // Neuman => the velocity remains same
                    pDat->u_g += 2 * pDat->cellKoef * cellData.GetDataAt(cell.GetCurrCellIndex()).u_g;

                    pDat->u_s += 2 * pDat->cellKoef * cellData.GetDataAt(cell.GetCurrCellIndex()).u_s;

                } else {

                    pDat->u_g += pDat->cellKoef * cellData.GetDataAt(cell.GetCurrCellIndex()).u_g;

                    pDat->u_s += pDat->cellKoef * cellData.GetDataAt(cell.GetCurrCellIndex()).u_s;

                }


            } else {

                if ((pDat->PointType & Type::INPUT) == Type::INPUT){

                    pDat->u_g = inFlow.u_g * FlowModulation(GetMesh().GetPoints().at(tmpPointIndex).GetX());

                    pDat->u_s = inFlow.u_s * FlowModulation(GetMesh().GetPoints().at(tmpPointIndex).GetX());

                } else {

                    if ((pDat->PointType & Type::INNER) == Type::INNER){

                        pDat->u_g += pDat->cellKoef * cellData.GetDataAt(cell.GetCurrCellIndex()).u_g;

                        pDat->u_s += pDat->cellKoef * cellData.GetDataAt(cell.GetCurrCellIndex()).u_s;


                    }
                }
            }
        }

    // returns true if the boundary is passed thru
    } while (!cell.GotoNextCellEdge());

}

*/





Type MultiphaseFlow::TypeOfCell(const MeshType::Cell &cell) {
    if (cell.getIndex() >= BOUNDARY_INDEX(size_t)){

        if (cell.getCenter()[1] <= 1e-5) {
            return Type::INFLOW;
        }

        if (cell.getCenter()[1] >= 34.349
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




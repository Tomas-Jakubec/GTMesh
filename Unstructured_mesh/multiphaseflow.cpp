#include "multiphaseflow.h"
double FlowData::rho_s = 0;
#include <stdexcept>





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


/*
void MultiphaseFlow::ComputeFlux(FVMethod::FluxComputeData &fcData)
{

    // first compute all variables interpolated
    // in the center of the edge
    if (fcData.GetCellRight().GetCellFlag() == INNER && fcData.GetCellLeft().GetCellFlag() == INNER) {

        FlowData &leftData = PrevState.GetDataAt(fcData.GetCellLeftIndex());

        FlowData &rightData = PrevState.GetDataAt(fcData.GetCellRightIndex());

        EdgeData &currEdgeData = EdgesData.at(fcData.GetCurrEdgeIndex());

        // Computation of fluxes of density and momentum of gaseous part
        ComputeFluxGas_inner(leftData, rightData, currEdgeData, fcData);


        // Compute flux of solid phase
        ComputeFluxSolid_inner(leftData, rightData, currEdgeData, fcData);

    } else {
        // Applying boundary conditions
        const Cell* innerCell = nullptr;
        const Cell* outerCell = nullptr;
        if (fcData.GetCellRight().GetCellFlag() == INNER ){
            innerCell = &fcData.GetCellRight();
            outerCell = &fcData.GetCellLeft();
        } else {
            innerCell = &fcData.GetCellLeft();
            outerCell = &fcData.GetCellRight();
        }

        FlowData *innerCellData = &PrevState.GetDataAt(innerCell->GetCellIndex());

        EdgeData *currEdgeData = &EdgesData.at(fcData.GetCurrEdgeIndex());

        switch(outerCell->GetCellFlag()){
        case INPUT : {


            ComputeFluxGas_inflow(*innerCellData, *innerCell, *currEdgeData, fcData);
            ComputeFluxSolid_inflow(*innerCellData, *innerCell, *currEdgeData, fcData);
            //ComputeFluxSolid_wall(*innerCellData, *currEdgeData, fcData);
        } break;

        case OUTPUT :{

            ComputeFluxGas_outflow(*innerCellData, *currEdgeData, fcData);
            ComputeFluxSolid_outflow(*innerCellData, *currEdgeData, fcData);
            //ComputeFluxSolid_wall(*innerCellData, *currEdgeData, fcData);
        } break;

        case WALL : {


                ComputeFluxGas_wall(*innerCellData, *currEdgeData, fcData);
                ComputeFluxSolid_wall(*innerCellData, *currEdgeData, fcData);

            } break;

        default : {
            throw std::runtime_error("unknown cell flag: " + std::to_string(outerCell->GetCellFlag()) +
                                     " in cell number: " + std::to_string(outerCell->GetCellIndex()));
            }
        }
    }

}
*/





/*
void MultiphaseFlow::ComputeSource(FVMethod::CellComputeData &ccData)
{
    FlowData& cellData = PrevState.GetDataAt(ccData.GetCurCellIndex());


    double invCellVol = 1.0 / (CellsVolumes.at(ccData.GetCurCellIndex()));

    cellData.fluxRho_g = 0;
    cellData.fluxRho_s = 0;
    cellData.fluxP_g = Vector(0.0,0.0);
    cellData.fluxP_s = Vector(0.0,0.0);

    size_t tmpEI = ccData.GetCurCell().GetCellEdgeIndex();
    do {
        if (ccData.GetCurCellIndex() == GetMesh().GetEdges().at(tmpEI).GetCellLeftIndex()){
            cellData.fluxRho_g += EdgesData.at(tmpEI).fluxRho_g;
            cellData.fluxRho_s += EdgesData.at(tmpEI).fluxRho_s;
            cellData.fluxP_g += EdgesData.at(tmpEI).fluxP_g;
            cellData.fluxP_s += EdgesData.at(tmpEI).fluxP_s;
        } else {
            cellData.fluxRho_g -= EdgesData.at(tmpEI).fluxRho_g;
            cellData.fluxRho_s -= EdgesData.at(tmpEI).fluxRho_s;
            cellData.fluxP_g -= EdgesData.at(tmpEI).fluxP_g;
            cellData.fluxP_s -= EdgesData.at(tmpEI).fluxP_s;
        }
        tmpEI = GetMesh().GetEdges().at(tmpEI).GetNextEdge(ccData.GetCurCellIndex());
    } while (tmpEI != ccData.GetCurCell().GetCellEdgeIndex());


    cellData.fluxP_g *= invCellVol;
    cellData.fluxRho_g *= invCellVol;
    cellData.fluxP_s *= invCellVol;
    cellData.fluxRho_s *= invCellVol;


    Vector drag = Beta_s(cellData) * (cellData.u_s - cellData.u_g);

    cellData.fluxP_g += (cellData.rho_g * Vector(0, -9.81) + drag);

    cellData.fluxP_s += ((rho_s - cellData.rho_g) * cellData.eps_s * Vector(0, -9.81) - drag);

}
*/






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



double MultiphaseFlow::FlowModulation(double x)
{
    double inFlowModulation = x;

    inFlowModulation = (- inFlowModulation * inFlowModulation + inFlowModulation - 0.0099) * 4.164931279;
    //inFlowModulation = -(inFlowModulation -1.9) * (inFlowModulation - 4.7) * 0.5102;//(-inFlowModulation * inFlowModulation + 6.6 * inFlowModulation - 8.93) * 0.5102;
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

        DBGVAR_CSV(face.getCellLeftIndex(), face.getCellRightIndex(), EXTRACTING_INDEX(size_t) & face.getCellLeftIndex(), EXTRACTING_INDEX(size_t) & face.getCellRightIndex());

        Vertex<2, double>& lv = (face.getCellLeftIndex() >= BOUNDARY_INDEX(size_t)) ?
                    mesh.getBoundaryCells().at(EXTRACTING_INDEX(size_t) & face.getCellLeftIndex()).getCenter() :
                    mesh.getCells().at(face.getCellLeftIndex()).getCenter();

        Vertex<2, double>& rv = (face.getCellRightIndex() >= BOUNDARY_INDEX(size_t)) ?
                    mesh.getBoundaryCells().at(EXTRACTING_INDEX(size_t) & face.getCellRightIndex()).getCenter() :
                    mesh.getCells().at(face.getCellRightIndex()).getCenter();

        meshData.at(face).LeftCellKoef = (face.getCenter() - lv).normEukleid() / dists.at(face);

        meshData.at(face).RightCellKoef = (face.getCenter() - rv).normEukleid() / dists.at(face);

    }

}











///**
// * @brief MultiphaseFlow::InitializePointData
// * Initializes point data this means point type
// * and first point velocity value
// */
//void MultiphaseFlow::InitializePointData()
//{
//    PointData ini;
//    ini.PointType = Type::INNER;
//    ini.cellKoef = 0;
//
//    PointData.resize(GetMesh().GetPoints().size(), ini);
//
//    /**
//    ****
//    **** Initialization of point types
//    **** according the type of the point
//    **** will decided which condition use
//    ****
//    **** if it is inner point then velocity
//    **** in the point is average of values
//    **** over cells neigboring with the point
//    ****
//    **** if it is outer point with
//    **** - wall condition
//    **** then its velocity is equal to zero
//    **** - input condition
//    **** then its velocity is equal to
//    **** inflow velocity
//    **** - output condition
//    **** then the velocity is equal to average
//    **** of the two neigboring cells velocities
//    ****
//    **** outer has higher priority than inner
//    **** wall condition has the highest
//    **** priority
//    ****
//    */
//    for (size_t i = 0; i < mesh.getBoundaryCells().size(); i++) {
//
//        Type cellType = TypeOfCell(mesh().getBoundaryCells().at(i));
//
//        size_t tmpEdgeIndex = GetMesh().GetBoundaryCells().at(i).GetCellEdgeIndex();
//
//        size_t tmpPointAIndex = GetMesh().GetEdges().at(tmpEdgeIndex).GetPointAIndex();
//
//        size_t tmpPointBIndex = GetMesh().GetEdges().at(tmpEdgeIndex).GetPointBIndex();
//
//        PointData.at(tmpPointAIndex).PointType = Type(PointData.at(tmpPointAIndex).PointType | cellType);
//
//        PointData.at(tmpPointAIndex).cellKoef++;
//
//        PointData.at(tmpPointBIndex).PointType = Type(PointData.at(tmpPointBIndex).PointType | cellType);
//
//        PointData.at(tmpPointBIndex).cellKoef++;
//    }
//    // Now all points are marked with its type
//
//    auto connections = MeshConnections<0,2>::connections(mesh);
//
//    for (const auto& vert : mesh.getVertices()){
//        meshData.at(vert). cellKoef = 1.0 / connections.at(vert).size();
//    }
//
//    DBGTRY(ActualizePointData());
//
//}

/*
void MultiphaseFlow::ActualizePointData() {

    // Initialize vector with zeros
//    #pragma omp for
    for(size_t i = 0; i < PointData.size(); i++) {
        PointData.at(i).u_g = Vector();
        PointData.at(i).u_s = Vector();
    }

    // Calculate vertex data of each inner cell
    for(auto& one_colour : CellToVertexColoring){
//        #pragma omp for
        for(size_t index = 0; index < one_colour.size(); index++){

            Mesh::MeshCell cell = GetMesh().GetMeshCell(one_colour.at(index));

            CalculateVertexData(cell, PrevState);

        }

    }
}

*/




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



/*

MultiphaseFlow::Type MultiphaseFlow::TypeOfCell(const MeshType::Cell &cell) {
    if (cell.getIndex() > BOUNDARY_INDEX){

        if (cell.getCenter()[1] <= 1e-5) {
            return Type::INPUT;
        }

        if (//cell.GetCenter().GetY() >= 34.349
            //cell.GetCenter().GetX() > 2 && (cell.GetCenter().GetY() < 32.0 && cell.GetCenter().GetY() > 30)
            //cell.GetCenter().GetX() > 8 && cell.GetCenter().GetY() >= 33.39
            cell.getCenter()[1] >= 1.999
           ) {
            DBGVAR(cell.getCenter());
            return Type::OUTPUT;
        }

        return Type::WALL;

    } else {

        return Type(cell.getFlag());

    }


    throw(std::runtime_error ("cell type not recognized " + std::to_string(cell.getIndex())));
}

*/





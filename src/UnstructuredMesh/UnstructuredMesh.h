#ifndef UNSTRUCTUREDMESH_H
#define UNSTRUCTUREDMESH_H
#include "MeshElements/MeshElements.h"
#include "vector"
#include <type_traits>
#include "MeshFunctions/MeshFunctions.h"
#include "MeshIO/MeshReader/FPMAMeshReader.h"
#include "MeshIO/MeshReader/VTKMeshReader.h"



/**
 * @brief The UnstructuredMesh class is the wrapper of MeshElements.
 * This class provides functionality of the mesh. The template parameters
 * are: <BR>
 * MeshDimension: the dimension of the mesh.<BR>
 * IndexType: the type of references in the mesh (integer type).<BR>
 * Real: the type of coordinates of the vertices.<BR>
 * Reserve: variadic parameter declaring reserved number of subelements of
 * elements with dimension higher than 1 and less than MeshDimension. The default
 * value is 0 which resolves in dynamic allocation of subelements.
 */
template <unsigned int MeshDimension, typename IndexType, typename Real, unsigned int ...Reserve>
class UnstructuredMesh : public MeshElements<MeshDimension, IndexType, Real, Reserve...>{

public:
    /**
     * @brief Initializes the centers of cells and faces of the mesh.
     * Using the method parameter the method of calculation can be chosen.
     * The default method is averaging. The METHOD_TESELATED calculates
     * the centers of faces in 3D by averaging using the measures of
     * partial triangles.
     */
    template<ComputationMethod Method = ComputationMethod::METHOD_DEFAULT>
    void initializeCenters(){
        auto centers = computeCenters<Method>(*this);

        for (auto& face : this->getFaces()){
            face.setCenter(centers[face]);
        }
        for (auto& cell : this->getCells()){
            cell.setCenter(centers[cell]);
        }
    }

    /**
     * @brief Calculates the measure of the elements with dimenson
     * higher than one.
     */
    template<ComputationMethod Method = ComputationMethod::METHOD_DEFAULT>
    MakeMeshDataContainer_t<Real, make_custom_integer_sequence_t<unsigned int, 1, MeshDimension>> computeElementMeasures() {
        return computeMeasures<Method>(*this);
    }

    /**
     * @brief Calculates the normal vectors of the faces in the mesh.
     */
    template<ComputationMethod Method = ComputationMethod::METHOD_DEFAULT>
    MeshDataContainer<Vector<MeshDimension, Real>, MeshDimension-1> computeFaceNormals() {
        return ::computeFaceNormals<Method>(*this);
    }

    /**
     * @brief Applies the passed function func on the elements of target dim.
     * The passed function must accept two parameters of IndexType. The first
     * will be used for the origin element index (where the loop starts). The
     * second is the reached element index. Note that the connected elements
     * might be visited more than once.
     * @param func the function to be applied
     * @example
     * @code
     * mesh.apply(
     *      [](IndexType oI, IndexType sI){
     *          std::cout << "origin: " << oI <<" connected: "<< sI;
     *      }
     * );
     * @endcode
     */
    template<unsigned int StartDim, unsigned int TargetDim, typename Functor>
    void apply(const Functor& func) {
        return MeshApply<StartDim, TargetDim>::apply(*this, func);
    }

    /**
     * @brief Applies the passed function func on the elements of target dim
     * connected to the element with index startElementIndex.
     * The passed function must accept two parameters of IndexType. The first
     * will be used for the origin element index (startElementIndex). The
     * second is the reached/connected element index. Note that the connected
     * elements might be visited more than once.
     * @param func the function to be applied
     * @example
     * @code
     * mesh.apply(
     *      42,
     *      [](IndexType oI, IndexType sI){
     *          std::cout << "origin: " << oI <<" connected: "<< sI;
     *      }
     * );
     * @endcode
     */
    template<unsigned int StartDim, unsigned int TargetDim, typename Functor>
    void apply(IndexType startElementIndex,const Functor& func) {
        return MeshApply<StartDim, TargetDim>::apply(startElementIndex, *this, func);
    }

    /**
     * @brief Returns the indexes of connected elements of
     * dimesnion TargetDim to elements with dimension
     * StartDim.
     */
    template<unsigned int StartDim, unsigned int TargetDim, Order ConnectionsOrder = ORDER_ASCEND>
    MeshDataContainer<std::vector<IndexType>, StartDim> connections() {
        return MeshConnections<StartDim, TargetDim, ConnectionsOrder>::connections(*this);
    }

    /**
     * @brief Returns the indexes of elements with dimesnion
     * ConnectedDim that neighbors with elements from StartDim.
     * The connection is determined over the elements from ConnectingDim.
     */
    template<unsigned int StartDim, unsigned int ConnectingDim, unsigned int ConnectedDim = StartDim, Order ConnectionsOrder = Order::ORDER_ASCEND>
    MeshDataContainer<std::vector<IndexType>, StartDim> neighborhood() {
        return MeshNeighborhood<StartDim, ConnectingDim, ConnectedDim, ConnectionsOrder>::neighbors(*this);
    }


    /**
     * @brief Colors the mesh elements of dimension StartDim with connection over
     * ConnectingDim.
     */
    template<unsigned int StartDim, unsigned int ConnectingDim,  ColoringMethod Method = ColoringMethod::METHOD_GREEDY>
    typename std::enable_if<Method == METHOD_GREEDY,MeshDataContainer<unsigned int, StartDim>>::type
    coloring() {
        return MeshColoring<StartDim, ConnectingDim, Method>::color(*this);
    }

    /**
     * @brief Colors the mesh elements of dimension StartDim with connection over
     * ConnectingDim. Moreover, this function allows to set a seed for the random
     * reselection of the colors. Random reselection provides more even distribution
     * of the colors.
     */
    template<unsigned int StartDim, unsigned int ConnectingDim,  ColoringMethod Method = ColoringMethod::METHOD_GREEDY>
    typename std::enable_if<Method == METHOD_RANDOM,MeshDataContainer<unsigned int, StartDim>>::type
    coloring(unsigned int seed = 1562315) {
        return MeshColoring<StartDim, ConnectingDim, Method>::color(*this, seed);
    }

/*
    Real CalculateFaceMeasureOverCellDist(IndexType edgeIndex, MeshDataContainer){

        const auto& edge = this->getEdges().at(edgeIndex);
        return CalculateEdgeMeasure(edgeIndex) / CalculateCellDist(edge.GetCellLeftIndex(), edge.GetCellRightIndex());

    }
  */
    /**
     * @brief Loads the mesh from a file passed as parameter filePath.
     * For mesh with dimension 3, there are 2 provided formats VTK and FPMA.
     */
    template<unsigned int MeshDim = MeshDimension>
    typename std::enable_if<MeshDim == 3,std::unique_ptr<MeshReader<MeshDimension>>>::type
    load(const std::string& filePath){

        typedef std::unique_ptr<MeshReader<MeshDimension>> retType;
        retType reader_ptr;

        std::ifstream file(filePath, std::ios::binary);
        if (!file.is_open()) {
            throw std::runtime_error("was not able to open file: \"" + filePath + "\"");
        }

        if (filePath.find(".vtk") != filePath.npos){
            DBGMSG("file recognized as VTK");
            auto reader = new VTKMeshReader<MeshDimension>();
            reader->loadFromStream(file, *this);
            reader_ptr = retType(reader);
        }

        if (filePath.find(".fpma") != filePath.npos){
            DBGMSG("file recognized as FPMA");
            auto reader = new FPMAMeshReader<MeshDimension>();
            reader->loadFromStream(file, *this);
            reader_ptr = retType(reader);
        }

        file.close();

        return reader_ptr;
    }

    /**
     * @brief Loads the mesh from a file passed as parameter filePath.
     * For mesh with dimension 2, only VTK format is supported so far.
     */
    template<unsigned int MeshDim = MeshDimension>
    typename std::enable_if<MeshDim == 2,std::unique_ptr<MeshReader<MeshDimension>>>::type
    load(const std::string& filePath){

        typedef std::unique_ptr<MeshReader<MeshDimension>> retType;
        retType reader_ptr;

        std::ifstream file(filePath, std::ios::binary);
        if (!file.is_open()) {
            throw std::runtime_error("was not able to open file: \"" + filePath + "\"");
        }

        if (filePath.find(".vtk") != filePath.npos){
            DBGMSG("file recognized as VTK");
            auto reader = new VTKMeshReader<MeshDimension>();
            reader->loadFromStream(file, *this);
            reader_ptr = retType(reader);
        }

        file.close();

        return reader_ptr;
    }
};

#endif // UNSTRUCTUREDMESH_H

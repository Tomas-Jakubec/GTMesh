#ifndef UNSTRUCTUREDMESH_H
#define UNSTRUCTUREDMESH_H
#include "MeshElements/MeshElements.h"
#include "vector"
#include <type_traits>
#include "MeshFunctions/MeshFunctions.h"
#include "MeshIO/MeshReader/FPMAMeshReader.h"
#include "MeshIO/MeshReader/VTKMeshReader.h"

#include "MeshIO/MeshWriter/FPMAMeshWriter.h"
#include "MeshIO/MeshWriter/VTKMeshWriter.h"


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
     * @brief Sets the centers up according to the ones passed by parameter.
     */
    void initializeCenters(MakeMeshDataContainer_t<Vertex<MeshDimension, Real>, make_custom_integer_sequence_t<unsigned int, 1,MeshDimension>> centers){
        for (auto& face : this->getFaces()){
            face.setCenter(centers[face]);
        }
        for (auto& cell : this->getCells()){
            cell.setCenter(centers[cell]);
        }
    }

    /**
     * @brief Calculates the centers of all elements in the mesh
     * with dimension higher or equal to one.
     * The METHOD_TESELATED calculates the centers of faces in
     * 3D by averaging using the measures of partial triangles.
     */
    template<ComputationMethod Method = ComputationMethod::METHOD_DEFAULT>
    auto computeElementCenters(){
        return computeCenters<Method>(*this);
    }

    /**
     * @brief Calculates the measure of the elements with dimenson
     * higher than one.
     */
    template<ComputationMethod Method = ComputationMethod::METHOD_DEFAULT, typename ..., unsigned int MD = MeshDimension, typename std::enable_if< (MD <= 3) , bool >::type = true>
    MakeMeshDataContainer_t<Real, make_custom_integer_sequence_t<unsigned int, 1, MeshDimension>> computeElementMeasures() const {
        return computeMeasures<Method>(*this);
    }

    /**
     * @brief Calculates the measure of the elements with dimenson
     * higher than one.
     */
    template<ComputationMethod Method = ComputationMethod::METHOD_DEFAULT, typename ..., unsigned int MD = MeshDimension, typename std::enable_if< (MD > 3) , bool >::type = true>
    MakeMeshDataContainer_t<Real, make_custom_integer_sequence_t<unsigned int, 1, MeshDimension>> computeElementMeasures() const {
        return computeMeasures<Method>(computeCenters<Method>(*this), *this);
    }

    /**
     * @brief Calculates the normal vectors of the faces in the mesh.
     */
    template<ComputationMethod Method = ComputationMethod::METHOD_DEFAULT, typename ..., unsigned int MD = MeshDimension, typename std::enable_if< (MD <= 3) , bool >::type = true>
    MeshDataContainer<Vector<MeshDimension, Real>, MeshDimension-1> computeFaceNormals() const {
        return ::computeFaceNormals<Method>(*this);
    }

    /**
     * @brief Calculates the normal vectors of the faces in the mesh.
     */
    template<ComputationMethod Method = ComputationMethod::METHOD_DEFAULT, typename ..., unsigned int MD = MeshDimension, typename std::enable_if< (MD > 3) , bool >::type = true>
    MeshDataContainer<Vector<MeshDimension, Real>, MeshDimension-1> computeFaceNormals() const {
        return ::computeFaceNormals<Method>(computeCenters<Method>(*this), *this);
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
    void apply(const Functor& func) const {
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
    void apply(IndexType startElementIndex,const Functor& func) const {
        return MeshApply<StartDim, TargetDim>::apply(startElementIndex, *this, func);
    }

    /**
     * @brief Returns the indexes of connected elements of
     * dimesnion TargetDim to elements with dimension
     * StartDim.
     */
    template<unsigned int StartDim, unsigned int TargetDim, Order ConnectionsOrder = ORDER_ASCEND>
    MeshDataContainer<std::vector<IndexType>, StartDim> connections() const {
        return MeshConnections<StartDim, TargetDim, ConnectionsOrder>::connections(*this);
    }

    /**
     * @brief Returns the indexes of elements with dimesnion
     * ConnectedDim that neighbors with elements from StartDim.
     * The connection is determined over the elements from ConnectingDim.
     */
    template<unsigned int StartDim, unsigned int ConnectingDim, unsigned int ConnectedDim = StartDim, Order ConnectionsOrder = Order::ORDER_ASCEND>
    MeshDataContainer<std::vector<IndexType>, StartDim> neighborhood() const {
        return MeshNeighborhood<StartDim, ConnectingDim, ConnectedDim, ConnectionsOrder>::neighbors(*this);
    }


    /**
     * @brief Colors the mesh elements of dimension StartDim with connection over
     * ConnectingDim.
     */
    template<unsigned int StartDim, unsigned int ConnectingDim,  ColoringMethod Method = ColoringMethod::METHOD_GREEDY>
    typename std::enable_if<Method == METHOD_GREEDY,MeshDataContainer<unsigned int, StartDim>>::type
    coloring() const {
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
    coloring(unsigned int seed = 1562315) const {
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
    template<typename...T, unsigned int MeshDim = MeshDimension, typename std::enable_if<MeshDim == 3, bool>::type = true>
    std::unique_ptr<MeshReader<MeshDimension>>
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
     * @brief Loads the mesh from a file passed as the parameter filePath.
     * For the mesh with dimension 2, only VTK format is supported so far.
     */
    template<typename...T, unsigned int MeshDim = MeshDimension, typename std::enable_if<MeshDim == 2, bool>::type = true>
    std::unique_ptr<MeshReader<MeshDimension>>
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


    /**
     * @brief Writes the mesh from a file passed as parameter filePath.
     * For mesh with dimension 3, there are 2 provided formats VTK and FPMA.
     * For the mesh with dimension 2, only VTK format is supported so far.
     */
    std::unique_ptr<MeshWriter<MeshDimension>>
    write(const std::string& filePath,
          const MeshReader<MeshDimension>& meshReader,
          const std::string& dataHeader = ""){
        return write(filePath, meshReader.getCellTypes(), dataHeader);
    }

    /**
     * @brief Writes the mesh from a file passed as parameter filePath.
     * For mesh with dimension 3, there are 2 provided formats VTK and FPMA.
     * For the mesh with dimension 2, only VTK format is supported so far.
     */
    void
    write(const std::string& filePath,
          std::unique_ptr<MeshWriter<MeshDimension>>& writer,
          const MeshReader<MeshDimension>& meshReader,
          const std::string& dataHeader = ""){
        return write(filePath, writer, meshReader.getCellTypes(), dataHeader);
    }
protected:
    template<unsigned int MeshDim = MeshDimension, typename std::enable_if<MeshDim == 3, bool>::type = true>
    auto PolytopeType() {
        return MeshReader<MeshDim>::elementType::ElementType::POLYHEDRON;
    };
    template<unsigned int MeshDim = MeshDimension, typename std::enable_if<MeshDim == 2, bool>::type = true>
    auto PolytopeType() {
        return MeshReader<MeshDim>::elementType::ElementType::POLYGON;
    };


public:
    /**
     * @brief Writes the mesh from a file passed as parameter filePath.
     * For mesh with dimension 3, there are 2 provided formats VTK and FPMA.
     * For the mesh with dimension 2, only VTK format is supported so far.
     * @return unique_ptr to an instance of MeshWriter which has indexed the mesh.
     */
    std::unique_ptr<MeshWriter<MeshDimension>>
    write(const std::string& filePath,
          const std::string& dataHeader = ""){
        using ET = typename MeshReader<MeshDimension>::elementType::ElementType;
        MeshDataContainer<ET, MeshDimension> cellTypes;
        cellTypes.allocateData(*this, PolytopeType());
        return write(filePath, cellTypes, dataHeader);
    }
private:
    /**
     * @brief Returns true if the filePath string ends with the string ext
     */
    bool checkExtension(const std::string& filePath, const std::string& ext){
        std::string fileExt = filePath.substr(filePath.find_last_of("."));
        return fileExt == ext;
    }
public:
    /**
     * @brief Writes the mesh from a file passed as parameter filePath.
     * For mesh with dimension 3, there are 2 provided formats VTK and FPMA.
     * For the mesh with dimension 2, only VTK format is supported so far.
     * @return unique_ptr to an instance of MeshWriter which has indexed the mesh.
     */
    template<typename...T, unsigned int MeshDim = MeshDimension, typename std::enable_if<MeshDim == 3, bool>::type = true>
    std::unique_ptr<MeshWriter<MeshDimension>>
    write(const std::string& filePath,
          const MeshDataContainer<typename MeshWriter<MeshDimension>::elementType::ElementType, MeshDim>& cellTypes,
          const std::string& dataHeader = ""){

        typedef std::unique_ptr<MeshWriter<MeshDimension>> retType;
        retType writer_ptr;

        if (checkExtension(filePath, ".vtk")){
            DBGMSG("file recognized as VTK");
            writer_ptr = std::make_unique<VTKMeshWriter<MeshDimension>>();
            write(filePath, writer_ptr, cellTypes, dataHeader);
        }

        if (checkExtension(filePath, ".fpma")){
            DBGMSG("file recognized as FPMA");
            writer_ptr = std::make_unique<FPMAMeshWriter<MeshDimension>>();
            write(filePath, writer_ptr, cellTypes, dataHeader);
        }

        return writer_ptr;
    }


    /**
     * @brief Writes the mesh to a file passed as parameter filePath.
     * For mesh with dimension 3, there are 2 provided formats VTK and FPMA.
     * This overload accepts a writer parameter which is a unique_ptr<MeshWriter<MeshDimension>>.
     * The writer is utilized to export the mesh. If the writer type and the format of
     * the export file does not match then the writer is corrected.
     * @param writer instance of MeshWriter (e.g. VTKMeshWriter).
     */
    template<typename...T, unsigned int MeshDim = MeshDimension, typename std::enable_if<MeshDim == 3, bool>::type = true>
    void
    write(const std::string& filePath,
          std::unique_ptr<MeshWriter<MeshDimension>>& writer,
          const MeshDataContainer<typename MeshWriter<MeshDimension>::elementType::ElementType, MeshDim>& cellTypes = MeshDataContainer<typename MeshWriter<MeshDimension>::elementType::ElementType, MeshDim>(),
          const std::string& dataHeader = ""){
        if (!writer) {
            writer = write(filePath, cellTypes, dataHeader);
        } else {

            std::ofstream file(filePath, std::ios::binary);
            if (!file.is_open()) {
                throw std::runtime_error("was not able to open file: \"" + filePath + "\"");
            }

            if (typeid (*writer.get()) == typeid (VTKMeshWriter<MeshDimension>)){
                if (!checkExtension(filePath, ".vtk")){
                    DBGMSG("The writer type and file name does not match!");
                    writer = write(filePath, cellTypes, dataHeader);
                } else {
                    DBGMSG("writer recognized as VTK");
                    auto writer_loc = dynamic_cast<VTKMeshWriter<MeshDimension>*>(writer.get());
                    writer_loc->writeHeader(file, dataHeader);
                    if (cellTypes.template getDataByPos<0>().size() == 0) {
                        MeshDataContainer<typename MeshWriter<MeshDimension>::elementType::ElementType, MeshDim> locCellTypes(*this, PolytopeType());
                        writer_loc->writeToStream(file, *this, locCellTypes);
                    } else {
                        writer_loc->writeToStream(file, *this, cellTypes);
                    }
                }
            }

            else if (typeid (*writer.get()).hash_code() == typeid (FPMAMeshWriter<MeshDimension>).hash_code()){
                if (!checkExtension(filePath, ".fpma")){
                    DBGMSG("The writer type and file name does not match!");
                    writer = write(filePath, cellTypes, dataHeader);
                } else {
                    DBGMSG("writer recognized as FPMA");
                    auto reader = dynamic_cast<FPMAMeshWriter<MeshDimension>*>(writer.get());
                    reader->writeToStream(file, *this);
                }
            }
            file.close();
        }
    }


    /**
     * @brief Writes the mesh to a file passed as parameter filePath.
     * For mesh with dimension 3, there are 2 provided formats VTK and FPMA.
     * This overload accepts a writer parameter which is a unique_ptr<MeshWriter<MeshDimension>>.
     * The writer is utilized to export the mesh. If the writer type and the format of
     * the export file does not match then the writer is corrected.
     * @param writer instance of MeshWriter (e.g. VTKMeshWriter).
     */
    template<typename...T, unsigned int MeshDim = MeshDimension, typename std::enable_if<MeshDim == 2, bool>::type = true>
    std::unique_ptr<MeshWriter<MeshDimension>>
    write(const std::string& filePath,
          const MeshDataContainer<typename MeshWriter<MeshDimension>::elementType::ElementType, MeshDim>& cellTypes,
          const std::string& dataHeader = ""){

        typedef std::unique_ptr<MeshWriter<MeshDimension>> retType;
        retType writer_ptr;

        if (checkExtension(filePath, ".vtk")){
            DBGMSG("file recognized as VTK");
            writer_ptr = std::make_unique<VTKMeshWriter<MeshDimension>>();
            write(filePath, writer_ptr, cellTypes, dataHeader);
        }

        return writer_ptr;
    }

    /**
     * @brief Writer the mesh to a file passed as the parameter filePath.
     * For the mesh with dimension 2, only VTK format is supported so far.
     */
    template<typename...T, unsigned int MeshDim = MeshDimension, typename std::enable_if<MeshDim == 2, bool>::type = true>
    void
    write(const std::string& filePath,
          std::unique_ptr<MeshWriter<MeshDimension>>& writer,
          const MeshDataContainer<typename MeshWriter<MeshDimension>::elementType::ElementType, MeshDim>& cellTypes = MeshDataContainer<typename MeshWriter<MeshDimension>::elementType::ElementType, MeshDim>(),
          const std::string& dataHeader = ""){
        if (!writer) {
            writer = write(filePath, cellTypes, dataHeader);
        } else {

            std::ofstream file(filePath, std::ios::binary);
            if (!file.is_open()) {
                throw std::runtime_error("was not able to open file: \"" + filePath + "\"");
            }

            if (typeid (*writer.get()) == typeid (VTKMeshWriter<MeshDimension>)){
                if (checkExtension(filePath, ".vtk")){
                    DBGMSG("writer recognized as VTK");
                    auto writer_loc = dynamic_cast<VTKMeshWriter<MeshDimension>*>(writer.get());
                    writer_loc->writeHeader(file, dataHeader);
                    if (cellTypes.template getDataByPos<0>().size() == 0) {
                        MeshDataContainer<typename MeshWriter<MeshDimension>::elementType::ElementType, MeshDim> locCellTypes(*this, PolytopeType());
                        writer_loc->writeToStream(file, *this, locCellTypes);
                    } else {
                        writer_loc->writeToStream(file, *this, cellTypes);
                    }
                } else {
                    DBGMSG("The writer type and file name does not match!");
                    writer = write(filePath, cellTypes, dataHeader);
                }
            }

            file.close();
        }

    }


};

#endif // UNSTRUCTUREDMESH_H

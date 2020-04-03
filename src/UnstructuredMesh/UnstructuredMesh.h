#ifndef UNSTRUCTUREDMESH_H
#define UNSTRUCTUREDMESH_H
#include "MeshElements/MeshElements.h"
#include "vector"
#include <type_traits>
#include "MeshFunctions/MeshFunctions.h"
#include "MeshIO/MeshReader/FPMAMeshReader.h"
#include "MeshIO/MeshReader/VTKMeshReader.h"



template <unsigned int MeshDimension, typename IndexType, typename Real, unsigned int ...Reserve>
class UnstructuredMesh : public MeshElements<MeshDimension, IndexType, Real, Reserve...>{

public:
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

    template<ComputationMethod Method = ComputationMethod::METHOD_DEFAULT>
    MakeMeshDataContainer_t<Real, make_custom_integer_sequence_t<unsigned int, 1, MeshDimension>> computeElementMeasures() {
        return computeMeasures<Method>(*this);
    }

    template<ComputationMethod Method = ComputationMethod::METHOD_DEFAULT>
    MeshDataContainer<Vector<MeshDimension, Real>, MeshDimension-1> computeFaceNormals() {
        return ::computeFaceNormals<Method>(*this);
    }

    template<unsigned int StartDim, unsigned int TargetDim, typename Functor>
    void apply(const Functor& func) {
        return MeshApply<StartDim, TargetDim>::apply(*this, func);
    }


    template<unsigned int StartDim, unsigned int TargetDim, typename Functor>
    void apply(IndexType startElementIndex,const Functor& func) {
        return MeshApply<StartDim, TargetDim>::apply(startElementIndex, *this, func);
    }

    template<unsigned int StartDim, unsigned int TargetDim, Order ConnectionsOrder = ORDER_ASCEND>
    MeshDataContainer<std::vector<IndexType>, StartDim> connections() {
        return MeshConnections<StartDim, TargetDim, ConnectionsOrder>::connections(*this);
    }

    template<unsigned int StartDim, unsigned int ConnectingDim, unsigned int ConnectedDim = StartDim, Order ConnectionsOrder = Order::ORDER_ASCEND>
    MeshDataContainer<std::vector<IndexType>, StartDim> neighborhood() {
        return MeshNeighborhood<StartDim, ConnectingDim, ConnectedDim, ConnectionsOrder>::neighbors(*this);
    }


    template<unsigned int StartDim, unsigned int ConnectingDim,  ColoringMethod Method = ColoringMethod::METHOD_GREEDY>
    typename std::enable_if<Method == METHOD_GREEDY,MeshDataContainer<unsigned int, StartDim>>::type
    coloring() {
        return ColorMesh<StartDim, ConnectingDim, Method>::color(*this);
    }

    template<unsigned int StartDim, unsigned int ConnectingDim,  ColoringMethod Method = ColoringMethod::METHOD_GREEDY>
    typename std::enable_if<Method == METHOD_RANDOM,MeshDataContainer<unsigned int, StartDim>>::type
    coloring(unsigned int seed = 1562315) {
        return ColorMesh<StartDim, ConnectingDim, Method>::color(*this, seed);
    }

/*
    Real CalculateFaceMeasureOverCellDist(IndexType edgeIndex, MeshDataContainer){

        const auto& edge = this->getEdges().at(edgeIndex);
        return CalculateEdgeMeasure(edgeIndex) / CalculateCellDist(edge.GetCellLeftIndex(), edge.GetCellRightIndex());

    }
  */
    // do not compile this if it is not used
    // it may cause problems for higher dimensions
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

    // there is only one implemented reader for mesh of dimension 2
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







    /**
      Jak udělat metodu, která by byla schopná prosházet elementy sítě v rozmezí nějakých dimenzí (obarvení grafu) (samostatná specializovaná třída volaná metodou)
      Jak veřejnit typy přístupné pomocí proměnné (using na typ)
      Jak zavést některé metody, které budou obecné a nemuset je přepisovat do specializací (asi to tam prostě okopírovat)
      Jak udělat procházení podelementů šablonově v závislosti na dimenzi elementu (výpočet některých vlastností sítě)
      Jak například vyrobit pomocné datové struktury jako centra elementů (asi šablonovým děděním)
        (mohl bych navrhnout třídu, která by vytvářela na míru vektory dlouhé jako dané vektory meshelement)


      Jak dobře koncipovat metody (počítání objemu...), zmínit i wrappery

      */



#endif // UNSTRUCTUREDMESH_H

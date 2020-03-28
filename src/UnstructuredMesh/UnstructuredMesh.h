#ifndef UNSTRUCTUREDMESH_H
#define UNSTRUCTUREDMESH_H
#include "MeshElements/MeshElements.h"
#include "vector"
#include <type_traits>
#include "MeshFunctions/MeshFunctions.h"
#include "MeshFunctions/ComputeCenters.h"
#include "MeshFunctions/ComputeMeasures.h"


template <unsigned int Dimension, typename IndexType, typename Real, unsigned int ...Reserve>
class UnstructuredMesh : public MeshElements<Dimension, IndexType, Real, Reserve...>{

public:
    template<ComputationMethod Method = ComputationMethod::DEFAULT>
    void initializeCenters(){
        auto centers = computeCenters<Method>(*this);

        for (auto& face : this->getFaces()){
            face.setCenter(centers[face]);
        }
        for (auto& cell : this->getCells()){
            cell.setCenter(centers[cell]);
        }
    }

    template<ComputationMethod Method = ComputationMethod::DEFAULT>
    MakeMeshDataContainer_t<Real, make_custom_integer_sequence_t<unsigned int, 1, Dimension>> computeElementMeasures() {
        return computeMeasures<Method>(*this);
    }

    template<ComputationMethod Method = ComputationMethod::DEFAULT>
    MeshDataContainer<Vector<Dimension, Real>, Dimension-1> computeFaceNormals() {
        return ::ComputeFaceNormals<Method>(*this);
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
public:
    Real computeCellDist(IndexType cellIndex1, IndexType cellIndex2){
        return (this->getCells().at(cellIndex1).GetCenter() - this->getCells().at(cellIndex2).GetCenter()).NormEukleid();
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

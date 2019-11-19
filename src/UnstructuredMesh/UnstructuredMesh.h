#ifndef UNSTRUCTUREDMESH_H
#define UNSTRUCTUREDMESH_H
#include "MeshElements/MeshElement.h"
#include "vector"
#include <type_traits>
#include "MeshFunctions/MeshFunctions.h"


template <unsigned int Dimension, typename IndexType, typename Real, unsigned int ...Reserve>
class UnstructuredMesh : public MeshElements<Dimension, IndexType, Real, Reserve...>{

public:
    void initializeCenters(){
        auto centers = ComputeCenters(*this);

        for (auto& face : this->getFaces()){
            face.setCenter(centers[face]);
        }
        for (auto& cell : this->getCells()){
            cell.setCenter(centers[cell]);
        }
    }

    MakeMeshDataContainer_t<Real, make_custom_integer_sequence_t<unsigned int, 1, Dimension>> computeElementMeasures() {
        return ComputeMeasures(*this);
    }

    MeshDataContainer<Vector<Dimension, Real>, Dimension-1> computeFaceNormals() {
        return ::ComputeFaceNormals(*this);
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
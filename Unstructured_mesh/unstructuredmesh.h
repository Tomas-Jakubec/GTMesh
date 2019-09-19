#ifndef UNSTRUCTUREDMESH_H
#define UNSTRUCTUREDMESH_H
#include "mesh_element.h"
#include "vector"
#include <type_traits>



template <unsigned int Dimension, typename IndexType, typename Real, unsigned int ...Reserve>
class UnstructuredMesh : public MeshElements<Dimension, IndexType, Real, Reserve...>{


public:
    Real CalculateCellDist(IndexType cellIndex1, IndexType cellIndex2){
        return (this->GetCells().at(cellIndex1).GetCenter() - this->GetCells().at(cellIndex2).GetCenter()).NormEukleid();
    }
};

template <typename IndexType, typename Real, unsigned int ...Reserve>
class UnstructuredMesh<2, IndexType, Real, Reserve...> : public MeshElements<2, IndexType, Real, Reserve...>{

public:


    Real CalculateCellDist(IndexType cellIndex1, IndexType cellIndex2){
        return (this->GetCells().at(cellIndex1).GetCenter() - this->GetCells().at(cellIndex2).GetCenter()).NormEukleid();
    }



    Real CalculateEdgeMeasure(IndexType index){
        {
            auto& edge = this->GetEdges().at(index);
            return (this->GetVertices().at(edge.GetVertexAIndex()) -
                    this->GetVertices().at(edge.GetVertexBIndex())).NormEukleid();
        }
    }

    Real CalculateFaceMeasure(IndexType index){
        {
            auto& edge = this->GetEdges().at(index);
            return (this->GetVertices().at(edge.GetVertexAIndex()) -
                    this->GetVertices().at(edge.GetVertexBIndex())).NormEukleid();
        }
    }

    Real CalculateCellMeasure(IndexType index){
        auto& cell = this->GetCells().at(index);
        IndexType tmpEdgeIndex = cell.GetBoundaryElementIndex();

        Vertex<2,Real> c = cell.GetCenter();

        double volume = 0;

        do {

            Vertex<2,Real> a = this->GetVertices().at(this->GetEdges().at(tmpEdgeIndex).GetVertexAIndex());
            Vertex<2,Real> b = this->GetVertices().at(this->GetEdges().at(tmpEdgeIndex).GetVertexBIndex());
            Real tmp = (c[0] - a[0]) * (b[1] - a[1]);
            tmp -= (c[1] - a[1]) * (b[0] - a[0]);
            volume += 0.5 * abs(tmp);
            tmpEdgeIndex = this->GetEdges().at(tmpEdgeIndex).GetNextBElem(cell.GetIndex());

        } while (tmpEdgeIndex != cell.GetBoundaryElementIndex());

        return volume;
    }
/*
template <typename data_t, unsigned int ...dimensions>
    struct data : public data<lowerdim, higherdim -1, data_t>{
        std::vector<data_t> vec;
    };
    template <unsigned int lowerdim, typename data_t>
        struct data<lowerdim, lowerdim, data_t>{
            std::vector<data_t> vec;
        };

*/

    void InitializeCenters(){
        auto& vertices = this->GetVertices();
        for(auto& edge : this->GetEdges()){
            edge.SetCenter((vertices.at(edge.GetVertexAIndex()) + vertices.at(edge.GetVertexBIndex())) * 0.5);
        }
        for(auto& cell : this->GetCells()){
            IndexType tmpEdge = cell.GetBoundaryElementIndex();
            Vertex<2, Real> tmpCenter = {0.0, 0.0};
            double counter = 0.0;
            do {
                tmpCenter += this->GetEdges().at(tmpEdge).GetCenter();
                counter++;
                tmpEdge = this->GetEdges().at(tmpEdge).GetNextBElem(cell.GetIndex());
            } while (tmpEdge != cell.GetBoundaryElementIndex());
            cell.SetCenter(tmpCenter / counter);
            DBGVAR(tmpCenter * (1.0 / counter));
        }
    }

    Real CalculateFaceMeasureOverCellDist(IndexType edgeIndex){

        const auto& edge = this->GetEdges().at(edgeIndex);
        return CalculateEdgeMeasure(edgeIndex) / CalculateCellDist(edge.GetCellLeftIndex(), edge.GetCellRightIndex());

    }

    /**
      Jak udělat metodu, která by byla schopná prosházet elementy sítě v rozmezí nějakých dimenzí (obarvení grafu) (samostatná specializovaná třída volaná metodou)
      Jak veřejnit typy přístupné pomocí proměnné (using na typ)
      Jak zavést některé metody, které budou obecné a nemuset je přepisovat do specializací (asi to tam prostě okopírovat)
      Jak udělat procházení podelementů šablonově v závislosti na dimenzi elementu (výpočet některých vlastností sítě)
      Jak například vyrobit pomocné datové struktury jako centra elementů (asi šablonovým děděním)
        (mohl bych navrhnout třídu, která by vytvářela na míru vektory dlouhé jako dané vektory meshelement)


      Jak dobře koncipovat metody (počítání objemu...), zmínit i wrappery

      */
};
/*
class CellWrap{
    UnstructuredMesh<2,size_t, double> *m;
    size_t cellIndex;

};
*/
#endif // UNSTRUCTUREDMESH_H
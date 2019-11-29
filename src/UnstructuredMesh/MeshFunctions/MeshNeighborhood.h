#ifndef MESHNEIGHBORHOOD_H
#define MESHNEIGHBORHOOD_H

#include "MeshConnections.h"

template <unsigned int StartDim, unsigned int ConnectingDim, unsigned int ConnectedDim = StartDim, Order order = Order::ORDER_ASCEND>
class MeshNegborhood{
    template<unsigned int MeshDimension, typename IndexType, typename Real, unsigned int ...Reserve>
    static MeshDataContainer<std::set<IndexType>, StartDim> neighbors(
                MeshElements<MeshDimension, IndexType, Real, Reserve...>& mesh
            ) {

        MeshDataContainer<std::set<IndexType>, StartDim> result;
        MeshDataContainer<std::set<IndexType>, StartDim> firstConnections = MeshConnections<StartDim, ConnectingDim, order>::connections(mesh);
        MeshDataContainer<std::set<IndexType>, StartDim> secondConnections = MeshConnections<ConnectingDim, ConnectedDim, order>::connections(mesh);

        for (IndexType elementIndex = 0; elementIndex < mesh.template getElements<StartDim>().size(); elementIndex++) {

            for (IndexType& firstConectedElem : firstConnections.template getDataByPos<0>().at(elementIndex)){

                for (IndexType& neighborIndex : secondConnections.template getDataByPos<0>().at(firstConectedElem)){

                    if (StartDim == ConnectedDim && elementIndex == neighborIndex) {
                        continue;
                    } else {
                        result.template getDataByPos<0>().at(elementIndex).insert(neighborIndex);
                    }
                }
            }
        }

        return result;
    }

};


template <unsigned int StartDim, unsigned int ConnectingDim, unsigned int ConnectedDim>
class MeshNegborhood<StartDim, ConnectingDim, ConnectedDim, Order::ORDER_ORIGINAL>{
    template<unsigned int MeshDimension, typename IndexType, typename Real, unsigned int ...Reserve>
    static MeshDataContainer<std::set<IndexType>, StartDim> neighbors(
                MeshElements<MeshDimension, IndexType, Real, Reserve...>& mesh
            ) {

        MeshDataContainer<std::map<IndexType, IndexType>, StartDim> tempResultMap(mesh);
        MeshDataContainer<std::vector<IndexType>, StartDim> firstConnections = MeshConnections<StartDim, ConnectingDim, Order::ORDER_ORIGINAL>::connections(mesh);
        MeshDataContainer<std::vector<IndexType>, StartDim> secondConnections = MeshConnections<ConnectingDim, ConnectedDim, Order::ORDER_ORIGINAL>::connections(mesh);

        for (IndexType elementIndex = 0; elementIndex < mesh.template getElements<StartDim>().size(); elementIndex++) {

            for (IndexType& firstConectedElem : firstConnections.template getDataByPos<0>().at(elementIndex)){

                for (IndexType& neighborIndex : secondConnections.template getDataByPos<0>().at(firstConectedElem)){

                    if (StartDim == ConnectedDim && elementIndex == neighborIndex) {
                        continue;
                    } else {
                        IndexType pos = tempResultMap.template getDataByPos<0>().at(elementIndex).size();
                        tempResultMap.template getDataByPos<0>().at(elementIndex).insert({neighborIndex, pos});
                    }
                }
            }
        }

        MeshDataContainer<std::vector<IndexType>, StartDim> result;

        for (IndexType i = 0; i < mesh.template getElements<StartDim>().size(); i++){
            //resize the vector at the position
            result.template getDataByPos<0>().at(i).resize(
                tempResultMap.template getDataByPos<0>().at(i).size()
            );

            for(std::pair<IndexType, IndexType>& mapElem : tempResultMap) {
                result.template getDataByPos<0>().at(i).at(mapElem.second) = mapElem.first;
            }
        }
        return result;
    }

};


#endif // MESHNEIGHBORHOOD_H

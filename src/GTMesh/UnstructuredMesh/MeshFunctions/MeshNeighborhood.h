#ifndef MESHNEIGHBORHOOD_H
#define MESHNEIGHBORHOOD_H

#include "MeshConnections.h"

template <unsigned int StartDim, unsigned int ConnectingDim, unsigned int ConnectedDim = StartDim, Order order = Order::ORDER_ASCEND>
class MeshNeighborhood{
public:
    template<unsigned int MeshDimension, typename IndexType, typename Real, unsigned int ...Reserve>
    static MeshDataContainer<std::vector<IndexType>, StartDim> neighbors(
                const MeshElements<MeshDimension, IndexType, Real, Reserve...>& mesh
            ) {

        MeshDataContainer<std::vector<IndexType>, StartDim> result(mesh);
        auto firstConnections = MeshConnections<StartDim, ConnectingDim, order>::connections(mesh);
        auto secondConnections = MeshConnections<ConnectingDim, ConnectedDim, order>::connections(mesh);

        for (IndexType elementIndex = 0; elementIndex < mesh.template getElements<StartDim>().size(); elementIndex++) {

            std::set<IndexType> tmpResultSet;

            for (IndexType& firstConectedElem : firstConnections.template getDataByPos<0>().at(elementIndex)){

                for (IndexType& neighborIndex : secondConnections.template getDataByPos<0>().at(firstConectedElem)){

                    if (StartDim == ConnectedDim && elementIndex == neighborIndex) {
                        continue;
                    } else {
                        tmpResultSet.insert(neighborIndex);
                    }
                }
            }
            result.template getDataByPos<0>()[elementIndex].insert(
                        result.template getDataByPos<0>()[elementIndex].begin(),
                        tmpResultSet.begin(),
                        tmpResultSet.end()
                        );
        }

        return result;
    }

};


template <unsigned int StartDim, unsigned int ConnectingDim, unsigned int ConnectedDim>
class MeshNeighborhood<StartDim, ConnectingDim, ConnectedDim, Order::ORDER_ORIGINAL>{
public:
    template<unsigned int MeshDimension, typename IndexType, typename Real, unsigned int ...Reserve>
    static MeshDataContainer<std::vector<IndexType>, StartDim> neighbors(
                const MeshElements<MeshDimension, IndexType, Real, Reserve...>& mesh
            ) {

        MeshDataContainer<std::vector<IndexType>, StartDim> result(mesh);
        auto firstConnections = MeshConnections<StartDim, ConnectingDim, Order::ORDER_ORIGINAL>::connections(mesh);
        auto secondConnections = MeshConnections<ConnectingDim, ConnectedDim, Order::ORDER_ORIGINAL>::connections(mesh);

        for (IndexType elementIndex = 0; elementIndex < mesh.template getElements<StartDim>().size(); elementIndex++) {

            std::map<IndexType, IndexType> tempResultMap;

            for (IndexType& firstConectedElem : firstConnections.template getDataByPos<0>().at(elementIndex)){

                for (IndexType& neighborIndex : secondConnections.template getDataByPos<0>().at(firstConectedElem)){

                    if (StartDim == ConnectedDim && elementIndex == neighborIndex) {
                        continue;
                    } else {
                        IndexType pos = tempResultMap.size();
                        tempResultMap.insert({neighborIndex, pos});
                    }
                }
            }

            result.template getDataByPos<0>().at(elementIndex).resize(
                tempResultMap.size()
            );

            for(std::pair<const IndexType, IndexType>& mapElem : tempResultMap) {
                result.template getDataByPos<0>().at(elementIndex).at(mapElem.second) = mapElem.first;
            }
        }

        return result;
    }

};


#endif // MESHNEIGHBORHOOD_H

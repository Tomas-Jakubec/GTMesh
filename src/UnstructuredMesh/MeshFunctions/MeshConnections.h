#ifndef MESHCONNECTIONS_H
#define MESHCONNECTIONS_H

#include "../MeshElements/MeshElement.h"
#include "../MeshDataContainer/MeshDataContainer.h"
#include "../../NumericStaticArray/Vector.h"
#include "../../Debug/Debug.h"
#include "MeshApply.h"
#include <valarray>
#include <set>
#include <map>


enum Order{
    ORDER_ASCEND,
    ORDER_ORIGINAL
};

template<unsigned int StartDim, unsigned int TargetDim, Order order = Order::ORDER_ASCEND>
struct MeshConnections {
    /**
     * @brief connections<HR>
     * Detects connections of mesh elements of StartDim to TargetDim.
     * Returns a MeshDataContainer of set<IndexType> allocated to StartDim elements.
     * The indexes are ordered in ascending way.
     * @param mesh
     * @return
     */
    template<unsigned int MeshDimension, typename IndexType, typename Real, unsigned int ...Reserve>
    static MeshDataContainer<std::set<IndexType>, StartDim> connections(
            MeshElements<MeshDimension, IndexType, Real, Reserve...>& mesh
            ) {
        MeshDataContainer<std::set<IndexType>, StartDim> result(mesh);
        MeshApply<StartDim, TargetDim, MeshDimension>::apply(mesh, [&result](IndexType ori, IndexType element){
            result.template getDataByPos<0>().at(ori).insert(element);
        });

        return result;
    }
};


template<unsigned int StartDim, unsigned int TargetDim>
struct MeshConnections<StartDim, TargetDim, Order::ORDER_ORIGINAL> {

    /**
     * @brief orderedConnections<HR>
     * This function returns connection in original sequence as in the mesh.
     * @param mesh
     * @return
     */
    template<unsigned int MeshDimension, typename IndexType, typename Real, unsigned int ...Reserve>
    static MeshDataContainer<std::vector<IndexType>, StartDim> connections(
            MeshElements<MeshDimension, IndexType, Real, Reserve...>& mesh
            ) {
        MeshDataContainer<std::map<IndexType, IndexType>, StartDim> tempMap(mesh);

        MeshApply<StartDim, TargetDim, MeshDimension>::apply(mesh, [&tempMap](IndexType ori, IndexType element){
            IndexType size = tempMap.template getDataByPos<0>().at(ori).size();
            tempMap.template getDataByPos<0>().at(ori).insert({element, size});
        });

        MeshDataContainer<std::vector<IndexType>, StartDim> result(mesh);
        for (IndexType i = 0; i < mesh.template getElements<StartDim>().size(); i++){
            //resize the vector at the position
            result.template getDataByPos<0>().at(i).resize(
                tempMap.template getDataByPos<0>().at(i).size()
            );

            for(std::pair<const IndexType, IndexType>& mapElem : tempMap.template getDataByPos<0>().at(i)) {
                result.template getDataByPos<0>().at(i).at(mapElem.second) = mapElem.first;
            }
        }
        return result;
    }

};

#endif // MESHCONNECTIONS_H

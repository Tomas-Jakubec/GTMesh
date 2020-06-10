#ifndef MESHCONNECTIONS_H
#define MESHCONNECTIONS_H

#include "../MeshElements/MeshElements.h"
#include "../MeshDataContainer/MeshDataContainer.h"
#include "../../NumericStaticArray/Vector.h"
#include "../../Debug/Debug.h"
#include "MeshApply.h"
#include "MeshFunctionsDefine.h"
#include <valarray>
#include <set>
#include <map>




template<unsigned int StartDim, unsigned int TargetDim, Order order = Order::ORDER_ASCEND>
struct MeshConnections {
    /**
     * @brief Detects connections of mesh elements of StartDim to TargetDim.
     * Returns a MeshDataContainer of vector<IndexType> allocated to StartDim elements.
     * The indexes are ordered in ascending way.
     * @param mesh
     * @return
     */
    template<unsigned int MeshDimension, typename IndexType, typename Real, unsigned int ...Reserve>
    static
    typename std::enable_if <
    (StartDim < TargetDim),
    MeshDataContainer<std::vector<IndexType>, StartDim>
    >::type
    connections(
            MeshElements<MeshDimension, IndexType, Real, Reserve...>& mesh
            ) {

        MeshDataContainer<std::set<IndexType>, StartDim> tmpSet(mesh);

        // Resolve the connections
        MeshApply<StartDim, TargetDim>::apply(mesh, [&tmpSet](IndexType orig, IndexType connected){
            tmpSet.template getDataByPos<0>().at(orig).insert(connected);
        });
        MeshDataContainer<std::vector<IndexType>, StartDim> res(mesh);
        // Copy the temporary data into result data
        for (IndexType i = IndexType(); i < res.template getDataByPos<0>().size(); i++) {
            res.template getDataByPos<0>()[i].insert(
                        res.template getDataByPos<0>()[i].begin(),
                        tmpSet.template getDataByPos<0>()[i].begin(),
                        tmpSet.template getDataByPos<0>()[i].end());
        }
        return res;
    }

    /**
     * @brief connections<HR>
     * Detects connections of mesh elements of StartDim to TargetDim.
     * Returns a MeshDataContainer of vector<IndexType> allocated to StartDim elements.
     * The indexes are ordered in ascending way.
     * @param mesh
     * @return
     */
    template<unsigned int MeshDimension, typename IndexType, typename Real, unsigned int ...Reserve>
    static
    typename std::enable_if <
    (StartDim >= TargetDim),
    MeshDataContainer<std::vector<IndexType>, StartDim>
    >::type
    connections(
            MeshElements<MeshDimension, IndexType, Real, Reserve...>& mesh
            ) {

        std::set<IndexType> tmpSet;
        MeshDataContainer<std::vector<IndexType>, StartDim> res(mesh);

        for(IndexType index = IndexType(); index < mesh.template getElements<StartDim>().size(); index++) {

            MeshApply<StartDim, TargetDim>::apply(index, mesh, [&tmpSet](IndexType, IndexType element){
                tmpSet.insert(element);
            });

            res.template getDataByPos<0>()[index].insert(
                        res.template getDataByPos<0>()[index].begin(),
                        tmpSet.begin(),
                        tmpSet.end());

            tmpSet.clear();
        }

        return res;
    }


    /**
     * @brief connections<HR>
     * Detects connections of mesh elements of StartDim to TargetDim.
     * Returns a MeshDataContainer of vector<IndexType> allocated to StartDim elements.
     * The indexes are ordered in ascending way.
     * @param mesh
     * @return
     */
    // template<unsigned int MeshDimension, typename IndexType, typename Real, unsigned int ...Reserve>
    // static
    // typename std::enable_if <
    // (StartDim == MeshDimension - 1) && (TargetDim == MeshDimension),
    // MeshDataContainer<std::array<IndexType, 2>, StartDim>
    // >::type
    // connections(
    //         MeshElements<MeshDimension, IndexType, Real, Reserve...>& mesh
    //         ) {
    //
    //
    //     MeshDataContainer<std::array<IndexType, 2>, StartDim> res(mesh);
    //
    //     for(IndexType index = IndexType(); index < mesh.template getElements<StartDim>().size(); index++) {
    //
    //         auto& face = mesh.template getElements<StartDim>()[index];
    //
    //         res.template getDataByPos<0>()[index] = std::array<IndexType, 2>{{face.getCellLeftIndex(), face.getCellRightIndex()}};
    //
    //
    //
    //     }
    //
    //     return res;
    // }
};


template<unsigned int StartDim, unsigned int TargetDim>
struct MeshConnections<StartDim, TargetDim, Order::ORDER_ORIGINAL> {

    /**
     * @brief This function returns connection in original sequence as in the mesh.
     * @param mesh
     * @return
     */
    template<unsigned int MeshDimension, typename IndexType, typename Real, unsigned int ...Reserve>
    static MeshDataContainer<std::vector<IndexType>, StartDim> connections(
            MeshElements<MeshDimension, IndexType, Real, Reserve...>& mesh
            ) {
        MeshDataContainer<std::map<IndexType, IndexType>, StartDim> tempMap(mesh);

        MeshApply<StartDim, TargetDim>::apply(
            mesh,
            [&tempMap](IndexType orig, IndexType connected){
                IndexType size = tempMap.template getDataByPos<0>().at(orig).size();
                tempMap.template getDataByPos<0>().at(orig).insert({connected, size});
            }
        );

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

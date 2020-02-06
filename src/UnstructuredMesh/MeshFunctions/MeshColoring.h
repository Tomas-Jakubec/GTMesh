#ifndef MESHCOLOURING_H
#define MESHCOLOURING_H

#include "MeshConnections.h"

template<unsigned int FromDim, unsigned int ToDim, bool Descend = true>
struct MeshColoring {

    template<unsigned int MeshDimension, typename IndexType, typename Real, unsigned int ...Reserve>
    static MeshDataContainer<unsigned int, FromDim> color(
            MeshElements<MeshDimension, IndexType, Real, Reserve...>& mesh
            ) {
        MeshDataContainer<unsigned int, FromDim> result(mesh);

        DBGMSG("starting the coloring procedure");
        unsigned int reserve = 16;
        MeshDataContainer<std::valarray<bool>, ToDim> attachedColours(mesh, std::valarray<bool>(false, reserve));



        for (IndexType elementIndex = 0; elementIndex < mesh.template getElements<FromDim>().size(); elementIndex++){
            std::valarray<bool> possibleColours(true,reserve);
            MeshApply<FromDim, ToDim>::apply(elementIndex, mesh,
                 [&possibleColours, &attachedColours](IndexType, IndexType element){
                     DBGTRY(possibleColours &= !attachedColours.template getDataByPos<0>().at(element);)
            });

            // Select the first possible colour
            unsigned int selectedColour = 0;
            while (!possibleColours[selectedColour]) {
                selectedColour++;
                if (selectedColour == possibleColours.size()){
                    reserve *= 2;
                    DBGVAR(reserve);
                    for (std::valarray<bool>& attColour : attachedColours.template getDataByPos<0>()){
                        std::valarray<bool> newAttColour(false, reserve);
                        for (size_t i = 0; i < attColour.size(); i++){
                            newAttColour[i] = attColour[i];
                        }
                        attColour.swap(newAttColour);
                    }
                    break;
                }
            }
            result.template getDataByPos<0>().at(elementIndex) = selectedColour;
            MeshApply<FromDim, ToDim>::apply(elementIndex, mesh,
                 [selectedColour, &attachedColours](IndexType, IndexType element){
                     DBGTRY(attachedColours.template getDataByPos<0>().at(element)[selectedColour] = true;)
            });

        }
        return result;
    }
};


template<unsigned int FromDim, unsigned int ToDim>
struct MeshColoring <FromDim, ToDim, false> {
    template<unsigned int MeshDimension, typename IndexType, typename Real, unsigned int ...Reserve>
    static MeshDataContainer<unsigned int, FromDim> color(
            MeshElements<MeshDimension, IndexType, Real, Reserve...>& mesh
            ) {
        // resulting container of colours
        MeshDataContainer<unsigned int, FromDim> result(mesh);

        DBGMSG("starting the coloring procedure");
        unsigned int reserve = 16;
        // allocates memory to the given dimension
        MeshDataContainer<std::valarray<bool>, ToDim> attachedColours(mesh, std::valarray<bool>(false, reserve));

        auto connections = MeshConnections<FromDim, ToDim>::connections(mesh);

        for (IndexType elementIndex = 0; elementIndex < mesh.template getElements<FromDim>().size(); elementIndex++){
            std::valarray<bool> possibleColours(true,reserve);
            for (IndexType element : connections.template getDataByPos<0>().at(elementIndex)){

                DBGTRY(possibleColours &= !attachedColours.template getDataByPos<0>().at(element);)

            }

            // Select the first possible colour
            unsigned int selectedColour = 0;
            while (!possibleColours[selectedColour]) {
                selectedColour++;
                if (selectedColour == possibleColours.size()){
                    reserve *= 2;

                    // If the number of colours exceeds the number of
                    // allocated bits, then allocate twice as much memory
                    for (std::valarray<bool>& attColour : attachedColours.template getDataByPos<0>()){
                        std::valarray<bool> newAttColour(false, reserve);
                        for (size_t i = 0; i < attColour.size(); i++){
                            newAttColour[i] = attColour[i];
                        }
                        attColour.swap(newAttColour);
                    }
                    break;
                }
            }

            result.template getDataByPos<0>().at(elementIndex) = selectedColour;

            for (IndexType element : connections.template getDataByPos<0>().at(elementIndex)){
                DBGTRY(attachedColours.template getDataByPos<0>().at(element)[selectedColour] = true;)
            }
        }
        return result;
    }
};


template <unsigned int FromDim, unsigned int ToDim>
struct ColorMesh{

template<unsigned int MeshDimension, typename IndexType, typename Real, unsigned int ...Reserve>
static MeshDataContainer<unsigned int, FromDim> color(
            MeshElements<MeshDimension, IndexType, Real, Reserve...>& mesh
        ){
    return MeshColoring<FromDim, ToDim, (FromDim > ToDim)>::color(mesh);
}
};


#endif // MESHCOLOURING_H

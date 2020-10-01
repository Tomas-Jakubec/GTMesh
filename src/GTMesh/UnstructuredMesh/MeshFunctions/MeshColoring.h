#ifndef MESHCOLOURING_H
#define MESHCOLOURING_H

#include "MeshConnections.h"
#include <random>

enum ColoringMethod{
    METHOD_GREEDY,
    METHOD_RANDOM
};
namespace Impl {


template<unsigned int ColoredDim, unsigned int ConnectingDim, ColoringMethod Method = METHOD_GREEDY, bool Descend = (ColoredDim > ConnectingDim)>
struct _MeshColoring {
    template<unsigned int MeshDimension, typename IndexType, typename Real, unsigned int ...Reserve>
    static MeshDataContainer<unsigned int, ColoredDim> color(
            const MeshElements<MeshDimension, IndexType, Real, Reserve...>& mesh
            ) {
        // resulting container of colours
        MeshDataContainer<unsigned int, ColoredDim> result(mesh);
        // Setup the initial reserve
        unsigned int reserve = 16;
        // allocates memory to the given dimension
        MeshDataContainer<std::valarray<bool>, ConnectingDim> attachedColors(mesh, std::valarray<bool>(false, reserve));

        auto connections = MeshConnections<ColoredDim, ConnectingDim>::connections(mesh);

        for (IndexType elementIndex = 0; elementIndex < mesh.template getElements<ColoredDim>().size(); elementIndex++){
            std::valarray<bool> freeColors(true,reserve);
            for (IndexType element : connections.template getDataByPos<0>().at(elementIndex)){

                DBGTRY(freeColors &= !attachedColors.template getDataByPos<0>().at(element);)

            }

            // Select the first possible colour
            unsigned int selectedColor = 0;
            while (!freeColors[selectedColor]) {
                selectedColor++;

                // If the number of colors exceeds the number of
                // allocated bits, then allocate twice as much memory
                if (selectedColor == freeColors.size()){
                    reserve *= 2;
                    for (std::valarray<bool>& attColor : attachedColors.template getDataByPos<0>()){
                        std::valarray<bool> newAttColor(false, reserve);
                        for (size_t i = 0; i < attColor.size(); i++){
                            newAttColor[i] = attColor[i];
                        }
                        attColor.swap(newAttColor);
                    }
                    break;
                }
            }

            result.template getDataByPos<0>().at(elementIndex) = selectedColor;

            for (IndexType element : connections.template getDataByPos<0>().at(elementIndex)){
                DBGTRY(attachedColors.template getDataByPos<0>().at(element)[selectedColor] = true;)
            }
        }
        return result;
    }

};


template<unsigned int ColoredDim, unsigned int ConnectingDim>
struct _MeshColoring <ColoredDim, ConnectingDim, METHOD_GREEDY, true> {
    template<unsigned int MeshDimension, typename IndexType, typename Real, unsigned int ...Reserve>
    static MeshDataContainer<unsigned int, ColoredDim> color(
            const MeshElements<MeshDimension, IndexType, Real, Reserve...>& mesh
            ) {
        MeshDataContainer<unsigned int, ColoredDim> result(mesh);

        DBGMSG("starting the coloring procedure");
        unsigned int reserve = 16;
        MeshDataContainer<std::valarray<bool>, ConnectingDim> attachedColors(mesh, std::valarray<bool>(false, reserve));



        for (IndexType elementIndex = 0; elementIndex < mesh.template getElements<ColoredDim>().size(); elementIndex++){
            std::valarray<bool> freeColors(true,reserve);
            MeshApply<ColoredDim, ConnectingDim>::apply(elementIndex, mesh,
                 [&freeColors, &attachedColors](IndexType, IndexType element){
                     DBGTRY(freeColors &= !attachedColors.template getDataByPos<0>().at(element);)
            });

            // Select the first possible colour
            unsigned int selectedColor = 0;
            while (!freeColors[selectedColor]) {
                selectedColor++;
                if (selectedColor == freeColors.size()){
                    reserve *= 2;
                    DBGVAR(reserve);
                    for (std::valarray<bool>& attColour : attachedColors.template getDataByPos<0>()){
                        std::valarray<bool> newAttColour(false, reserve);
                        for (size_t i = 0; i < attColour.size(); i++){
                            newAttColour[i] = attColour[i];
                        }
                        attColour.swap(newAttColour);
                    }
                    break;
                }
            }
            result.template getDataByPos<0>().at(elementIndex) = selectedColor;
            MeshApply<ColoredDim, ConnectingDim>::apply(elementIndex, mesh,
                 [selectedColor, &attachedColors](IndexType, IndexType element){
                     DBGTRY(attachedColors.template getDataByPos<0>().at(element)[selectedColor] = true;)
            });

        }
        return result;
    }
};


template<unsigned int ColoredDim, unsigned int ConnectingDim, bool Descend>
struct _MeshColoring <ColoredDim, ConnectingDim, METHOD_RANDOM, Descend> {
    template<unsigned int MeshDimension, typename IndexType, typename Real, unsigned int ...Reserve>
    static MeshDataContainer<unsigned int, ColoredDim> color(
            const MeshElements<MeshDimension, IndexType, Real, Reserve...>& mesh,
            unsigned int seed = 1562315
            ) {
        // calculating the reserve
        auto result = _MeshColoring<ColoredDim, ConnectingDim, METHOD_GREEDY>::color(mesh);

        // get maximal color index
        auto max = std::max_element(result.template getDataByPos<0>().begin(), result.template getDataByPos<0>().end());

        unsigned int reserve = *max + 1;
        // setup the random machine
        std::mt19937 generator(seed);

        DBGMSG("starting the recoloring procedure");
        // allocates memory to the given dimension
        MeshDataContainer<std::valarray<bool>, ConnectingDim> attachedColors(mesh, std::valarray<bool>(false, reserve));

        auto connections = MeshConnections<ColoredDim, ConnectingDim>::connections(mesh);

        for (IndexType elementIndex = 0; elementIndex < mesh.template getElements<ColoredDim>().size(); elementIndex++){

            for (IndexType element : connections.template getDataByPos<0>().at(elementIndex)){
                DBGTRY(attachedColors.template getDataByPos<0>().at(element)[result.template getDataByPos<0>()[elementIndex]] = true;)
            }
        }

        for (IndexType elementIndex = 0; elementIndex < mesh.template getElements<ColoredDim>().size(); elementIndex++){
            std::valarray<bool> freeColors(true,reserve);
            for (IndexType element : connections.template getDataByPos<0>().at(elementIndex)){

                DBGTRY(freeColors &= !attachedColors.template getDataByPos<0>().at(element);)

            }

            unsigned int nPossibleColors = 0;
            for (bool& free : freeColors){
                nPossibleColors+=static_cast<unsigned int>(free);
            }
            if (nPossibleColors > 0){
                // Select the color randomly
                unsigned int selectedColor = 0;
                std::uniform_int_distribution<unsigned int> distribution(0, nPossibleColors);

                unsigned int rnd = distribution(generator);
                // if the random index of new color is equal to
                // the number of possible colors,
                // then the color of the element remains the same
                if (rnd == nPossibleColors) continue;

                while (rnd != 0 || !freeColors[selectedColor]) {
                    if(freeColors[selectedColor]){
                        rnd--;
                    }
                    selectedColor++;
                }

                result.template getDataByPos<0>().at(elementIndex) = selectedColor;

                for (IndexType element : connections.template getDataByPos<0>().at(elementIndex)){
                    DBGTRY(attachedColors.template getDataByPos<0>().at(element)[result.template getDataByPos<0>()[elementIndex]] = false;)
                    DBGTRY(attachedColors.template getDataByPos<0>().at(element)[selectedColor] = true;)
                }
            }
        }
        return result;
    }
};

}

template <unsigned int ColoredDim, unsigned int ConnectingDim, ColoringMethod Method = METHOD_GREEDY>
struct MeshColoring{

template<unsigned int MeshDimension, typename IndexType, typename Real, unsigned int ...Reserve>
static MeshDataContainer<unsigned int, ColoredDim> color(
            const MeshElements<MeshDimension, IndexType, Real, Reserve...>& mesh
        ){
    return Impl::_MeshColoring<ColoredDim, ConnectingDim, Method>::color(mesh);
}
};


template <unsigned int ColoredDim, unsigned int ConnectingDim>
struct MeshColoring<ColoredDim, ConnectingDim, METHOD_RANDOM>{

template<unsigned int MeshDimension, typename IndexType, typename Real, unsigned int ...Reserve>
static MeshDataContainer<unsigned int, ColoredDim> color(
            const MeshElements<MeshDimension, IndexType, Real, Reserve...>& mesh,
            unsigned int seed = 1562315
        ){
    return Impl::_MeshColoring<ColoredDim, ConnectingDim, METHOD_RANDOM>::color(mesh, seed);
}
};

#endif // MESHCOLOURING_H

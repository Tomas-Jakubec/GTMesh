#ifndef MESHAPPLY_H
#define MESHAPPLY_H

#include "../MeshElements/MeshElements.h"
#include <functional>

namespace Impl {

template <unsigned int CurrentDimension, unsigned int StartDimension, unsigned int TargetDimension, unsigned int MeshDimension, bool End, bool Descend>
struct MeshRun {

    template<typename Functor, typename IndexType, typename Real, unsigned int ...Reserve>
    static void run(const MeshElements<MeshDimension, IndexType, Real, Reserve...>& mesh,
                    IndexType origElementIndex,
                    IndexType subelementIndex,
                    Functor fun){


        const auto& currentElement = mesh.template getElements<CurrentDimension>().at(subelementIndex);
        const auto& subelementContainer = currentElement.getSubelements();
        for (unsigned int i = 0; i < subelementContainer.size(); i++){
            MeshRun< CurrentDimension - 1, StartDimension, TargetDimension, MeshDimension, TargetDimension == CurrentDimension - 1, Descend>
                    ::run(mesh, origElementIndex, subelementContainer[i], fun);

        }


    }
};





template <unsigned int StartDimension, unsigned int TargetDimension, unsigned int MeshDimension, bool Descend>
struct MeshRun<MeshDimension, StartDimension, TargetDimension, MeshDimension, false, Descend> {

    template<typename Functor, typename IndexType, typename Real, unsigned int ...Reserve>
    static void run(const MeshElements<MeshDimension, IndexType, Real, Reserve...>& mesh,
                    IndexType origElementIndex,
                    IndexType subelementIndex,
                    Functor fun){

        const auto& cell = mesh.getCells().at(subelementIndex);
        IndexType tmpFace = cell.getBoundaryElementIndex();
        do {
            MeshRun<MeshDimension - 1, StartDimension, TargetDimension, MeshDimension, TargetDimension == MeshDimension - 1, Descend>
                    ::run(mesh, origElementIndex, tmpFace, fun);
            tmpFace = mesh.getFaces().at(tmpFace).getNextBElem(subelementIndex);
        } while (tmpFace != cell.getBoundaryElementIndex());

    }
};





template <unsigned int StartDimension, unsigned int TargetDimension, unsigned int MeshDimension, bool Descend>
struct MeshRun<1, StartDimension, TargetDimension, MeshDimension, false, Descend> {

    template<typename Functor, typename IndexType, typename Real, unsigned int ...Reserve>
    static void run(const MeshElements<MeshDimension, IndexType, Real, Reserve...>& mesh,
                    IndexType origElementIndex,
                    IndexType subelementIndex,
                    Functor fun){

        const auto& edge = mesh.getEdges().at(subelementIndex);
        MeshRun<0, StartDimension, TargetDimension, MeshDimension, TargetDimension == 0, Descend>
                ::run(mesh, origElementIndex, edge.getVertexAIndex(), fun);
        MeshRun<0, StartDimension, TargetDimension, MeshDimension, TargetDimension == 0, Descend>
                ::run(mesh, origElementIndex, edge.getVertexBIndex(), fun);
    }
};


#define INVALID_FUNCTOR_DESCR \
"The Functor fun must be a function with void return type and two arguments of IndexType,\
 the first is index of StartDimension element and the second is the index of the\
 TargetDimension element"


template <unsigned int CurrentDimension,unsigned int StartDimension, unsigned int TargetDimension, unsigned int MeshDimension>
struct MeshRun<CurrentDimension, StartDimension, TargetDimension, MeshDimension, true, true> {

    template<typename Functor, typename IndexType, typename Real, unsigned int ...Reserve>
    static void run(const MeshElements<MeshDimension, IndexType, Real, Reserve...>& ,
                    IndexType origElementIndex,
                    IndexType subelementIndex,
                    Functor fun){
        static_assert (std::is_assignable<std::function<void(IndexType, IndexType)>,Functor>::value,
                       INVALID_FUNCTOR_DESCR);

            fun(origElementIndex, subelementIndex);

    }
};

template <unsigned int CurrentDimension,unsigned int StartDimension, unsigned int TargetDimension, unsigned int MeshDimension>
struct MeshRun<CurrentDimension, StartDimension, TargetDimension, MeshDimension, true, false> {

    template<typename Functor, typename IndexType, typename Real, unsigned int ...Reserve>
    static void run(const MeshElements<MeshDimension, IndexType, Real, Reserve...>& ,
                    IndexType origElementIndex,
                    IndexType subelementIndex,
                    Functor fun){
        static_assert (std::is_assignable<std::function<void(IndexType, IndexType)>,Functor>::value,
                       INVALID_FUNCTOR_DESCR);

            fun(subelementIndex, origElementIndex);

    }
};
} // namespace Impl


template <unsigned int StartDimension, unsigned int TargetDimension>
struct MeshApply {
    static constexpr unsigned int hDim = (StartDimension > TargetDimension) ? StartDimension : TargetDimension;
    static constexpr unsigned int lDim = (StartDimension > TargetDimension) ? TargetDimension : StartDimension;

    template<unsigned int MeshDimension, typename Functor, typename IndexType, typename Real, unsigned int ...Reserve>
    static void apply(const MeshElements<MeshDimension, IndexType, Real, Reserve...>& mesh,
                      Functor f) {
        for (IndexType currElement = 0; currElement < mesh.template getElements<hDim>().size(); currElement++){
            Impl::MeshRun<
                    hDim,
                    hDim,
                    lDim,
                    MeshDimension,
                    StartDimension == TargetDimension,
                    (StartDimension > TargetDimension)>::run(mesh, currElement, currElement, f);
        }
    }

    template<unsigned int MeshDimension, typename Functor, typename IndexType, typename Real, unsigned int ...Reserve>
    static void apply(IndexType elementIndex,
                      const MeshElements<MeshDimension, IndexType, Real, Reserve...>& mesh,
                      Functor f) {
        static_assert (StartDimension >= TargetDimension, "It is possible \
to iterate over connected elements of a single element \
for StartDimesnion > TargetDimesnion only.");

            Impl::MeshRun<
                    hDim,
                    hDim,
                    lDim,
                    MeshDimension,
                    StartDimension == TargetDimension,
                    (StartDimension > TargetDimension)>::run(mesh, elementIndex, elementIndex, f);

    }
};


#endif // MESHAPPLY_H

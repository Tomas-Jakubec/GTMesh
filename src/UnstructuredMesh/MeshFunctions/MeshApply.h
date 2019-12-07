#ifndef MESHAPPLY_H
#define MESHAPPLY_H

#include "../MeshElements/MeshElement.h"
#include <functional>



template <unsigned int CurrentDimension, unsigned int StartDimension, unsigned int TargetDimension, unsigned int MeshDimension, bool End, bool Descend>
struct MeshRun {

    template<typename Func, typename IndexType, typename Real, unsigned int ...Reserve>
    static void run(const MeshElements<MeshDimension, IndexType, Real, Reserve...>& mesh,
                    IndexType origElementIndex,
                    IndexType index,
                    Func fun){


        auto currentElement = mesh.template getElements<CurrentDimension>().at(index);
        for (auto& sube : currentElement.getSubelements())
        MeshRun< CurrentDimension - 1, StartDimension, TargetDimension, MeshDimension, TargetDimension == CurrentDimension - 1, Descend>::run(mesh, origElementIndex, sube.index, fun);


    }
};





template <unsigned int StartDimension, unsigned int TargetDimension, unsigned int MeshDimension, bool Descend>
struct MeshRun<MeshDimension, StartDimension, TargetDimension, MeshDimension, false, Descend> {

    template<typename Func, typename IndexType, typename Real, unsigned int ...Reserve>
    static void run(const MeshElements<MeshDimension, IndexType, Real, Reserve...>& mesh,
                    IndexType origElementIndex,
                    IndexType index,
                    Func fun){

        auto& cell = mesh.getCells().at(index);
        IndexType tmpFace = cell.getBoundaryElementIndex();
        do {
            MeshRun<MeshDimension - 1, StartDimension, TargetDimension, MeshDimension, TargetDimension == MeshDimension - 1, Descend>::run(mesh, origElementIndex, tmpFace, fun);
            tmpFace = mesh.getFaces().at(tmpFace).getNextBElem(cell.getIndex());
        } while (tmpFace != cell.getBoundaryElementIndex());

    }
};





template <unsigned int StartDimension, unsigned int TargetDimension, unsigned int MeshDimension, bool Descend>
struct MeshRun<1, StartDimension, TargetDimension, MeshDimension, false, Descend> {

    template<typename Func, typename IndexType, typename Real, unsigned int ...Reserve>
    static void run(const MeshElements<MeshDimension, IndexType, Real, Reserve...>& mesh,
                    IndexType origElementIndex,
                    IndexType index,
                    Func fun){

        auto& edge = mesh.getEdges().at(index);
        MeshRun<0, StartDimension, TargetDimension, MeshDimension, TargetDimension == 0, Descend>::run(mesh, origElementIndex, edge.getVertexAIndex(), fun);
        MeshRun<0, StartDimension, TargetDimension, MeshDimension, TargetDimension == 0, Descend>::run(mesh, origElementIndex, edge.getVertexBIndex(), fun);
    }
};





template <unsigned int CurrentDimension,unsigned int StartDimension, unsigned int TargetDimension, unsigned int MeshDimension>
struct MeshRun<CurrentDimension, StartDimension, TargetDimension, MeshDimension, true, true> {

    template<typename Functor, typename IndexType, typename Real, unsigned int ...Reserve>
    static void run(const MeshElements<MeshDimension, IndexType, Real, Reserve...>& ,
                    IndexType origElementIndex,
                    IndexType index,
                    Functor fun){
        static_assert (std::is_assignable<std::function<void(IndexType, IndexType)>,Functor>::value,
                       "The Functor fun must be a function with void return type and two arguments of IndexType,\
 the first is index of StartDimension element and the second is the index of the TargetDimension element");

            fun(origElementIndex, index);

    }
};

template <unsigned int CurrentDimension,unsigned int StartDimension, unsigned int TargetDimension, unsigned int MeshDimension>
struct MeshRun<CurrentDimension, StartDimension, TargetDimension, MeshDimension, true, false> {

    template<typename Functor, typename IndexType, typename Real, unsigned int ...Reserve>
    static void run(const MeshElements<MeshDimension, IndexType, Real, Reserve...>& ,
                    IndexType origElementIndex,
                    IndexType index,
                    Functor fun){
        static_assert (std::is_assignable<std::function<void(IndexType, IndexType)>,Functor>::value,
                       "The Functor fun must be a function with void return type and two arguments of IndexType,\
 the first is index of StartDimension element and the second is the index of the TargetDimension element");

            fun(index, origElementIndex);

    }
};



template <unsigned int StartDimension, unsigned int TargetDimension, unsigned int MeshDimension>
struct MeshApply {
    template<typename Functor, typename IndexType, typename Real, unsigned int ...Reserve>
    static void apply(const MeshElements<MeshDimension, IndexType, Real, Reserve...>& mesh,
                      Functor f) {
        for (auto& startElement : mesh.template getElements<(StartDimension > TargetDimension) ? StartDimension : TargetDimension>()){
            MeshRun<
                    (StartDimension > TargetDimension) ? StartDimension : TargetDimension,
                    (StartDimension > TargetDimension) ? StartDimension : TargetDimension,
                    (StartDimension > TargetDimension) ? TargetDimension : StartDimension,
                    MeshDimension,
                    StartDimension == TargetDimension,
                    (StartDimension > TargetDimension)>::run(mesh, startElement.getIndex(), startElement.getIndex(), f);
        }
    }

    template<typename Functor, typename IndexType, typename Real, unsigned int ...Reserve>
    static void apply(IndexType elementIndex,
                      const MeshElements<MeshDimension, IndexType, Real, Reserve...>& mesh,
                      Functor f) {
            MeshRun<
                    (StartDimension > TargetDimension) ? StartDimension : TargetDimension,
                    (StartDimension > TargetDimension) ? StartDimension : TargetDimension,
                    (StartDimension > TargetDimension) ? TargetDimension : StartDimension,
                    MeshDimension,
                    StartDimension == TargetDimension,
                    (StartDimension > TargetDimension)>::run(mesh, elementIndex, elementIndex, f);

    }
};


#endif // MESHAPPLY_H

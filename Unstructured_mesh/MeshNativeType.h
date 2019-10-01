#ifndef MESHNATIVETYPE_H
#define MESHNATIVETYPE_H
template <unsigned int MeshDimension>
struct MeshNativeType{

};


template<>
struct MeshNativeType<2>{
    enum ElementType{
        LINE = 1,
        TRIANGLE,
        QUAD,
        POLYGON
    };
};

#endif // MESHNATIVETYPE_H

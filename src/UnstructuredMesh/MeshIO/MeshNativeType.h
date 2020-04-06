#ifndef MESHNATIVETYPE_H
#define MESHNATIVETYPE_H
template <unsigned int MeshDimension>
struct MeshNativeType{

};


template<>
struct MeshNativeType<2>{
    enum ElementType{
        LINE = 100,
        TRIANGLE = 200,
        QUAD,
        POLYGON
    };
};

template<>
struct MeshNativeType<3>{
    enum ElementType{
        TETRA = 300,
        HEXAHEDRON,
        WEDGE,
        PYRAMID,
        N_PYRAMID,
        POLYHEDRON
    };
};

#endif // MESHNATIVETYPE_H

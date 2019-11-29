#ifndef VECTOR_H
#define VECTOR_H

#include "Vertex.h"

template <unsigned int Dim, typename Real = double>
using Vector = Vertex<Dim, Real>;


template <typename Real>
Vector<3, Real> vectorProduct(const Vector<3, Real>& v1, const Vector<3, Real>& v2){
    Vector<3,Real> res = {};
    res[0] = v1[1]*v2[2] - v1[2]*v2[1];
    res[1] = v1[2]*v2[0] - v1[0]*v2[2];
    res[2] = v1[0]*v2[1] - v1[1]*v2[0];
    return res;
}



#endif // VECTOR_H

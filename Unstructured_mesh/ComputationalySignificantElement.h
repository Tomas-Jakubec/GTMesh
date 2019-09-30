#ifndef COMPUTATIONALY_SIGNIFICANT_ELEMENT_H
#define COMPUTATIONALY_SIGNIFICANT_ELEMENT_H
#include "Vertex.h"

template <unsigned int MeshDim, typename Real>
class ComputationallySignificantElement
{
protected:
    Vertex<MeshDim, Real> center;
    int flag;
public:
    ComputationallySignificantElement() {
        center = {};
        flag = int();
    }

    Vertex<MeshDim, Real>& getCenter(){
        return center;
    }

    const Vertex<MeshDim, Real>& getCenter() const {
        return center;
    }

    void setCenter(const Vertex<MeshDim, Real>& v) {
        center = v;
    }

    int& getFlag() {
        return flag;
    }

    const int& getFlag() const {
        return flag;
    }
};


#endif // COMPUTATIONALY_SIGNIFICANT_ELEMENT_H

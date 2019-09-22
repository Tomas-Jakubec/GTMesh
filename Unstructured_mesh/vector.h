#ifndef VECTOR_H
#define VECTOR_H

#include <math.h>
#include <iostream>
#include <initializer_list>
#include "inline_array_operations.h"



template <unsigned int Dim, typename Real = double>
class Vector {
    /**
     * @brief Coordinates
     */
    Real Coordinates[Dim] = {};

public:
    Vector(){}
    Vector(std::initializer_list<Real> l){
        *this = l;
    }

    Vector<Dim, Real>& operator =(std::initializer_list<Real> l);

    void SetCoordinate(Real coord, unsigned int pos){
        Coordinates[pos] = coord;
    }

    Real& operator[](unsigned int pos){
        return Coordinates[pos];
    }

    const Real& operator[](unsigned int pos) const {
        return Coordinates[pos];
    }

    Real NormEukleid();

    inline Real SumOfSquares() {
        return  inlineScalarProduct<Dim, Real>::computation(Coordinates, Coordinates);
    }

    Vector<Dim, Real> operator-(const Vector<Dim, Real>&) const;

    Vector<Dim, Real> operator+(const Vector<Dim, Real>&) const;

    Vector<Dim, Real> operator*(const Real&) const;

    Vector<Dim, Real> operator/(const Real&) const;
    /**
     * @brief operator *
     *
     * Scalar product of two vectors
     * @param v
     * @return
     */
    Real operator*(const Vector& v);

    Vector<Dim, Real>& operator+=(const Vector<Dim, Real>&);
    Vector<Dim, Real>& operator-=(const Vector<Dim, Real>&);
    Vector<Dim, Real>& operator*=(const Real&);
    Vector<Dim, Real>& operator/=(const Real&);


    bool operator==(const Vector<Dim, Real>&) const;
    bool operator!=(const Vector<Dim, Real>&) const;
    //friend Vertex<Dim, Real>;
};


template<unsigned int Dim, typename Real>
Vector<Dim, Real>& Vector<Dim, Real>::operator =(std::initializer_list<Real> l){
    unsigned int i = 0;

    for(Real x : l){
        if (i < Dim){
            Coordinates[i] = x;
        }else{
            break;
        }
        i++;
    }
    if (i < Dim){
        for (; i < Dim; i++) {
            Coordinates[i] = Real();
        }
    }
    return *this;
}

template <typename Real>
Vector<3, Real> VectorProduct(const Vector<3, Real>& v1, const Vector<3, Real>& v2){
    Vector<3,Real> res = {};
    res[0] = v1[1]*v2[2] - v1[2]*v2[1];
    res[1] = v1[2]*v2[0] - v1[0]*v2[2];
    res[2] = v1[0]*v2[1] - v1[1]*v2[0];
    return res;
}




/*
** Calculates the Eucleid norm of the point
*/
template <unsigned int Dim, typename Real>
Real Vector<Dim, Real>::NormEukleid(){
    return sqrt(SumOfSquares());
}




/*
** Overloaded operators
*/
//subtracting two points
template <unsigned int Dim, typename Real>
Vector<Dim, Real> Vector<Dim, Real>::operator -(const Vector<Dim, Real>& v) const {
    Vector<Dim, Real> res;
    inlineSubtraction<Dim, Real>::computation(res.Coordinates, this->Coordinates, v.Coordinates);
    return res;
}

//addition of two points
template <unsigned int Dim, typename Real>
Vector<Dim, Real> Vector<Dim, Real>::operator +(const Vector<Dim, Real>& v) const {
    Vector<Dim, Real> res;
    inlineAddition<Dim, Real>::computation(res.Coordinates, this->Coordinates, v.Coordinates);
    return res;
}

//multiplying coordinates with real number
template <unsigned int Dim, typename Real>
Vector<Dim, Real> Vector<Dim, Real>::operator *(const Real& x) const {
    Vector<Dim, Real> res;
    inlineMultiplication<Dim, Real>::computation(res.Coordinates, this->Coordinates, x);
    return res;
}


//division
template <unsigned int Dim, typename Real>
Vector<Dim, Real> Vector<Dim, Real>::operator /(const Real& x) const {
    return this->operator*(Real(1.0)/x);
}




template<unsigned int Dim, typename Real>
Real Vector<Dim, Real>::operator*(const Vector &v)
{
    return inlineScalarProduct<Dim, Real>::computation(Coordinates, v.Coordinates);
}

// Adds value of coordinates of another Point
template <unsigned int Dim, typename Real>
Vector<Dim, Real>& Vector<Dim, Real>::operator +=(const Vector<Dim, Real>& v){
    inlineAddition<Dim, Real>::computation(Coordinates, v.Coordinates);
    return *this;
}

// Subtracts value of coordinates of another Point
template <unsigned int Dim, typename Real>
Vector<Dim, Real>& Vector<Dim, Real>::operator -=(const Vector<Dim, Real>& v){
    inlineSubtraction<Dim, Real>::computation(Coordinates, v.Coordinates);
    return *this;
}


// Adds value of coordinates of another Point
template <unsigned int Dim, typename Real>
Vector<Dim, Real>& Vector<Dim, Real>::operator *=(const Real& x){
    inlineMultiplication<Dim, Real>::computation(Coordinates, x);
    return *this;
}

// Subtracts value of coordinates of another Point
template <unsigned int Dim, typename Real>
Vector<Dim, Real>& Vector<Dim, Real>::operator /=(const Real& x){
    this->operator*=(Real(1.0)/x);
    return *this;
}

// Compares two points wether they are the same
template <unsigned int Dim, typename Real>
bool Vector<Dim, Real>::operator ==(const Vector<Dim, Real>& v) const {
    for(unsigned int i = 0; i < Dim; i++) {
        if(this->operator[](i) != v[i]){
            return false;
        }
    }
    return true;
}

// Compares two points wether they are not the same
template <unsigned int Dim, typename Real>
bool Vector<Dim, Real>::operator !=(const Vector<Dim, Real>& v) const {
    return !(*this == v);
}

template <unsigned int Dim, typename Real = double>
std::ostream& operator <<(std::ostream& ost, const Vector<Dim,Real>& v) {
    for (unsigned int i = 0; i < Dim; i++) {
        ost << v[i] << ' ';
    }

    return ost;
}


template <unsigned int Dim, typename Real = double>
std::istream& operator >>(std::istream& ist, Vector<Dim,Real>& v) {
    for (unsigned int i = 0; i < Dim; i++) {
        ist >> v[i];
    }
    return ist;
}



#endif // VECTOR_H

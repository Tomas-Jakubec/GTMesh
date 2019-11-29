#ifndef GRAMMSCHMIDT_H
#define GRAMMSCHMIDT_H


#include "Vector.h"
#include <array>

/**
 * @brief GrammSchmidt
 * Gramm-Schmidt process without normalization. Can be used to
 * calculate volume in any dimension.
 * @param vectors [in / out] vector containing the vectors and
 * after the process the vectors are changed.
 */
template <unsigned int NumVec, unsigned int Dimension,typename IndexType, typename Real>
void grammSchmidt(std::array<Vector<Dimension, Real>, NumVec>& vectors){
    /*
     * Vector of inverse suquare of norm
     */
    std::array<Real, NumVec> invSumOfSquares;
    invSumOfSquares.at(0) = 1.0 / vectors.at(0).sumOfSquares();
    std::array<Real, NumVec> coef;

    for (IndexType i = 1; i < vectors.size(); i++) {
        /*
         * Coefitiens of scalar products.
         * The coefitients are computed in advance for better stability
         */
        for (IndexType j = 0; j < i; j++) {

            coef.at(j) = (vectors.at(i)*(vectors.at(j))) * invSumOfSquares.at(j);

        }
        for (IndexType j = 0; j < i; j++) {

            vectors.at(i) -= (vectors.at(j) * coef.at(j));
        }

        invSumOfSquares.at(i) = 1.0 / vectors.at(i).sumOfSquares();
    }
}

#endif // GRAMMSCHMIDT_H

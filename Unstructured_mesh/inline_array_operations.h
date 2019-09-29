#ifndef INLINE_ARRAY_OPERATIONS_H
#define INLINE_ARRAY_OPERATIONS_H


template <unsigned int N, typename Real>
struct inlineScalarProduct {
    static inline Real computation(const Real *x, const Real *y){
        return inlineScalarProduct<N-1, Real>::computation(x, y) + x[N-1] * y[N-1];
    }
};


template <typename Real>
struct inlineScalarProduct<1, Real>
{
    static inline double computation(const Real *x, const Real *y){
        return x[0] * y[0];
    }
};




template <unsigned int N, typename Real>
struct inlineAddition{
    static inline void computation(Real *res, const Real *x, const Real *y){
        inlineAddition<N-1, Real>::computation(res, x, y);
        res[N-1] = x[N-1] + y[N-1];
    }
    static inline void computation(Real *x, const Real *y){
        inlineAddition<N-1, Real>::computation(x, y);
        x[N-1] += y[N-1];
    }
};

template <typename Real>
struct inlineAddition<1, Real>{
    static inline void computation(Real *res, const Real *x, const Real *y){
        res[0] = x[0] + y[0];
    }
    static inline void computation(Real *x, const Real *y){
        x[0] += y[0];
    }
};


template <unsigned int N, typename Real>
struct inlineSubtraction{
    static inline void computation(Real *res, const Real *x, const Real *y){
        inlineSubtraction<N-1, Real>::computation(res, x, y);
        res[N-1] = x[N-1] - y[N-1];
    }
    static inline void computation(Real *x, const Real *y){
        inlineSubtraction<N-1, Real>::computation(x, y);
        x[N-1] -= y[N-1];
    }
};

template <typename Real>
struct inlineSubtraction<1, Real>{
    static inline void computation(Real *res, const Real *x, const Real *y){
        res[0] = x[0] - y[0];
    }
    static inline void computation(Real *x, const Real *y){
        x[0] -= y[0];
    }
};



template <unsigned int N, typename Real>
struct inlineMultiplication{
    static inline void computation(Real *res, const Real *x, const Real& alpha){
        inlineMultiplication<N-1, Real>::computation(res, x, alpha);
        res[N-1] = x[N-1] * alpha;
    }
    static inline void computation(Real *x, const Real alpha){
        inlineMultiplication<N-1, Real>::computation(x, alpha);
        x[N-1] *= alpha;
    }
};

template <typename Real>
struct inlineMultiplication<1, Real>{
    static inline void computation(Real *res, const Real *x, const Real& alpha){
        res[0] = x[0] * alpha;
    }
    static inline void computation(Real *x, const Real alpha){
        x[0] *= alpha;
    }
};



#endif // INLINE_ARRAY_OPERATIONS_H

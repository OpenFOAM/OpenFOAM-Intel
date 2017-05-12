#ifndef __Expression_H
#define __Expression_H

namespace Foam
{

namespace Expression
{

template <class T> struct eTraits { typedef typename T::elemType elemType; typedef typename T::boundaryExpr boundaryType; enum Properties { HasMesh=1 }; };

template <class T>
struct eIface
{
    typedef typename eTraits<T>::elemType elemType;
    typedef typename eTraits<T>::boundaryType boundaryType;
    const T& param_;
    eIface(const T& p) : param_(p) {}
    dimensionSet dimensions() const { return param_.dimensions(); }
    elemType get(const label i) const { return param_.get(i); }
    const boundaryType& getBoundary(const label i) const { return param_.getBoundary(i); }
};

template <class T1, int N1, class T2, int N2> struct mesh_finder { static const fvMesh& get(const T1&, const T2&); };
template <class T1, class T2, int N2> struct mesh_finder<T1, 1, T2, N2> { static const fvMesh& get(const T1& t1, const T2&) { return t1.mesh(); } };
template <class T1, int N1, class T2> struct mesh_finder<T1, N1, T2, 1> { static const fvMesh& get(const T1&, const T2& t2) { return t2.mesh(); } };
template <class T1, class T2> struct mesh_finder<T1, 1, T2, 1> { static const fvMesh& get(const T1& t1, const T2&) { return t1.mesh(); } };

template <class T1, class T2>
const fvMesh& getMesh(const T1& t1, const T2& t2)
{
    return mesh_finder<T1, eTraits<T1>::HasMesh, T2, eTraits<T2>::HasMesh>::get(t1, t2);
}

template <class A, class B> struct eReturnTraits { };
template <> struct eReturnTraits<scalar, scalar> { typedef scalar retType; };
template <> struct eReturnTraits<scalar, symmTensor> { typedef symmTensor retType; };



#define REGISTER_TYPE(Type, ElemType) \
\
template <> struct eTraits<Type> { typedef ElemType elemType; typedef fvPatchField<ElemType> boundaryType; enum Properties { HasMesh=1 }; }; \
template <> struct eTraits<fvPatchField<ElemType> > { typedef ElemType elemType; typedef fvPatchField<ElemType> boundaryType; enum Properties { HasMesh=1 }; }; \
template <> struct eTraits<ElemType> { typedef ElemType elemType; typedef ElemType boundaryType; enum Properties { HasMesh=0 }; }; \
\
template <> dimensionSet eIface<Type>::dimensions() const { return param_.dimensions(); } \
template <> typename eIface<Type>::elemType eIface<Type>::get(const label i) const { return param_[i]; } \
template <> const typename eIface<Type>::boundaryType& eIface<Type>::getBoundary(const label i) const { return param_.boundaryField()[i]; } \
\
template <> dimensionSet eIface<fvPatchField<ElemType> >::dimensions() const { return dimless; } \
template <> typename eIface<fvPatchField<ElemType> >::elemType eIface<fvPatchField<ElemType> >::get(const label i) const { return param_[i]; } \
template <> const typename eIface<fvPatchField<ElemType> >::boundaryType& eIface<fvPatchField<ElemType> >::getBoundary(const label) const { return param_; } \
\
template <> dimensionSet eIface<ElemType>::dimensions() const { return dimless; } \
template <> typename eIface<ElemType>::elemType eIface<ElemType>::get(const label) const { return param_; } \
template <> const typename eIface<ElemType>::boundaryType& eIface<ElemType>::getBoundary(const label) const { return param_; }


REGISTER_TYPE(volScalarField, scalar)
REGISTER_TYPE(volVectorField, vector)
REGISTER_TYPE(volTensorField, tensor)
REGISTER_TYPE(volSymmTensorField, symmTensor)

template <> struct eTraits<dimensionedScalar> { typedef scalar elemType; typedef dimensionedScalar boundaryType; enum Properties { HasMesh=0 }; };

template <> dimensionSet eIface<dimensionedScalar>::dimensions() const { return param_.dimensions(); }
template <> typename eIface<dimensionedScalar>::elemType eIface<dimensionedScalar>::get(const label) const { return param_.value(); }
template <> const typename eIface<dimensionedScalar>::boundaryType& eIface<dimensionedScalar>::getBoundary(const label) const { return param_; }

template <> struct eTraits<geometricOneField> { typedef scalar elemType; typedef geometricOneField boundaryType; enum Properties { HasMesh=0 }; };

template <> dimensionSet eIface<geometricOneField>::dimensions() const { return dimless; }
template <> typename eIface<geometricOneField>::elemType eIface<geometricOneField>::get(const label) const { return 1.0; }
template <> const typename eIface<geometricOneField>::boundaryType& eIface<geometricOneField>::getBoundary(const label) const { return param_; }


#define BINARY_OPERATION(opFunc, op) \
\
template <class lhs_t, class rhs_t> \
class opFunc \
{ \
    const lhs_t& lhs_; \
    const rhs_t& rhs_; \
public: \
    typedef typename eReturnTraits<typename eTraits<lhs_t>::elemType, typename eTraits<rhs_t>::elemType>::retType elemType; \
    typedef opFunc <typename eTraits<lhs_t>::boundaryType, typename eTraits<rhs_t>::boundaryType> boundaryExpr; \
    opFunc (const lhs_t& lhs, const rhs_t& rhs) : lhs_(lhs), rhs_(rhs) {} \
    dimensionSet dimensions() const { return eIface<lhs_t>(lhs_).dimensions() op eIface<rhs_t>(rhs_).dimensions(); } \
    elemType get(const label i) const { return eIface<lhs_t>(lhs_).get(i) op eIface<rhs_t>(rhs_).get(i); } \
    boundaryExpr getBoundary(const label i) const { return boundaryExpr(eIface<lhs_t>(lhs_).getBoundary(i), eIface<rhs_t>(rhs_).getBoundary(i)); } \
    const fvMesh& mesh() const { return getMesh(lhs_, rhs_); } \
}; \
\
template <class lhs_t, class rhs_t> \
opFunc <lhs_t, rhs_t> operator op(const lhs_t& lhs, const rhs_t& rhs) \
{ \
    return opFunc <lhs_t, rhs_t>(lhs, rhs); \
}

#define BINARY_OPERATION_RET_TYPE(opFunc, op, retType) \
\
template <class lhs_t, class rhs_t> \
class opFunc \
{ \
    const lhs_t& lhs_; \
    const rhs_t& rhs_; \
public: \
    typedef retType elemType; \
    typedef opFunc <typename eTraits<lhs_t>::boundaryType, typename eTraits<rhs_t>::boundaryType> boundaryExpr; \
    opFunc (const lhs_t& lhs, const rhs_t& rhs) : lhs_(lhs), rhs_(rhs) {} \
    dimensionSet dimensions() const { return eIface<lhs_t>(lhs_).dimensions() op eIface<rhs_t>(rhs_).dimensions(); } \
    elemType get(const label i) const { return eIface<lhs_t>(lhs_).get(i) op eIface<rhs_t>(rhs_).get(i); } \
    boundaryExpr getBoundary(const label i) const { return boundaryExpr(eIface<lhs_t>(lhs_).getBoundary(i), eIface<rhs_t>(rhs_).getBoundary(i)); } \
    const fvMesh& mesh() const { return getMesh(lhs_, rhs_); } \
}; \
\
template <class lhs_t, class rhs_t> \
opFunc <lhs_t, rhs_t> operator op(const lhs_t& lhs, const rhs_t& rhs) \
{ \
    return opFunc <lhs_t, rhs_t>(lhs, rhs); \
}

#define BINARY_FUNCTION(opFunc, func) \
\
template <class lhs_t, class rhs_t> \
class opFunc \
{ \
    const lhs_t& lhs_; \
    const rhs_t& rhs_; \
public: \
    typedef typename eReturnTraits<typename eTraits<lhs_t>::elemType, typename eTraits<rhs_t>::elemType>::retType elemType; \
    typedef opFunc <typename eTraits<lhs_t>::boundaryType, typename eTraits<rhs_t>::boundaryType> boundaryExpr; \
    opFunc (const lhs_t& lhs, const rhs_t& rhs) : lhs_(lhs), rhs_(rhs) {} \
    dimensionSet dimensions() const { return func (eIface<lhs_t>(lhs_).dimensions(), eIface<rhs_t>(rhs_).dimensions()); } \
    elemType get(const label i) const { return func (eIface<lhs_t>(lhs_).get(i), eIface<rhs_t>(rhs_).get(i)); } \
    boundaryExpr getBoundary(const label i) const { return boundaryExpr(eIface<lhs_t>(lhs_).getBoundary(i), eIface<rhs_t>(rhs_).getBoundary(i)); } \
    const fvMesh& mesh() const { return getMesh(lhs_, rhs_); } \
}; \
\
template <class lhs_t, class rhs_t> \
opFunc <lhs_t, rhs_t> func(const lhs_t& lhs, const rhs_t& rhs) \
{ \
    return opFunc <lhs_t, rhs_t>(lhs, rhs); \
}

#define UNARY_FUNCTION(opFunc, func, retType) \
\
template <class T> \
class opFunc \
{ \
    const T& v_; \
public: \
    typedef retType elemType; \
    typedef opFunc <typename eTraits<T>::boundaryType> boundaryExpr; \
    opFunc (const T& v) : v_(v) {} \
    dimensionSet dimensions() const { return func (eIface<T>(v_).dimensions()); } \
    elemType get(const label i) const { return func (eIface<T>(v_).get(i)); } \
    boundaryExpr getBoundary(const label i) const { return boundaryExpr(eIface<T>(v_).getBoundary(i)); } \
    const fvMesh& mesh() const { return v_.mesh(); } \
}; \
\
template <class T> \
opFunc <T> func(const T& v) \
{ \
    return opFunc <T>(v); \
}

#define UNARY_FUNCTION_NO_DIM_CHANGE(opFunc, func, retType) \
\
template <class T> \
class opFunc \
{ \
    const T& v_; \
public: \
    typedef retType elemType; \
    typedef opFunc <typename eTraits<T>::boundaryType> boundaryExpr; \
    opFunc (const T& v) : v_(v) {} \
    dimensionSet dimensions() const { return eIface<T>(v_).dimensions(); } \
    elemType get(const label i) const { return func (eIface<T>(v_).get(i)); } \
    boundaryExpr getBoundary(const label i) const { return boundaryExpr(eIface<T>(v_).getBoundary(i)); } \
    const fvMesh& mesh() const { return v_.mesh(); } \
}; \
\
template <class T> \
opFunc <T> func(const T& v) \
{ \
    return opFunc <T>(v); \
}

BINARY_OPERATION(multiply_t, *)
BINARY_OPERATION(divide_t, /)
BINARY_OPERATION(add_t, +)
BINARY_OPERATION(subtract_t, -)

BINARY_OPERATION_RET_TYPE(doubleDotProduct_t, &&, scalar)
BINARY_OPERATION_RET_TYPE(dotProduct_t, &, scalar)

BINARY_FUNCTION(minimum_t, min)
BINARY_FUNCTION(maximum_t, max)

UNARY_FUNCTION_NO_DIM_CHANGE(symm_t, symm, symmTensor)
UNARY_FUNCTION_NO_DIM_CHANGE(twoSymm_t, twoSymm, symmTensor)
UNARY_FUNCTION_NO_DIM_CHANGE(dev_t, dev, symmTensor)
UNARY_FUNCTION(magSqr_t, magSqr, scalar)

UNARY_FUNCTION_NO_DIM_CHANGE(tanh_t, tanh, scalar)
UNARY_FUNCTION(pow4_t, pow4, scalar)
UNARY_FUNCTION(sqrt_t, sqrt, scalar)
UNARY_FUNCTION(sqr_t, sqr, scalar)

template <class boundaryField, class expr_t>
void evalB(boundaryField& res, const expr_t& expr)
{
    const label sz = res.size();
    for(int i = 0; i < sz; ++i)
    {
        res[i] = expr.get(i);
    }
}


template <class expr_t>
void eval(volScalarField& res, const expr_t& expr)
{
    const label sz = res.size();
    for(int i = 0; i < sz; ++i)
    {
        res[i] = expr.get(i);
    }

    const label nbf = res.boundaryField().size();
    for(int b = 0; b < nbf; ++b)
    {
        evalB(res.boundaryField()[b], expr.getBoundary(b));
    }

    res.dimensions().reset(expr.dimensions());
}


template <class expr_t>
tmp<volScalarField> eval(const char* name, const expr_t& expr)
{
    tmp<volScalarField> res(
        new volScalarField(
            IOobject
            (
                name,
                expr.mesh().time().timeName(),
                expr.mesh()
            ),
            expr.mesh(),
            dimensionedScalar("eval", expr.dimensions(), 0.0)
        )
    );

    eval(res(), expr);

    return res;
}

template <class expr_t>
tmp<volScalarField> eval(const expr_t& expr)
{
    return eval("eval", expr);
}


} // end of namespace Expression

} // end of namespace Foam

#endif

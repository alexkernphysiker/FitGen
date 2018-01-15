// this file is distributed under
// LGPLv3 license
#ifndef ____mVYpymgQ
#define ____mVYpymgQ
#if __cplusplus<201100L
#error c++>=11 is needed for using math_h headers
#endif
#include <math.h>
#include <functional>
#include <math_h/functions.h>
#include "fit.h"
namespace Genetic
{
template<int a, int b>
struct max2 {
    enum {val = (a > b) ? a : b};
};
template<int a, int b, int c>
struct max3 {
    enum {val = (max2<a, b>::val > c) ? max2<a, b>::val : c};
};
template<int a, int b, int c, int d>
struct max4 {
    enum {val = (max3<a, b, c>::val > d) ? max3<a, b, c>::val : d};
};
template<int a, int b, int c, int d, int e>
struct max5 {
    enum {val = (max4<a, b, c, d>::val > e) ? max4<a, b, c, d>::val : e};
};
template<int a, int b, int c, int d, int e, int f>
struct max6 {
    enum {val = (max5<a, b, c, d, e>::val > f) ? max5<a, b, c, d, e>::val : f};
};
enum Recuring {first, second, third, fourth, fifth, sixth};
template<int value, Recuring recurring = first>
class Const: public virtual IParamFunc
{
public:
    Const() {}
    virtual ~Const() {}
    virtual double operator()(const ParamSet &, const ParamSet &)const override
    {
        return double(value);
    }
    enum {ParamCount = 0, ArgCount = 0};
};
template<int x_index, Recuring recurring = first>
class Arg: public virtual IParamFunc
{
public:
    Arg() {}
    virtual ~Arg() {}
    virtual double operator()(const ParamSet &X, const ParamSet &)const override
    {
        return X[x_index];
    }
    enum {ParamCount = 0, ArgCount = x_index + 1};
};
template<int p_index, Recuring recurring = first>
class Par: public virtual IParamFunc
{
public:
    Par() {}
    virtual ~Par() {}
    virtual double operator()(const ParamSet &, const ParamSet &P)const override
    {
        return P[p_index];
    }
    enum {ParamCount = p_index + 1, ArgCount = 0};
};
template<class FUNC1, int p_index, unsigned int power, Recuring recurring = first>
class PolynomFunc: public virtual FUNC1
{
public:
    PolynomFunc() {}
    virtual ~PolynomFunc() {}
    virtual double operator()(const ParamSet &X, const ParamSet &P)const override
    {
        return MathTemplates::Polynom<power, p_index, double, ParamSet>(FUNC1::operator()(X, P), P);
    }
    enum {ParamCount = max2 < FUNC1::ParamCount, p_index + power + 1 >::val, ArgCount = FUNC1::ArgCount};
};
template<double(func)(const double &), class FUNC, Recuring recurring = first>
class Func: public virtual FUNC
{
public:
    Func(): FUNC() {}
    virtual ~Func() {}
    virtual double operator()(const ParamSet &X, const ParamSet &P)const override
    {
        return func(FUNC::operator()(X, P));
    }
    enum {ParamCount = FUNC::ParamCount, ArgCount = FUNC::ArgCount};
};
template<double(func)(const double &, const double &), class FUNC1, class FUNC2, Recuring recurring = first>
class Func2: public virtual FUNC1, public virtual FUNC2
{
public:
    Func2(): FUNC1(), FUNC2() {}
    virtual ~Func2() {}
    virtual double operator()(const ParamSet &X, const ParamSet &P)const override
    {
        return func(FUNC1::operator()(X, P), FUNC2::operator()(X, P));
    }
    enum {
        ParamCount = max2<FUNC1::ParamCount, FUNC2::ParamCount>::val,
        ArgCount = max2<FUNC1::ArgCount, FUNC2::ArgCount>::val
    };
};
template<double(func)(const double &, const double &, const double &), class FUNC1, class FUNC2, class FUNC3, Recuring recurring = first>
class Func3: public virtual FUNC1, public virtual FUNC2, public virtual FUNC3
{
public:
    Func3(): FUNC1(), FUNC2(), FUNC3() {}
    virtual ~Func3() {}
    virtual double operator()(const ParamSet &X, const ParamSet &P)const override
    {
        return func(FUNC1::operator()(X, P), FUNC2::operator()(X, P), FUNC3::operator()(X, P));
    }
    enum {
        ParamCount = max3<FUNC1::ParamCount, FUNC2::ParamCount, FUNC3::ParamCount>::val,
        ArgCount = max3<FUNC1::ArgCount, FUNC2::ArgCount, FUNC3::ArgCount>::val
    };
};
template<double(func)(const double &, const double &, const double &, const double &), class FUNC1, class FUNC2, class FUNC3, class FUNC4, Recuring recurring = first>
class Func4: public virtual FUNC1, public virtual FUNC2, public virtual FUNC3, public virtual FUNC4
{
public:
    Func4(): FUNC1(), FUNC2(), FUNC3(), FUNC4() {}
    virtual ~Func4() {}
    virtual double operator()(const ParamSet &X, const ParamSet &P)const override
    {
        return func(FUNC1::operator()(X, P), FUNC2::operator()(X, P), FUNC3::operator()(X, P),
                    FUNC4::operator()(X, P));
    }
    enum {
        ParamCount = max4<FUNC1::ParamCount, FUNC2::ParamCount, FUNC3::ParamCount, FUNC4::ParamCount>::val,
        ArgCount = max4<FUNC1::ArgCount, FUNC2::ArgCount, FUNC3::ArgCount, FUNC4::ArgCount>::val
    };
};
template<double(func)(const double &, const double &, const double &, const double &, const double &), class FUNC1, class FUNC2, class FUNC3, class FUNC4, class FUNC5, Recuring recurring = first>
class Func5: public virtual FUNC1, public virtual FUNC2, public virtual FUNC3, public virtual FUNC4, public virtual FUNC5
{
public:
    Func5(): FUNC1(), FUNC2(), FUNC3(), FUNC4(), FUNC5() {}
    virtual ~Func5() {}
    virtual double operator()(const ParamSet &X, const ParamSet &P)const override
    {
        return func(FUNC1::operator()(X, P), FUNC2::operator()(X, P), FUNC3::operator()(X, P),
                    FUNC4::operator()(X, P), FUNC5::operator()(X, P));
    }
    enum {
        ParamCount = max5<FUNC1::ParamCount, FUNC2::ParamCount, FUNC3::ParamCount, FUNC4::ParamCount, FUNC5::ParamCount>::val,
        ArgCount = max5<FUNC1::ArgCount, FUNC2::ArgCount, FUNC3::ArgCount, FUNC4::ArgCount, FUNC5::ArgCount>::val
    };
};
template<double(func)(const double &, const double &, const double &, const double &, const double &, const double &), class FUNC1, class FUNC2, class FUNC3, class FUNC4, class FUNC5, class FUNC6, Recuring recurring = first>
class Func6: public virtual FUNC1, public virtual FUNC2, public virtual FUNC3, public virtual FUNC4, public virtual FUNC5, public virtual FUNC6
{
public:
    Func6(): FUNC1(), FUNC2(), FUNC3(), FUNC4(), FUNC5(), FUNC6() {}
    virtual ~Func6() {}
    virtual double operator()(const ParamSet &X, const ParamSet &P)const override
    {
        return func(FUNC1::operator()(X, P), FUNC2::operator()(X, P), FUNC3::operator()(X, P),
                    FUNC4::operator()(X, P), FUNC5::operator()(X, P), FUNC6::operator()(X, P));
    }
    enum {
        ParamCount = max6<FUNC1::ParamCount, FUNC2::ParamCount, FUNC3::ParamCount, FUNC4::ParamCount, FUNC5::ParamCount, FUNC6::ParamCount>::val,
        ArgCount = max6<FUNC1::ArgCount, FUNC2::ArgCount, FUNC3::ArgCount, FUNC4::ArgCount, FUNC5::ArgCount, FUNC6::ArgCount>::val
    };
};
template<class FUNC, int x_index, int p_index, Recuring recurring = first>
class ArgShift: public virtual FUNC
{
public:
    ArgShift(): FUNC() {} virtual ~ArgShift() {}
    virtual double operator()(const ParamSet &X, const ParamSet &P)const override
    {
        ParamSet x = X;
        x(x_index) += P[p_index];
        return FUNC::operator()((x), (P));
    }
    enum {ParamCount = FUNC::ParamCount, ArgCount = FUNC::ArgCount};
};
template<class FUNC, int x_index, int p_index, Recuring recurring = first>
class ArgScale: public virtual FUNC
{
public:
    ArgScale(): FUNC() {}
    virtual ~ArgScale() {}
    virtual double operator()(const ParamSet &X, const ParamSet &P)const override
    {
        ParamSet x = X;
        x(x_index) *= P[p_index];
        return FUNC::operator()((x), (P));
    }
    enum {ParamCount = FUNC::ParamCount, ArgCount = FUNC::ArgCount};
};
template<class FUNC, Recuring recurring = first>
class Minus: public virtual FUNC
{
public:
    Minus(): FUNC() {}
    virtual ~Minus() {}
    virtual double operator()(const ParamSet &X, const ParamSet &P)const override
    {
        return -FUNC::operator()((X), (P));
    }
    enum {
        ParamCount = FUNC::ParamCount,
        ArgCount = FUNC::ArgCount
    };
};
template<class FUNC1, class FUNC2, Recuring recurring = first>
class Add: public virtual FUNC1, public virtual FUNC2
{
public:
    Add(): FUNC1(), FUNC2() {}
    virtual ~Add() {}
    virtual double operator()(const ParamSet &X, const ParamSet &P)const override
    {
        return FUNC1::operator()((X), (P)) + FUNC2::operator()((X), (P));
    }
    enum {
        ParamCount = max2<FUNC1::ParamCount, FUNC2::ParamCount>::val,
        ArgCount = max2<FUNC1::ArgCount, FUNC2::ArgCount>::val
    };
};
template<class FUNC1, class FUNC2, Recuring recurring = first>class Add2: public Add<FUNC1, FUNC2, recurring> {};
template<class FUNC1, class FUNC2, class FUNC3, Recuring recurring = first>
class Add3: public virtual FUNC1, public virtual FUNC2, public virtual FUNC3
{
public:
    Add3(): FUNC1(), FUNC2(), FUNC3() {}
    virtual ~Add3() {}
    virtual double operator()(const ParamSet &X, const ParamSet &P)const override
    {
        return FUNC1::operator()((X), (P)) + FUNC2::operator()((X), (P)) + FUNC3::operator()((X), (P));
    }
    enum {
        ParamCount = max3<FUNC1::ParamCount, FUNC2::ParamCount, FUNC3::ParamCount>::val,
        ArgCount = max3<FUNC1::ArgCount, FUNC2::ArgCount, FUNC3::ArgCount>::val
    };
};
template<class FUNC1, class FUNC2, class FUNC3, class FUNC4, Recuring recurring = first>
class Add4: public virtual FUNC1, public virtual FUNC2, public virtual FUNC3, public virtual FUNC4
{
public:
    Add4(): FUNC1(), FUNC2(), FUNC3(), FUNC4() {}
    virtual ~Add4() {}
    virtual double operator()(const ParamSet &X, const ParamSet &P)const override
    {
        return FUNC1::operator()((X), (P)) + FUNC2::operator()((X), (P)) + FUNC3::operator()((X), (P))
               + FUNC4::operator()((X), (P));
    }
    enum {
        ParamCount = max4<FUNC1::ParamCount, FUNC2::ParamCount, FUNC3::ParamCount
        , FUNC4::ParamCount>::val,
        ArgCount = max4<FUNC1::ArgCount, FUNC2::ArgCount, FUNC3::ArgCount
        , FUNC4::ArgCount>::val
    };
};
template<class FUNC1, class FUNC2, class FUNC3, class FUNC4, class FUNC5, Recuring recurring = first>
class Add5: public virtual FUNC1, public virtual FUNC2, public virtual FUNC3
    , public virtual FUNC4, public virtual FUNC5
{
public:
    Add5(): FUNC1(), FUNC2(), FUNC3(), FUNC4(), FUNC5() {}
    virtual ~Add5() {}
    virtual double operator()(const ParamSet &X, const ParamSet &P)const override
    {
        return FUNC1::operator()((X), (P)) + FUNC2::operator()((X), (P)) + FUNC3::operator()((X), (P))
               + FUNC4::operator()((X), (P)) + FUNC5::operator()((X), (P));
    }
    enum {
        ParamCount = max5<FUNC1::ParamCount, FUNC2::ParamCount, FUNC3::ParamCount
        , FUNC4::ParamCount, FUNC5::ParamCount>::val,
        ArgCount = max5<FUNC1::ArgCount, FUNC2::ArgCount, FUNC3::ArgCount
        , FUNC4::ArgCount, FUNC5::ArgCount>::val
    };
};
template<class FUNC1, class FUNC2, class FUNC3, class FUNC4, class FUNC5, class FUNC6, Recuring recurring = first>
class Add6: public virtual FUNC1, public virtual FUNC2, public virtual FUNC3
    , public virtual FUNC4, public virtual FUNC5, public virtual FUNC6
{
public:
    Add6(): FUNC1(), FUNC2(), FUNC3(), FUNC4(), FUNC5(), FUNC6() {}
    virtual ~Add6() {}
    virtual double operator()(const ParamSet &X, const ParamSet &P)const override
    {
        return FUNC1::operator()((X), (P)) + FUNC2::operator()((X), (P)) + FUNC3::operator()((X), (P))
               + FUNC4::operator()((X), (P)) + FUNC5::operator()((X), (P)) + FUNC6::operator()((X), (P));
    }
    enum {
        ParamCount = max6<FUNC1::ParamCount, FUNC2::ParamCount, FUNC3::ParamCount
        , FUNC4::ParamCount, FUNC5::ParamCount, FUNC6::ParamCount>::val,
        ArgCount = max6<FUNC1::ArgCount, FUNC2::ArgCount, FUNC3::ArgCount
        , FUNC4::ArgCount, FUNC5::ArgCount, FUNC6::ArgCount>::val
    };
};

template<class FUNC1, class FUNC2, Recuring recurring = first>
class Sub: public virtual FUNC1, public virtual FUNC2
{
public:
    Sub(): FUNC1(), FUNC2() {}
    virtual ~Sub() {}
    virtual double operator()(const ParamSet &X, const ParamSet &P)const override
    {
        return FUNC1::operator()((X), (P)) - FUNC2::operator()((X), (P));
    }
    enum {
        ParamCount = max2<FUNC1::ParamCount, FUNC2::ParamCount>::val,
        ArgCount = max2<FUNC1::ArgCount, FUNC2::ArgCount>::val
    };
};
template<class FUNC1, class FUNC2, Recuring recurring = first>
class Mul: public virtual FUNC1, public virtual FUNC2
{
public:
    Mul(): FUNC1(), FUNC2() {}
    virtual ~Mul() {}
    virtual double operator()(const ParamSet &X, const ParamSet &P)const override
    {
        return FUNC1::operator()((X), (P)) * FUNC2::operator()((X), (P));
    }
    enum {
        ParamCount = max2<FUNC1::ParamCount, FUNC2::ParamCount>::val,
        ArgCount = max2<FUNC1::ArgCount, FUNC2::ArgCount>::val
    };
};
template<class FUNC1, class FUNC2, Recuring recurring = first>class Mul2: public Mul<FUNC1, FUNC2, recurring> {};
template<class FUNC1, class FUNC2, class FUNC3, Recuring recurring = first>
class Mul3: public virtual FUNC1, public virtual FUNC2, public virtual FUNC3
{
public:
    Mul3(): FUNC1(), FUNC2(), FUNC3() {}
    virtual ~Mul3() {}
    virtual double operator()(const ParamSet &X, const ParamSet &P)const override
    {
        return FUNC1::operator()((X), (P)) * FUNC2::operator()((X), (P)) * FUNC3::operator()((X), (P));
    }
    enum {
        ParamCount = max3<FUNC1::ParamCount, FUNC2::ParamCount, FUNC3::ParamCount>::val,
        ArgCount = max3<FUNC1::ArgCount, FUNC2::ArgCount, FUNC3::ArgCount>::val
    };
};
template<class FUNC1, class FUNC2, class FUNC3, class FUNC4, Recuring recurring = first>
class Mul4: public virtual FUNC1, public virtual FUNC2, public virtual FUNC3, public virtual FUNC4
{
public:
    Mul4(): FUNC1(), FUNC2(), FUNC3(), FUNC4() {}
    virtual ~Mul4() {}
    virtual double operator()(const ParamSet &X, const ParamSet &P)const override
    {
        return FUNC1::operator()((X), (P)) * FUNC2::operator()((X), (P)) * FUNC3::operator()((X), (P))
               * FUNC4::operator()((X), (P));
    }
    enum {
        ParamCount = max4<FUNC1::ParamCount, FUNC2::ParamCount, FUNC3::ParamCount
        , FUNC4::ParamCount>::val,
        ArgCount = max4<FUNC1::ArgCount, FUNC2::ArgCount, FUNC3::ArgCount
        , FUNC4::ArgCount>::val
    };
};
template<class FUNC1, class FUNC2, class FUNC3, class FUNC4, class FUNC5, Recuring recurring = first>
class Mul5: public virtual FUNC1, public virtual FUNC2, public virtual FUNC3
    , public virtual FUNC4, public virtual FUNC5
{
public:
    Mul5(): FUNC1(), FUNC2(), FUNC3(), FUNC4(), FUNC5() {}
    virtual ~Mul5() {}
    virtual double operator()(const ParamSet &X, const ParamSet &P)const override
    {
        return FUNC1::operator()((X), (P)) * FUNC2::operator()((X), (P)) * FUNC3::operator()((X), (P))
               * FUNC4::operator()((X), (P)) * FUNC5::operator()((X), (P));
    }
    enum {
        ParamCount = max5<FUNC1::ParamCount, FUNC2::ParamCount, FUNC3::ParamCount
        , FUNC4::ParamCount, FUNC5::ParamCount>::val,
        ArgCount = max5<FUNC1::ArgCount, FUNC2::ArgCount, FUNC3::ArgCount
        , FUNC4::ArgCount, FUNC5::ArgCount>::val
    };
};
template<class FUNC1, class FUNC2, class FUNC3, class FUNC4, class FUNC5, class FUNC6, Recuring recurring = first>
class Mul6: public virtual FUNC1, public virtual FUNC2, public virtual FUNC3
    , public virtual FUNC4, public virtual FUNC5, public virtual FUNC6
{
public:
    Mul6(): FUNC1(), FUNC2(), FUNC3(), FUNC4(), FUNC5(), FUNC6() {}
    virtual ~Mul6() {}
    virtual double operator()(const ParamSet &X, const ParamSet &P)const override
    {
        return FUNC1::operator()((X), (P)) * FUNC2::operator()((X), (P)) * FUNC3::operator()((X), (P))
               * FUNC4::operator()((X), (P)) * FUNC5::operator()((X), (P)) * FUNC6::operator()((X), (P));
    }
    enum {
        ParamCount = max6<FUNC1::ParamCount, FUNC2::ParamCount, FUNC3::ParamCount
        , FUNC4::ParamCount, FUNC5::ParamCount, FUNC6::ParamCount>::val,
        ArgCount = max6<FUNC1::ArgCount, FUNC2::ArgCount, FUNC3::ArgCount
        , FUNC4::ArgCount, FUNC5::ArgCount, FUNC6::ArgCount>::val
    };
};




template<class FUNC1, class FUNC2, Recuring recurring = first>
class Div: public virtual FUNC1, public virtual FUNC2
{
public:
    Div(): FUNC1(), FUNC2() {}
    virtual ~Div() {}
    virtual double operator()(const ParamSet &X, const ParamSet &P)const override
    {
        return FUNC1::operator()((X), (P)) / FUNC2::operator()((X), (P));
    }
    enum {
        ParamCount = max2<FUNC1::ParamCount, FUNC2::ParamCount>::val,
        ArgCount = max2<FUNC1::ArgCount, FUNC2::ArgCount>::val
    };
};
template<class FUNC1, class FUNC2, Recuring recurring = first>
class Power: public virtual FUNC1, public virtual FUNC2
{
public:
    Power(): FUNC1(), FUNC2() {}
    virtual ~Power() {}
    virtual double operator()(const ParamSet &X, const ParamSet &P)const override
    {
        return pow(FUNC1::operator()((X), (P)), FUNC2::operator()((X), (P)));
    }
    enum {
        ParamCount = max2<FUNC1::ParamCount, FUNC2::ParamCount>::val,
        ArgCount = max2<FUNC1::ArgCount, FUNC2::ArgCount>::val
    };
};

template<class FUNC1, Recuring recurring = first>
class Sqrt: public virtual FUNC1
{
public:
    Sqrt(): FUNC1() {}
    virtual ~Sqrt() {}
    virtual double operator()(const ParamSet &X, const ParamSet &P)const override
    {
        return sqrt(FUNC1::operator()((X), (P)));
    }
    enum {
        ParamCount = FUNC1::ParamCount, ArgCount = FUNC1::ArgCount
    };
};
template<class FUNC1, Recuring recurring = first>
class Sqr: public virtual FUNC1
{
public:
    Sqr(): FUNC1() {}
    virtual ~Sqr() {}
    virtual double operator()(const ParamSet &X, const ParamSet &P)const override
    {
        return pow(FUNC1::operator()((X), (P)), 2);
    }
    enum {
        ParamCount = FUNC1::ParamCount, ArgCount = FUNC1::ArgCount
    };
};
template<class FUNC1, Recuring recurring = first>
class Pow3: public virtual FUNC1
{
public:
    Pow3(): FUNC1() {}
    virtual ~Pow3() {}
    virtual double operator()(const ParamSet &X, const ParamSet &P)const override
    {
        return pow(FUNC1::operator()((X), (P)), 3);
    }
    enum {
        ParamCount = FUNC1::ParamCount, ArgCount = FUNC1::ArgCount
    };
};
template<class FUNC1, Recuring recurring = first>
class Pow4: public virtual FUNC1
{
public:
    Pow4(): FUNC1() {}
    virtual ~Pow4() {}
    virtual double operator()(const ParamSet &X, const ParamSet &P)const override
    {
        return pow(FUNC1::operator()((X), (P)), 4);
    }
    enum {
        ParamCount = FUNC1::ParamCount, ArgCount = FUNC1::ArgCount
    };
};
template<class FUNC1, Recuring recurring = first>
class Pow5: public virtual FUNC1
{
public:
    Pow5(): FUNC1() {}
    virtual ~Pow5() {}
    virtual double operator()(const ParamSet &X, const ParamSet &P)const override
    {
        return pow(FUNC1::operator()((X), (P)), 5);
    }
    enum {
        ParamCount = FUNC1::ParamCount, ArgCount = FUNC1::ArgCount
    };
};
template<class FUNC1, Recuring recurring = first>
class Pow6: public virtual FUNC1
{
public:
    Pow6(): FUNC1() {}
    virtual ~Pow6() {}
    virtual double operator()(const ParamSet &X, const ParamSet &P)const override
    {
        return pow(FUNC1::operator()((X), (P)), 6);
    }
    enum {
        ParamCount = FUNC1::ParamCount, ArgCount = FUNC1::ArgCount
    };
};
}
#endif

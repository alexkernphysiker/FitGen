#ifndef EQUATION_H
#define EQUATION_H
#include <math.h>
#include "fit_gen.h"
namespace Fit {
//include templates for creating complicated functions double (func)(ParamSet&)
#define use_num_type double
#define use_indexer_type ParamSet&
#include <math_h/wrap1.h>
#undef use_num_type
#undef use_indexer_type
// "empty" IParamFunc for giving as a parameter when it's not needed
class NoParamFunc:public IParamFunc{
public:
	NoParamFunc(){}virtual ~NoParamFunc(){}
	virtual double operator()(ParamSet&,ParamSet&) override{return 0;}
	virtual bool CorrectParams(ParamSet&) override{return true;}
};
// IOptimalityFunction for solving equation func(P)=0
template<double (func)(ParamSet&)>
class Equation:public IOptimalityFunction{
public:Equation(){}virtual ~Equation(){}
	virtual double operator()(ParamSet &P, IParamFunc&)override{
		double res=func(P);if(res>=0)return res; else return -res;
	}
};
// IOptimalityFunction for solving equation func1(P)=func2(P)
template<double (func1)(ParamSet&),double (func2)(ParamSet&)>
class Equation2:public IOptimalityFunction{
public:Equation2(){}virtual ~Equation2(){}
	virtual double operator()(ParamSet &P, IParamFunc&)override{
		double res=func1(P)-func2(P);if(res>=0)return res; else return -res;
	}
};
// IOptimalityFunction for searching minimum of func(P)
template<double (func)(ParamSet&)>
class SearchMin:public IOptimalityFunction{
public:SearchMin(){}virtual ~SearchMin(){}
	virtual double operator()(ParamSet &P, IParamFunc&)override{return func(P);}
};
}
#endif // EQUATION_H
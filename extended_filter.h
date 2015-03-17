#ifndef TXCBIZOCGJTYMMMX
#define TXCBIZOCGJTYMMMX
#include "fit_gen.h"
namespace Fit{
	using namespace std;
	enum condition{EQ,NE,GR,NG,LE,NL};
	inline bool TakeCondition(double a,condition c, double b){
		switch(c){
			case EQ:return a==b;
			case NE:return a!=b;
			case GR:return a>b;
			case NG:return a<=b;
			case LE:return a<b;
			case NL:return a>=b;
		}
	}
	inline double ZERO(ParamSet&){return 0;}
	template<double (func1)(ParamSet&),condition c,double (func2)(ParamSet&)>
	class FuncCondition:IParamCheck{
	public:
		FuncCondition(){}
		virtual ~FuncCondition(){}
		virtual bool CorrectParams(ParamSet params)override{
			return TakeCondition(func1(params),c,func2(params));
		}
	};
	template<class Func1,class Func2>
	class FilterCondition:public IParamCheck{
	private:
		Func1 func1;
		Func2 func2;
		condition c;
	public:
		FilterCondition(Func1 f1, condition cond,Func2 f2){
			func1=f1;
			func2=f2;
			c=cond;
		}
		virtual ~FilterCondition(){}
		virtual bool CorrectParams(ParamSet params)override{
			return TakeCondition(func1(params),c,func2(params));
		}
	};
	template<class FUNC>
	class ParamWrap{
	private:
		shared_ptr<FUNC> func;
		ParamSet X;
	public:
		ParamWrap(){
			func=make_shared<FUNC>();
			X=parZeros(FUNC::ParamCount);
		}
		ParamWrap(const ParamWrap &C){
			func=C.func;
			X=C.X;
		}
		~ParamWrap(){}
		double operator()(ParamSet P){
			return func->operator()(X,P);
		}
	};
	template<class FUNC1,condition c, class FUNC2>
	shared_ptr<IParamCheck> Condition(){
		return make_shared<FilterCondition<ParamWrap<FUNC1>,ParamWrap<FUNC2>>>(ParamWrap<FUNC1>(),c,ParamWrap<FUNC2>());
	}
}
#endif
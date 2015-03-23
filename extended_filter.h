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
	class ParamWrap:FUNC{
	protected:
		ParamSet X;
	public:
		enum{ParamCount=FUNC::ParamCount,ArgCount=FUNC::ArgCount};
		ParamWrap():FUNC(){
			X=parZeros(ArgCount);
		}
		ParamWrap(const ParamWrap &C){
			X=C.X;
		}
		virtual ~ParamWrap(){}
		virtual double operator()(ParamSet P){
			return FUNC::operator()(X,P);
		}
	};
	template<class FUNC,int arg0>
	class ParamWrap1:public ParamWrap<FUNC>{
	public:
		ParamWrap1():ParamWrap<FUNC>(){}
		ParamWrap1(const ParamWrap1 &C):ParamWrap<FUNC>(C){}
		virtual ~ParamWrap1(){}
		virtual double operator()(ParamSet P)override{
			ParamWrap<FUNC>::X.Set(0,P[arg0]);
			return ParamWrap<FUNC>::operator()(P);
		}
	};
	template<class FUNC,int arg0,int arg1>
	class ParamWrap2:public ParamWrap1<FUNC,arg0>{
	public:
		ParamWrap2():ParamWrap1<FUNC,arg0>(){}
		ParamWrap2(const ParamWrap2 &C):ParamWrap1<FUNC,arg0>(C){}
		virtual ~ParamWrap2(){}
		virtual double operator()(ParamSet P)override{
			ParamWrap<FUNC>::X.Set(0,P[arg1]);
			return ParamWrap1<FUNC,arg0>::operator()(P);
		}
	};
	template<class FUNC1,condition c, class FUNC2>
	shared_ptr<IParamCheck> Condition(){
		return make_shared<FilterCondition<FUNC1,FUNC2>>(FUNC1(),c,FUNC2());
	}
}
#endif
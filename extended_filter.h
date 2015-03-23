#ifndef TXCBIZOCGJTYMMMX
#define TXCBIZOCGJTYMMMX
#include "fit_gen.h"
#include "paramfunc.h"
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
		ParamWrap(const ParamWrap &C):X(C.X){}
		virtual ~ParamWrap(){}
		virtual double operator()(ParamSet P){
			return FUNC::operator()(X,P);
		}
	};
	template<class FUNC,double (arg0)(ParamSet&)>
	class ParamWrap1:public ParamWrap<FUNC>{
	public:
		ParamWrap1():ParamWrap<FUNC>(){}
		ParamWrap1(const ParamWrap1 &C):ParamWrap<FUNC>(C){}
		virtual ~ParamWrap1(){}
		virtual double operator()(ParamSet P)override{
			ParamWrap<FUNC>::X.Set(0,arg0(P));
			return ParamWrap<FUNC>::operator()(P);
		}
	};
	template<class FUNC,double (arg0)(ParamSet&),double (arg1)(ParamSet&)>
	class ParamWrap2:public ParamWrap1<FUNC,arg0>{
	public:
		ParamWrap2():ParamWrap1<FUNC,arg0>(){}
		ParamWrap2(const ParamWrap2 &C):ParamWrap1<FUNC,arg0>(C){}
		virtual ~ParamWrap2(){}
		virtual double operator()(ParamSet P)override{
			ParamWrap<FUNC>::X.Set(1,arg1(P));
			return ParamWrap1<FUNC,arg0>::operator()(P);
		}
	};
	template<class FUNC,double (arg0)(ParamSet&),double (arg1)(ParamSet&),
		double (arg2)(ParamSet&)>
	class ParamWrap3:public ParamWrap2<FUNC,arg0,arg1>{
	public:
		ParamWrap3():ParamWrap2<FUNC,arg0,arg1>(){}
		ParamWrap3(const ParamWrap3 &C):ParamWrap2<FUNC,arg0,arg1>(C){}
		virtual ~ParamWrap3(){}
		virtual double operator()(ParamSet P)override{
			ParamWrap<FUNC>::X.Set(2,arg2(P));
			return ParamWrap2<FUNC,arg0,arg1>::operator()(P);
		}
	};
	template<class FUNC,double (arg0)(ParamSet&),double (arg1)(ParamSet&),
		double (arg2)(ParamSet&),double (arg3)(ParamSet&)>
	class ParamWrap4:public ParamWrap3<FUNC,arg0,arg1,arg2>{
	public:
		ParamWrap4():ParamWrap3<FUNC,arg0,arg1,arg2>(){}
		ParamWrap4(const ParamWrap4 &C):ParamWrap3<FUNC,arg0,arg1,arg2>(C){}
		virtual ~ParamWrap4(){}
		virtual double operator()(ParamSet P)override{
			ParamWrap<FUNC>::X.Set(3,arg3(P));
			return ParamWrap3<FUNC,arg0,arg1,arg2>::operator()(P);
		}
	};
	template<class FUNC,double (arg0)(ParamSet&),double (arg1)(ParamSet&),
		double (arg2)(ParamSet&),double (arg3)(ParamSet&),
		double (arg4)(ParamSet&)>
	class ParamWrap5:public ParamWrap4<FUNC,arg0,arg1,arg2,arg3>{
	public:
		ParamWrap5():ParamWrap4<FUNC,arg0,arg1,arg2,arg3>(){}
		ParamWrap5(const ParamWrap5 &C):ParamWrap4<FUNC,arg0,arg1,arg2,arg3>(C){}
		virtual ~ParamWrap5(){}
		virtual double operator()(ParamSet P)override{
			ParamWrap<FUNC>::X.Set(4,arg4(P));
			return ParamWrap4<FUNC,arg0,arg1,arg2,arg3>::operator()(P);
		}
	};
	template<class FUNC,double (arg0)(ParamSet&),double (arg1)(ParamSet&),
		double (arg2)(ParamSet&),double (arg3)(ParamSet&),
		double (arg4)(ParamSet&),double (arg5)(ParamSet&)>
	class ParamWrap6:public ParamWrap5<FUNC,arg0,arg1,arg2,arg3,arg4>{
	public:
		ParamWrap6():ParamWrap5<FUNC,arg0,arg1,arg2,arg3,arg4>(){}
		ParamWrap6(const ParamWrap6 &C):ParamWrap5<FUNC,arg0,arg1,arg2,arg3,arg4>(C){}
		virtual ~ParamWrap6(){}
		virtual double operator()(ParamSet P)override{
			ParamWrap<FUNC>::X.Set(5,arg5(P));
			return ParamWrap5<FUNC,arg0,arg1,arg2,arg3,arg4>::operator()(P);
		}
	};
	template<class FUNC1,condition c, class FUNC2>
	shared_ptr<IParamCheck> Condition(){
		return make_shared<FilterCondition<FUNC1,FUNC2>>(FUNC1(),c,FUNC2());
	}
}
#endif
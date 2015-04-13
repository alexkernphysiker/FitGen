#ifndef CNTMUXAHBUIYZXIG
#define CNTMUXAHBUIYZXIG
#include "abstract.h"
#include "paramfunc.h"
namespace Genetic{
	using namespace std;
	template<class CONDITION=function<bool(ParamSet&)>>
	class Filter:public IParamCheck{
	private:
		CONDITION condition;
	public:
		Filter(CONDITION c){
			condition=c;
		}
		virtual ~Filter(){}
		virtual bool operator()(ParamSet&P)override{
			return condition(P);
		}
	};
	
	class Above:public IParamCheck{
	public:
		Above();
		Above(ParamSet v);
		virtual ~Above(){}
		virtual bool operator()(ParamSet&P)override;
		Above &operator<<(double value);
	private:
		ParamSet m_data;
	};
	class Below:public IParamCheck{
	public:
		Below();
		Below(ParamSet v);
		virtual ~Below(){}
		virtual bool operator()(ParamSet&P)override;
		Below &operator<<(double value);
	private:
		ParamSet m_data;
	};
	template<class Filter>
	inline shared_ptr<Filter> operator<<(shared_ptr<Filter> filter, double value){
		filter->operator<<(value);
		return filter;
	}
	
	class AbstractFilterMulti:public IParamCheck{
	public:
		virtual ~AbstractFilterMulti();
		int Count();
		IParamCheck &Get(int i);
		AbstractFilterMulti &Add(std::shared_ptr<IParamCheck> val);
	protected:
		vector<shared_ptr<IParamCheck>> m_data;
	};
	class And:public AbstractFilterMulti{
	public:
		And(){}
		virtual ~And(){}
		virtual bool operator()(ParamSet&P)override;
	};
	class Or:public AbstractFilterMulti{
	public:
		Or(){}
		virtual ~Or(){}
		virtual bool operator()(ParamSet&P)override;
	};
	template<class Filter>
	inline shared_ptr<Filter> operator<<(shared_ptr<Filter> filter, std::shared_ptr<IParamCheck> value){
		filter->Add(value);
		return filter;
	}
	
	enum condition_enum{EQ,NE,GR,NG,LE,NL};
	inline bool TakeCondition(double a,condition_enum c, double b){
		switch(c){
			case EQ:return a==b;
			case NE:return a!=b;
			case GR:return a>b;
			case NG:return a<=b;
			case LE:return a<b;
			case NL:return a>=b;
		}
	}
	template<double (func1)(ParamSet&),condition_enum c,double (func2)(ParamSet&)>
	class filterCondition:public IParamCheck{
	public:
		filterCondition(){}
		virtual ~filterCondition(){}
		virtual bool operator()(ParamSet params)override{
			return TakeCondition(func1(params),c,func2(params));
		}
	};
	template<double (func1)(ParamSet&),condition_enum c,double (func2)(ParamSet&)>
	inline shared_ptr<IParamCheck> condition(){
		return make_shared<filterCondition<func1,c,func2>>();
	}
	template<class Func1,class Func2>
	class FilterCondition:public IParamCheck{
	private:
		Func1 func1;
		Func2 func2;
		condition_enum c;
	public:
		FilterCondition(Func1 f1, condition_enum cond,Func2 f2){
			func1=f1;
			func2=f2;
			c=cond;
		}
		virtual ~FilterCondition(){}
		virtual bool operator()(ParamSet&P)override{
			return TakeCondition(func1(P),c,func2(P));
		}
	};
	template<class FUNC>
	class Wrap:FUNC{
	protected:
		ParamSet X;
	public:
		enum{ParamCount=FUNC::ParamCount,ArgCount=FUNC::ArgCount};
		Wrap():FUNC(){
			X=parZeros(ArgCount);
		}
		Wrap(const Wrap &C):X(C.X){}
		virtual ~Wrap(){}
		virtual double operator()(ParamSet&P){
			return FUNC::operator()(X,P);
		}
	};
	template<class FUNC,double (arg0)(ParamSet&)>
	class ParamWrap1:public Wrap<FUNC>{
	public:
		ParamWrap1():Wrap<FUNC>(){}
		ParamWrap1(const ParamWrap1 &C):Wrap<FUNC>(C){}
		virtual ~ParamWrap1(){}
		virtual double operator()(ParamSet &P)override{
			Wrap<FUNC>::X.Set(0,arg0(P));
			return Wrap<FUNC>::operator()(P);
		}
	};
	template<class FUNC,double (arg0)(ParamSet&),double (arg1)(ParamSet&)>
	class ParamWrap2:public ParamWrap1<FUNC,arg0>{
	public:
		ParamWrap2():ParamWrap1<FUNC,arg0>(){}
		ParamWrap2(const ParamWrap2 &C):ParamWrap1<FUNC,arg0>(C){}
		virtual ~ParamWrap2(){}
		virtual double operator()(ParamSet &P)override{
			Wrap<FUNC>::X.Set(1,arg1(P));
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
		virtual double operator()(ParamSet &P)override{
			Wrap<FUNC>::X.Set(2,arg2(P));
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
		virtual double operator()(ParamSet &P)override{
			Wrap<FUNC>::X.Set(3,arg3(P));
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
		virtual double operator()(ParamSet &P)override{
			Wrap<FUNC>::X.Set(4,arg4(P));
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
		virtual double operator()(ParamSet &P)override{
			Wrap<FUNC>::X.Set(5,arg5(P));
			return ParamWrap5<FUNC,arg0,arg1,arg2,arg3,arg4>::operator()(P);
		}
	};
	template<class FUNC1,condition_enum c, class FUNC2>
	inline shared_ptr<IParamCheck> Condition(){
		return make_shared<FilterCondition<FUNC1,FUNC2>>(FUNC1(),c,FUNC2());
	}
}
#endif

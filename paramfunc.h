#ifndef ____mVYpymgQ
#define ____mVYpymgQ
#include <math.h>
#include "fit_gen.h"
#include "math_h/functions.h"
namespace Fit{
	using namespace std;
	#define use_num_type double
	#define use_indexer_type ParamSet&
	#include "math_h/wrap2.h"
	#undef use_num_type
	#undef use_indexer_type
	template<class FUNC>
	class PARAMFUNC:public virtual IParamFunc{
	private:
		FUNC Func;
	public:
		PARAMFUNC(FUNC func):Func(func){}
		virtual ~PARAMFUNC(){}
		virtual double operator()(ParamSet& X,ParamSet& P)override{return Func(X,P);}
		virtual bool CorrectParams(ParamSet&)override{return true;}
	};
	template<class FUNC,class FILTER>
	class PARAM_FUNC:public virtual PARAMFUNC<FUNC>{
	private:
		FILTER Filter;
	public:
		PARAM_FUNC(FUNC func,FILTER filter):PARAMFUNC<FUNC>(func),Filter(filter){}
		virtual ~PARAM_FUNC(){}
		virtual bool CorrectParams(ParamSet &P)override{return Filter(P);}
	};
	template<double (func)(ParamSet&,ParamSet&)>
	class ParamFunc:public virtual IParamFunc{
	public:
		ParamFunc(){}virtual ~ParamFunc(){}
		virtual double operator()(ParamSet &X, ParamSet &P) override{return func(X,P);}
		virtual bool CorrectParams(ParamSet &) override{return true;}
	};
	template<double (func)(ParamSet&,ParamSet&),bool (filter)(ParamSet&)>
	class Param_Func:public virtual IParamFunc{
	public:
		Param_Func(){}virtual ~Param_Func(){}
		virtual double operator()(ParamSet &X, ParamSet &P) override{return func(X,P);}
		virtual bool CorrectParams(ParamSet &P) override{return filter(P);}
	};
	
	template<int value>
	class Const:public virtual IParamFunc{
	public:
		Const(){}virtual ~Const(){}
		virtual double operator()(ParamSet &, ParamSet &) override{return double(value);}
		virtual bool CorrectParams(ParamSet &) override{return true;}
	};
	template<int x_index>
	class Arg:public virtual IParamFunc{
	public:
		Arg(){}virtual ~Arg(){}
		virtual double operator()(ParamSet &X, ParamSet &) override{return X[x_index];}
		virtual bool CorrectParams(ParamSet &) override{return true;}
	};
	template<int p_index>
	class Par:public virtual IParamFunc{
	public:
		Par(){}virtual ~Par(){}
		virtual double operator()(ParamSet &, ParamSet &P) override{return P[p_index];}
		virtual bool CorrectParams(ParamSet &) override{return true;}
	};
	template<int x_index,int p_index,unsigned int power>
	class PolynomFunc:public virtual IParamFunc{
	public:
		PolynomFunc(){}
		virtual ~PolynomFunc(){}
		virtual double operator()(ParamSet &X, ParamSet &P) override{
			return Polynom<power,double,ParamSet,p_index>(X[x_index],P);
		}
		virtual bool CorrectParams(ParamSet&) override{return true;}
	};
	template<double (func)(double), class FUNC>
	class Func:public virtual FUNC{
		public:Func():FUNC(){}virtual ~Func(){}
		virtual double operator()(ParamSet &X, ParamSet &P) override{
			return func(FUNC::operator()(X,P));
		}
		virtual bool CorrectParams(ParamSet& P) override{return FUNC::CorrectParams(P);}
	};
	template<double (func)(double,double), class FUNC1, class FUNC2>
	class Func2:public virtual FUNC1,public virtual FUNC2{
		public:Func2():FUNC1(),FUNC2(){}virtual ~Func2(){}
		virtual double operator()(ParamSet &X, ParamSet &P) override{
			return func(FUNC1::operator()(X,P),FUNC2::operator()(X,P));
		}
		virtual bool CorrectParams(ParamSet& P) override{
			return FUNC1::CorrectParams(P)&&FUNC2::CorrectParams(P);
		}
	};
	template<double (func)(double,double,double), class FUNC1, class FUNC2, class FUNC3>
	class Func3:public virtual FUNC1,public virtual FUNC2,public virtual FUNC3{
		public:Func3():FUNC1(),FUNC2(),FUNC3(){}virtual ~Func3(){}
		virtual double operator()(ParamSet &X, ParamSet &P) override{
			return func(FUNC1::operator()(X,P),FUNC2::operator()(X,P),FUNC3::operator()(X,P));
		}
		virtual bool CorrectParams(ParamSet& P) override{
			return FUNC1::CorrectParams(P)&&FUNC2::CorrectParams(P)&&FUNC3::CorrectParams(P);
		}
	};
	template<double (func)(double,double,double,double), class FUNC1, class FUNC2, class FUNC3, class FUNC4>
	class Func4:public virtual FUNC1,public virtual FUNC2,public virtual FUNC3,public virtual FUNC4{
		public:Func4():FUNC1(),FUNC2(),FUNC3(),FUNC4(){}virtual ~Func4(){}
		virtual double operator()(ParamSet &X, ParamSet &P) override{
			return func(FUNC1::operator()(X,P),FUNC2::operator()(X,P),FUNC3::operator()(X,P),FUNC4::operator()(X,P));
		}
		virtual bool CorrectParams(ParamSet& P) override{
			return FUNC1::CorrectParams(P)&&FUNC2::CorrectParams(P)&&FUNC3::CorrectParams(P)&&FUNC4::CorrectParams(P);
		}
	};
	template<double (func)(double,double,double,double,double), class FUNC1, class FUNC2, class FUNC3, class FUNC4, class FUNC5>
	class Func5:public virtual FUNC1,public virtual FUNC2,public virtual FUNC3,public virtual FUNC4,public virtual FUNC5{
		public:Func5():FUNC1(),FUNC2(),FUNC3(),FUNC4(),FUNC5(){}virtual ~Func5(){}
		virtual double operator()(ParamSet &X, ParamSet &P) override{
			return func(FUNC1::operator()(X,P),FUNC2::operator()(X,P),FUNC3::operator()(X,P),FUNC4::operator()(X,P),FUNC5::operator()(X,P));
		}
		virtual bool CorrectParams(ParamSet& P) override{
			return FUNC1::CorrectParams(P)&&FUNC2::CorrectParams(P)&&FUNC3::CorrectParams(P)&&FUNC4::CorrectParams(P)&&FUNC5::CorrectParams(P);
		}
	};
	template<double (func)(double,double,double,double,double,double), class FUNC1, class FUNC2, class FUNC3, class FUNC4, class FUNC5, class FUNC6>
	class Func6:public virtual FUNC1,public virtual FUNC2,public virtual FUNC3,public virtual FUNC4,public virtual FUNC5,public virtual FUNC6{
		public:Func6():FUNC1(),FUNC2(),FUNC3(),FUNC4(),FUNC5(),FUNC6(){}virtual ~Func6(){}
		virtual double operator()(ParamSet &X, ParamSet &P) override{
			return func(FUNC1::operator()(X,P),FUNC2::operator()(X,P),FUNC3::operator()(X,P),FUNC4::operator()(X,P),FUNC5::operator()(X,P),FUNC6::operator()(X,P));
		}
		virtual bool CorrectParams(ParamSet& P) override{
			return FUNC1::CorrectParams(P)&&FUNC2::CorrectParams(P)&&FUNC3::CorrectParams(P)&&FUNC4::CorrectParams(P)&&FUNC5::CorrectParams(P)&&FUNC6::CorrectParams(P);
		}
	};
	template<class FUNC,int x_index, int p_index>
	class ArgShift:public virtual FUNC{
	public:
		ArgShift():FUNC(){}virtual ~ArgShift(){}
		virtual double operator()(ParamSet &X, ParamSet &P) override{
			ParamSet x=X;
			x.Set(x_index,x[x_index]+P[p_index]);
			return FUNC::operator()(x,P);
		}
		virtual bool CorrectParams(ParamSet& P) override{return FUNC::CorrectParams(P);}
	};
	template<class FUNC,int x_index, int p_index>
	class ArgScale:public virtual FUNC{
	public:
		ArgScale():FUNC(){}virtual ~ArgScale(){}
		virtual double operator()(ParamSet &X, ParamSet &P) override{
			ParamSet x=X;
			x.Set(x_index,x[x_index]*P[p_index]);
			return FUNC::operator()(x,P);
		}
		virtual bool CorrectParams(ParamSet& P) override{return FUNC::CorrectParams(P);}
	};
	template<class FUNC1,class FUNC2>
	class Add:public virtual FUNC1, public virtual FUNC2{
	public:
		Add():FUNC1(),FUNC2(){}virtual ~Add(){}
		virtual double operator()(ParamSet &X, ParamSet &P) override{
			return FUNC1::operator()(X,P)+FUNC2::operator()(X,P);
		}
		virtual bool CorrectParams(ParamSet &P) override{
			return FUNC1::CorrectParams(P)&&FUNC2::CorrectParams(P);
		}
	};
	template<class FUNC1,class FUNC2>
	class Sub:public virtual FUNC1, public virtual FUNC2{
	public:
		Sub():FUNC1(),FUNC2(){}virtual ~Sub(){}
		virtual double operator()(ParamSet &X, ParamSet &P) override{
			return FUNC1::operator()(X,P)-FUNC2::operator()(X,P);
		}
		virtual bool CorrectParams(ParamSet &P) override{
			return FUNC1::CorrectParams(P)&&FUNC2::CorrectParams(P);
		}
	};
	template<class FUNC1,class FUNC2>
	class Mul:public virtual FUNC1, public virtual FUNC2{
	public:
		Mul():FUNC1(),FUNC2(){}virtual ~Mul(){}
		virtual double operator()(ParamSet &X, ParamSet &P) override{
			return FUNC1::operator()(X,P)*FUNC2::operator()(X,P);
		}
		virtual bool CorrectParams(ParamSet &P) override{
			return FUNC1::CorrectParams(P)&&FUNC2::CorrectParams(P);
		}
	};
	template<class FUNC1,class FUNC2>
	class Div:public virtual FUNC1, public virtual FUNC2{
	public:
		Div():FUNC1(),FUNC2(){}virtual ~Div(){}
		virtual double operator()(ParamSet &X, ParamSet &P) override{
			return FUNC1::operator()(X,P)/FUNC2::operator()(X,P);
		}
		virtual bool CorrectParams(ParamSet &P) override{
			return FUNC1::CorrectParams(P)&&FUNC2::CorrectParams(P);
		}
	};
	template<class FUNC1,class FUNC2>
	class Power:public virtual FUNC1, public virtual FUNC2{
	public:
		Power():FUNC1(),FUNC2(){}virtual ~Power(){}
		virtual double operator()(ParamSet &X, ParamSet &P) override{
			return pow(FUNC1::operator()(X,P),FUNC2::operator()(X,P));
		}
		virtual bool CorrectParams(ParamSet &P) override{
			return FUNC1::CorrectParams(P)&&FUNC2::CorrectParams(P);
		}
	};
}
#endif
#ifndef ____mVYpymgQ
#define ____mVYpymgQ
#include <math.h>
#include "fit_gen.h"
#include "math_h/functions.h"
namespace Fit{
	using namespace std;
	template<int a,int b>
	struct max2{enum{val=(a>b)?a:b};};
	template<int a,int b,int c>
	struct max3{enum{val=(max2<a,b>::val>c)?max2<a,b>::val:c};};
	template<int a,int b,int c,int d>
	struct max4{enum{val=(max3<a,b,c>::val>d)?max3<a,b,c>::val:d};};
	template<int a,int b,int c,int d,int e>
	struct max5{enum{val=(max4<a,b,c,d>::val>e)?max4<a,b,c,d>::val:e};};
	template<int a,int b,int c,int d,int e,int f>
	struct max6{enum{val=(max5<a,b,c,d,e>::val>f)?max5<a,b,c,d,e>::val:f};};
	
	template<int value>
	class Const:public virtual IParamFunc{
	public:
		Const(){}
		virtual ~Const(){}
		virtual double operator()(ParamSet, ParamSet) override{
			return double(value);
		}
		virtual bool CorrectParams(ParamSet) override{
			return true;
		}
		enum{ParamCount=0,ArgCount=0};
	};
	template<int x_index>
	class Arg:public virtual IParamFunc{
	public:
		Arg(){}
		virtual ~Arg(){}
		virtual double operator()(ParamSet X, ParamSet) override{
			return X[x_index];
		}
		virtual bool CorrectParams(ParamSet) override{
			return true;
		}
		enum{ParamCount=0,ArgCount=x_index+1};
	};
	template<int p_index>
	class Par:public virtual IParamFunc{
	public:
		Par(){}
		virtual ~Par(){}
		virtual double operator()(ParamSet, ParamSet P) override{
			return P[p_index];
		}
		virtual bool CorrectParams(ParamSet) override{
			return true;
		}
		enum{ParamCount=p_index+1,ArgCount=0};
	};
	template<int x_index,int p_index,unsigned int power>
	class PolynomFunc:public virtual IParamFunc{
	public:
		PolynomFunc(){}
		virtual ~PolynomFunc(){}
		virtual double operator()(ParamSet X, ParamSet P) override{
			return Polynom<power,double,ParamSet,p_index>(X[x_index],P);
		}
		virtual bool CorrectParams(ParamSet) override{
			return true;
		}
		enum{ParamCount=p_index+power+1,ArgCount=x_index+1};
	};
	template<double (func)(double), class FUNC>
	class Func:public virtual FUNC{
		public:Func():FUNC(){}
		virtual ~Func(){}
		virtual double operator()(ParamSet X, ParamSet P) override{
			return func(FUNC::operator()(X,P));
		}
		virtual bool CorrectParams(ParamSet P) override{
			return FUNC::CorrectParams(P);
		}
		enum{ParamCount=FUNC::ParamCount,ArgCount=FUNC::ArgCount};
	};
	template<double (func)(double,double), class FUNC1, class FUNC2>
	class Func2:public virtual FUNC1,public virtual FUNC2{
		public:Func2():FUNC1(),FUNC2(){}
		virtual ~Func2(){}
		virtual double operator()(ParamSet X, ParamSet P) override{
			return func(FUNC1::operator()(X,P),FUNC2::operator()(X,P));
		}
		virtual bool CorrectParams(ParamSet P) override{
			return FUNC1::CorrectParams(P)&&FUNC2::CorrectParams(P);
		}
		enum{
			ParamCount=max2<FUNC1::ParamCount,FUNC2::ParamCount>::val,
			ArgCount=max2<FUNC1::ArgCount,FUNC2::ArgCount>::val
		};
	};
	template<double (func)(double,double,double), class FUNC1, class FUNC2, class FUNC3>
	class Func3:public virtual FUNC1,public virtual FUNC2,public virtual FUNC3{
		public:Func3():FUNC1(),FUNC2(),FUNC3(){}
		virtual ~Func3(){}
		virtual double operator()(ParamSet X, ParamSet P) override{
			return func(FUNC1::operator()(X,P),FUNC2::operator()(X,P),FUNC3::operator()(X,P));
		}
		virtual bool CorrectParams(ParamSet P) override{
			return FUNC1::CorrectParams(P)&&FUNC2::CorrectParams(P)&&FUNC3::CorrectParams(P);
		}
		enum{
			ParamCount=max3<FUNC1::ParamCount,FUNC2::ParamCount,FUNC3::ParamCount>::val,
			ArgCount=max3<FUNC1::ArgCount,FUNC2::ArgCount,FUNC3::ArgCount>::val
		};
	};
	template<double (func)(double,double,double,double), class FUNC1, class FUNC2, class FUNC3, class FUNC4>
	class Func4:public virtual FUNC1,public virtual FUNC2,public virtual FUNC3,public virtual FUNC4{
		public:Func4():FUNC1(),FUNC2(),FUNC3(),FUNC4(){}
		virtual ~Func4(){}
		virtual double operator()(ParamSet X, ParamSet P) override{
			return func(FUNC1::operator()(X,P),FUNC2::operator()(X,P),FUNC3::operator()(X,P),FUNC4::operator()(X,P));
		}
		virtual bool CorrectParams(ParamSet P) override{
			return FUNC1::CorrectParams(P)&&FUNC2::CorrectParams(P)&&FUNC3::CorrectParams(P)&&FUNC4::CorrectParams(P);
		}
		enum{
			ParamCount=max4<FUNC1::ParamCount,FUNC2::ParamCount,FUNC3::ParamCount,FUNC4::ParamCount>::val,
			ArgCount=max4<FUNC1::ArgCount,FUNC2::ArgCount,FUNC3::ArgCount,FUNC4::ArgCount>::val
		};
	};
	template<double (func)(double,double,double,double,double), class FUNC1, class FUNC2, class FUNC3, class FUNC4, class FUNC5>
	class Func5:public virtual FUNC1,public virtual FUNC2,public virtual FUNC3,public virtual FUNC4,public virtual FUNC5{
		public:Func5():FUNC1(),FUNC2(),FUNC3(),FUNC4(),FUNC5(){}
		virtual ~Func5(){}
		virtual double operator()(ParamSet X, ParamSet P) override{
			return func(FUNC1::operator()(X,P),FUNC2::operator()(X,P),FUNC3::operator()(X,P),FUNC4::operator()(X,P),FUNC5::operator()(X,P));
		}
		virtual bool CorrectParams(ParamSet P) override{
			return FUNC1::CorrectParams(P)&&FUNC2::CorrectParams(P)&&FUNC3::CorrectParams(P)&&FUNC4::CorrectParams(P)&&FUNC5::CorrectParams(P);
		}
		enum{
			ParamCount=max5<FUNC1::ParamCount,FUNC2::ParamCount,FUNC3::ParamCount,FUNC4::ParamCount,FUNC5::ParamCount>::val,
			ArgCount=max5<FUNC1::ArgCount,FUNC2::ArgCount,FUNC3::ArgCount,FUNC4::ArgCount,FUNC5::ArgCount>::val
		};
	};
	template<double (func)(double,double,double,double,double,double), class FUNC1, class FUNC2, class FUNC3, class FUNC4, class FUNC5, class FUNC6>
	class Func6:public virtual FUNC1,public virtual FUNC2,public virtual FUNC3,public virtual FUNC4,public virtual FUNC5,public virtual FUNC6{
		public:Func6():FUNC1(),FUNC2(),FUNC3(),FUNC4(),FUNC5(),FUNC6(){}
		virtual ~Func6(){}
		virtual double operator()(ParamSet X, ParamSet P) override{
			return func(FUNC1::operator()(X,P),FUNC2::operator()(X,P),FUNC3::operator()(X,P),FUNC4::operator()(X,P),FUNC5::operator()(X,P),FUNC6::operator()(X,P));
		}
		virtual bool CorrectParams(ParamSet P) override{
			return FUNC1::CorrectParams(P)&&FUNC2::CorrectParams(P)&&FUNC3::CorrectParams(P)&&FUNC4::CorrectParams(P)&&FUNC5::CorrectParams(P)&&FUNC6::CorrectParams(P);
		}
		enum{
			ParamCount=max6<FUNC1::ParamCount,FUNC2::ParamCount,FUNC3::ParamCount,FUNC4::ParamCount,FUNC5::ParamCount,FUNC6::ParamCount>::val,
			ArgCount=max6<FUNC1::ArgCount,FUNC2::ArgCount,FUNC3::ArgCount,FUNC4::ArgCount,FUNC5::ArgCount,FUNC6::ArgCount>::val
		};
	};
	template<class FUNC,int x_index, int p_index>
	class ArgShift:public virtual FUNC{
	public:
		ArgShift():FUNC(){}virtual ~ArgShift(){}
		virtual double operator()(ParamSet X, ParamSet P) override{
			ParamSet x=X;
			x.Set(x_index,x[x_index]+P[p_index]);
			return FUNC::operator()(x,P);
		}
		virtual bool CorrectParams(ParamSet P) override{
			return FUNC::CorrectParams(P);
		}
		enum{ParamCount=FUNC::ParamCount,ArgCount=FUNC::ArgCount};
	};
	template<class FUNC,int x_index, int p_index>
	class ArgScale:public virtual FUNC{
	public:
		ArgScale():FUNC(){}
		virtual ~ArgScale(){}
		virtual double operator()(ParamSet X, ParamSet P) override{
			ParamSet x=X;
			x.Set(x_index,x[x_index]*P[p_index]);
			return FUNC::operator()(x,P);
		}
		virtual bool CorrectParams(ParamSet P) override{
			return FUNC::CorrectParams(P);
		}
		enum{ParamCount=FUNC::ParamCount,ArgCount=FUNC::ArgCount};
	};
	template<class FUNC>
	class Minus:public virtual FUNC{
	public:
		Minus():FUNC(){}
		virtual ~Minus(){}
		virtual double operator()(ParamSet X, ParamSet P) override{
			return -FUNC::operator()(X,P);
		}
		virtual bool CorrectParams(ParamSet P) override{
			return FUNC::CorrectParams(P);
		}
		enum{
			ParamCount=FUNC::ParamCount,
			ArgCount=FUNC::ArgCount
		};
	};
	template<class FUNC1,class FUNC2>
	class Add:public virtual FUNC1, public virtual FUNC2{
	public:
		Add():FUNC1(),FUNC2(){}
		virtual ~Add(){}
		virtual double operator()(ParamSet X, ParamSet P) override{
			return FUNC1::operator()(X,P)+FUNC2::operator()(X,P);
		}
		virtual bool CorrectParams(ParamSet P) override{
			return FUNC1::CorrectParams(P)&&FUNC2::CorrectParams(P);
		}
		enum{
			ParamCount=max2<FUNC1::ParamCount,FUNC2::ParamCount>::val,
			ArgCount=max2<FUNC1::ArgCount,FUNC2::ArgCount>::val
		};
	};
	template<class FUNC1,class FUNC2>
	class Sub:public virtual FUNC1, public virtual FUNC2{
	public:
		Sub():FUNC1(),FUNC2(){}
		virtual ~Sub(){}
		virtual double operator()(ParamSet X, ParamSet P) override{
			return FUNC1::operator()(X,P)-FUNC2::operator()(X,P);
		}
		virtual bool CorrectParams(ParamSet P) override{
			return FUNC1::CorrectParams(P)&&FUNC2::CorrectParams(P);
		}
		enum{
			ParamCount=max2<FUNC1::ParamCount,FUNC2::ParamCount>::val,
			ArgCount=max2<FUNC1::ArgCount,FUNC2::ArgCount>::val
		};
	};
	template<class FUNC1,class FUNC2>
	class Mul:public virtual FUNC1, public virtual FUNC2{
	public:
		Mul():FUNC1(),FUNC2(){}
		virtual ~Mul(){}
		virtual double operator()(ParamSet X, ParamSet P) override{
			return FUNC1::operator()(X,P)*FUNC2::operator()(X,P);
		}
		virtual bool CorrectParams(ParamSet P) override{
			return FUNC1::CorrectParams(P)&&FUNC2::CorrectParams(P);
		}
		enum{
			ParamCount=max2<FUNC1::ParamCount,FUNC2::ParamCount>::val,
			ArgCount=max2<FUNC1::ArgCount,FUNC2::ArgCount>::val
		};
	};
	template<class FUNC1,class FUNC2>
	class Div:public virtual FUNC1, public virtual FUNC2{
	public:
		Div():FUNC1(),FUNC2(){}
		virtual ~Div(){}
		virtual double operator()(ParamSet X, ParamSet P) override{
			return FUNC1::operator()(X,P)/FUNC2::operator()(X,P);
		}
		virtual bool CorrectParams(ParamSet P) override{
			return FUNC1::CorrectParams(P)&&FUNC2::CorrectParams(P);
		}
		enum{
			ParamCount=max2<FUNC1::ParamCount,FUNC2::ParamCount>::val,
			ArgCount=max2<FUNC1::ArgCount,FUNC2::ArgCount>::val
		};
	};
	template<class FUNC1,class FUNC2>
	class Power:public virtual FUNC1, public virtual FUNC2{
	public:
		Power():FUNC1(),FUNC2(){}
		virtual ~Power(){}
		virtual double operator()(ParamSet X, ParamSet P) override{
			return pow(FUNC1::operator()(X,P),FUNC2::operator()(X,P));
		}
		virtual bool CorrectParams(ParamSet P) override{
			return FUNC1::CorrectParams(P)&&FUNC2::CorrectParams(P);
		}
		enum{
			ParamCount=max2<FUNC1::ParamCount,FUNC2::ParamCount>::val,
			ArgCount=max2<FUNC1::ArgCount,FUNC2::ArgCount>::val
		};
	};
}
#endif
#ifndef RGCFVQNCAQKWQXFU
#define RGCFVQNCAQKWQXFU
#include "fit_gen.h"
#include "fitexception.h"
#include "math_h/randomfunc.h"
namespace Fit{
	using namespace std;
	template<class FITGEN>
	class Crossing:public FITGEN{
	private:
		double P;
	public:
		Crossing(shared_ptr<IParamFunc> function, shared_ptr<IOptimalityFunction> optimality):
			FITGEN(function,optimality){
				P=0;
		}
		virtual ~Crossing(){}
		double CrossingProbability(){
			return P;
		}
		void SetCrossingProbability(double val){
			if((val<0)||(val>1))
				throw new FitException("Invalid crossing probability value");
			P=val;
		}
	protected:
		virtual ParamSet born(ParamSet C)override{
			auto X=FITGEN::born(C);
			if(P>0){
				auto C=FITGEN::born(this->Parameters(rand()%this->PopulationSize()));
				for(int i=0; i<this->ParamCount();i++)
					if(RandomUniformly(0.0,1.0)<P)
						X.Set(i,C[i]);
			}
			return X;
		}
	};
	template<class FITGEN>
	class AbsoluteRandomMutations:public FITGEN{
	private:
		ParamSet P;
	public:
		AbsoluteRandomMutations(shared_ptr<IParamFunc> function, shared_ptr<IOptimalityFunction> optimality):
			FITGEN(function,optimality){}
		virtual ~AbsoluteRandomMutations(){}
		ParamSet AbsoluteMutations(){
			return P;
		}
		void SetAbsoluteMutations(ParamSet p){
			for(double v:p)
				if(v<0)
					throw new FitException("Wrong random mutation coefficient. Cannot be negative");
			P=p;
		}
	protected:
		virtual ParamSet born(ParamSet C)override{
			auto res=FITGEN::born(C);
			for(int i=0;i<res.Count();i++){
				res.Set(i,res[i]+RandomGauss(P[i]));
			}
			return res;
		}
	};
	template<class FITGEN>
	class RelativeRandomMutations:public FITGEN{
	private:
		ParamSet P;
	public:
		RelativeRandomMutations(shared_ptr<IParamFunc> function, shared_ptr<IOptimalityFunction> optimality):
			FITGEN(function,optimality){}
		virtual ~RelativeRandomMutations(){}
		ParamSet RelativeMutations(){
			return P;
		}
		void SetRelativeMutations(ParamSet p){
			for(double v:p)
				if(v<0)
					throw new FitException("Wrong random mutation coefficient. Cannot be negative");
				P=p;
		}
	protected:
		virtual ParamSet born(ParamSet C)override{
			auto res=FITGEN::born(C);
			for(int i=0;i<res.Count();i++){
				res.Set(i,res[i]+RandomGauss(P[i]*C[i]));
			}
			return res;
		}
	};
}
#endif
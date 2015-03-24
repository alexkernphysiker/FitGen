#ifndef RGCFVQNCAQKWQXFU
#define RGCFVQNCAQKWQXFU
#include "fit_gen.h"
#include "fitexception.h"
#include "math_h/randomfunc.h"
namespace Fit{
	using namespace std;
	template<class FITGEN=AbstractGenetic>
	class DifferentialRandomMutations: public FITGEN{
	private:
		double M;
	public:
		DifferentialRandomMutations(shared_ptr<IParamFunc> function, 
			shared_ptr<IOptimalityFunction> optimality, 
			unsigned int threads_count
		):FITGEN(function,optimality,threads_count),M(0.5){}
		virtual ~DifferentialRandomMutations(){}
		double MutationCoefficient(){
			return M;
		}
		void SetMutationCoefficient(double val){
			if(val<=0)
				throw new FitException("Invalid mutation coefficient value");
			M=val;
		}
	protected:
		virtual void mutations(ParamSet &C)override{
			FITGEN::mutations(C);
			auto A=AbstractGenetic::Parameters(rand()%AbstractGenetic::PopulationSize());
			auto B=AbstractGenetic::Parameters(rand()%AbstractGenetic::PopulationSize());
			for(int i=0; i<AbstractGenetic::ParamCount();i++)
				C.Set(i,C[i]+M*(A[i]-B[i]));
		}
	};
	template<class FITGEN=AbstractGenetic>
	class Crossing:public FITGEN{
	private:
		double P;
	public:
		Crossing(shared_ptr<IParamFunc> function, 
			shared_ptr<IOptimalityFunction> optimality,
			unsigned int threads_count
		):FITGEN(function,optimality,threads_count),P(0.1){}
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
		virtual void mutations(ParamSet &C)override{
			FITGEN::mutations(C);
			if(P>0){
				auto X=AbstractGenetic::Parameters(rand()%AbstractGenetic::PopulationSize());
				FITGEN::mutations(X);
				for(int i=0; i<AbstractGenetic::ParamCount();i++)
					if(RandomUniformly(0.0,1.0)<P)
						C.Set(i,X[i]);
			}
		}
	};
	template<class FITGEN=AbstractGenetic>
	class AbsoluteRandomMutations:public FITGEN{
	private:
		ParamSet M;
	public:
		AbsoluteRandomMutations(shared_ptr<IParamFunc> function, 
			shared_ptr<IOptimalityFunction> optimality,
			unsigned int threads_count
		):FITGEN(function,optimality,threads_count){}
		virtual ~AbsoluteRandomMutations(){}
		ParamSet AbsoluteMutations(){
			return M;
		}
		void SetAbsoluteMutations(ParamSet p){
			for(double v:p)
				if(v<0)
					throw new FitException("Wrong random mutation coefficient. Cannot be negative");
			M=p;
		}
	protected:
		virtual void mutations(ParamSet &C)override{
			FITGEN::mutations(C);
			for(int i=0;i<AbstractGenetic::ParamCount();i++){
				C.Set(i,C[i]+RandomGauss(M[i]));
			}
		}
	};
	template<class FITGEN=AbstractGenetic>
	class RelativeRandomMutations:public FITGEN{
	private:
		ParamSet M;
	public:
		RelativeRandomMutations(shared_ptr<IParamFunc> function, 
			shared_ptr<IOptimalityFunction> optimality,
			unsigned int threads_count
		):FITGEN(function,optimality,threads_count){}
		virtual ~RelativeRandomMutations(){}
		ParamSet RelativeMutations(){
			return M;
		}
		void SetRelativeMutations(ParamSet p){
			for(double v:p)
				if(v<0)
					throw new FitException("Wrong random mutation coefficient. Cannot be negative");
				M=p;
		}
	protected:
		virtual void mutations(ParamSet &C)override{
			FITGEN::mutations(C);
			for(int i=0;i<AbstractGenetic::ParamCount();i++){
				C.Set(i,C[i]*(1+RandomGauss(M[i])));
			}
		}
	};
}
#endif
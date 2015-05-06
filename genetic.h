#ifndef RGCFVQNCAQKWQXFU
#define RGCFVQNCAQKWQXFU
#include "abstract.h"
#include "genetic_exception.h"
#include "math_h/randomfunc.h"
namespace Genetic{
	using namespace std;
	template<class FITGEN=AbstractGenetic>
	class DifferentialMutations: public FITGEN{
	private:
		double M;
	public:
		DifferentialMutations(shared_ptr<IOptimalityFunction> optimality)
			:FITGEN(optimality),M(0.5){}
		virtual ~DifferentialMutations(){}
		double MutationCoefficient(){
			return M;
		}
		void SetMutationCoefficient(double val){
			if(val<=0)
				throw new GeneticException("Invalid mutation coefficient value");
			M=val;
		}
	protected:
		virtual void mutations(ParamSet &C)override{
			FITGEN::mutations(C);
			auto A=AbstractGenetic::Parameters(RandomUniformlyI(0,AbstractGenetic::PopulationSize()-1));
			auto B=AbstractGenetic::Parameters(RandomUniformlyI(0,AbstractGenetic::PopulationSize()-1));
			for(int i=0; i<C.Count();i++)
				C.Set(i,C[i]+M*(A[i]-B[i]));
		}
	};
	template<class FITGEN=AbstractGenetic>
	class Crossing:public FITGEN{
	private:
		double P;
	public:
		Crossing(shared_ptr<IOptimalityFunction> optimality)
			:FITGEN(optimality),P(0){}
		virtual ~Crossing(){}
		double CrossingProbability(){
			return P;
		}
		void SetCrossingProbability(double val){
			if((val<0)||(val>1))
				throw new GeneticException("Invalid crossing probability value");
			P=val;
		}
	protected:
		virtual void mutations(ParamSet &C)override{
			FITGEN::mutations(C);
			if(RandomUniformlyR(0.0,1.0)<P){
				auto X=AbstractGenetic::Parameters(RandomUniformlyI(0,AbstractGenetic::PopulationSize()-1));
				FITGEN::mutations(X);
				for(int i=0; i<C.Count();i++)
					if(rand()%2==1)
						C.Set(i,X[i]);
			}
		}
	};
	template<class FITGEN=AbstractGenetic>
	class AbsoluteMutations:public FITGEN{
	private:
		ParamSet M;
		double P;
	public:
		AbsoluteMutations(shared_ptr<IOptimalityFunction> optimality)
			:FITGEN(optimality),P(0){}
		virtual ~AbsoluteMutations(){}
		ParamSet AbsoluteMutationCoeficients(){
			return M;
		}
		void SetAbsoluteMutationCoeficients(ParamSet p){
			for(double v:p)
				if(v<0)
					throw new GeneticException("Wrong random mutation coefficient. Cannot be negative");
			M=p;
		}
		double AbsoluteMutationsProbability(){
			return P;
		}
		void SetAbsoluteMutationsProbability(double val){
			if((val<0)||(val>1))
				throw new GeneticException("Invalid crossing probability value");
			P=val;
		}
	protected:
		virtual void mutations(ParamSet &C)override{
			FITGEN::mutations(C);
			if(RandomUniformlyR(0.0,1.0)<P)
				for(int i=0;i<C.Count();i++)
					C.Set(i,C[i]+RandomGauss(M[i]));
		}
	};
	template<class FITGEN=AbstractGenetic>
	class RelativeMutations:public FITGEN{
	private:
		ParamSet M;
		double P;
	public:
		RelativeMutations(shared_ptr<IOptimalityFunction> optimality)
			:FITGEN(optimality),P(0){}
		virtual ~RelativeMutations(){}
		ParamSet RelativeMutationCoefficients(){
			return M;
		}
		void SetRelativeMutationCoefficients(ParamSet p){
			for(double v:p)
				if(v<0)
					throw new GeneticException("Wrong random mutation coefficient. Cannot be negative");
				M=p;
		}
		double RelativeMutationsProbability(){
			return P;
		}
		void SetRelativeMutationsProbability(double val){
			if((val<0)||(val>1))
				throw new GeneticException("Invalid crossing probability value");
			P=val;
		}
	protected:
		virtual void mutations(ParamSet &C)override{
			FITGEN::mutations(C);
			if(RandomUniformlyR(0.0,1.0)<P)
				for(int i=0;i<AbstractGenetic::ParamCount();i++)
					C.Set(i,C[i]*(1+RandomGauss(M[i])));
		}
	};
	template<class FITGEN>
	class ExactCopying:public FITGEN{
	private:
		double P;
	public:
		ExactCopying(shared_ptr<IOptimalityFunction> optimality)
			:FITGEN(optimality),P(0){}
		virtual ~ExactCopying(){}
		void SetExactCopyingProbability(double value){
			if((value<0)||(value>1))
				throw;
			P=value;
		}
		double ExactCopyingProbability(){
			return P;
		}
	protected:
		virtual void mutations(ParamSet &C)override{
			if(RandomUniformlyR(0.0,1.0)>=P)
				FITGEN::mutations(C);
		}
	};
}
#endif
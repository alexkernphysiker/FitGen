// this file is distributed under 
// GPL v 3.0 license
#ifndef RGCFVQNCAQKWQXFU
#define RGCFVQNCAQKWQXFU
#include <random>
#include "abstract.h"
#include "genetic_exception.h"
namespace Genetic{
	using namespace std;
	template<class FITGEN=AbstractGenetic>
	class DifferentialMutations: public virtual FITGEN{
	private:
		double M;
		std::default_random_engine generator;
	public:
		DifferentialMutations():FITGEN(),M(0.5){}
		virtual ~DifferentialMutations(){}
		double MutationCoefficient(){
			return M;
		}
		void SetMutationCoefficient(double val){
			if(val<0)
				throw GeneticException("DifferentialMutations: mutation coefficient should be a positive value");
			M=val;
		}
	protected:
		virtual void mutations(ParamSet &C)override{
			FITGEN::mutations(C);
			std::uniform_int_distribution<int> randomelement(0,AbstractGenetic::PopulationSize()-1);
			auto A=AbstractGenetic::Parameters(randomelement(generator));
			auto B=AbstractGenetic::Parameters(randomelement(generator));
			for(int i=0; i<C.Count();i++)
				C.Set(i,C[i]+M*(A[i]-B[i]));
		}
	};
	template<class FITGEN=AbstractGenetic>
	class Crossing:public virtual FITGEN{
	private:
		double P;
		std::default_random_engine generator;
	public:
		Crossing():FITGEN(),P(0){}
		virtual ~Crossing(){}
		double CrossingProbability(){
			return P;
		}
		void SetCrossingProbability(double val){
			if((val<0)||(val>1))
				throw GeneticException("Crossing: probability value should fit the condition 0<=P<=1");
			P=val;
		}
	protected:
		virtual void mutations(ParamSet &C)override{
			FITGEN::mutations(C);
			std::uniform_real_distribution<double> Prob(0,1);
			if(Prob(generator)<P){
				std::uniform_int_distribution<int> randomelement(0,AbstractGenetic::PopulationSize()-1);
				auto X=AbstractGenetic::Parameters(randomelement(generator));
				FITGEN::mutations(X);
				for(int i=0; i<C.Count();i++)
					if(Prob(generator)<0.5)
						C.Set(i,X[i]);
			}
		}
	};
	template<class FITGEN=AbstractGenetic>
	class AbsoluteMutations:public virtual FITGEN{
	private:
		double P;
		std::vector<std::normal_distribution<double>> distr;
		std::default_random_engine generator;
	public:
		AbsoluteMutations():FITGEN(),P(0){}
		virtual ~AbsoluteMutations(){}
		ParamSet AbsoluteMutationCoefficients(){
			ParamSet P;
			for(auto&d:distr)
				P<<d.stddev();
			return P;
		}
		void SetAbsoluteMutationCoefficients(ParamSet p){
			distr.clear();
			for(double v:p){
				if(v<0)
					throw GeneticException("AbsoluteMutations: mutation coefficient cannot be negative");
				distr.push_back(std::normal_distribution<double>(0,v));
			}
		}
		double AbsoluteMutationsProbability(){
			return P;
		}
		void SetAbsoluteMutationsProbability(double val){
			if((val<0)||(val>1))
				throw GeneticException("AbsoluteMutations: probability value should fit the condition 0<=P<=1");
			P=val;
		}
	protected:
		virtual void mutations(ParamSet &C)override{
			FITGEN::mutations(C);
			std::uniform_real_distribution<double> Prob(0,1);
			if(Prob(generator)<P)
				for(int i=0;i<AbstractGenetic::ParamCount();i++)
					C.Set(i,C[i]+distr[i](generator));
		}
	};
	template<class FITGEN=AbstractGenetic>
	class RelativeMutations:public virtual FITGEN{
	private:
		double P;
		std::vector<std::normal_distribution<double>> distr;
		std::default_random_engine generator;
	public:
		RelativeMutations():FITGEN(),P(0){}
		virtual ~RelativeMutations(){}
		ParamSet RelativeMutationCoefficients(){
			ParamSet P;
			for(auto&d:distr)
				P<<d.stddev();
			return P;
		}
		void SetRelativeMutationCoefficients(ParamSet p){
			distr.clear();
			for(double v:p){
				if(v<0)
					throw GeneticException("AbsoluteMutations: mutation coefficient cannot be negative");
				distr.push_back(std::normal_distribution<double>(0,v));
			}
		}
		double RelativeMutationsProbability(){
			return P;
		}
		void SetRelativeMutationsProbability(double val){
			if((val<0)||(val>1))
				throw GeneticException("RelativeMutations: probability value should fit the condition 0<=P<=1");
			P=val;
		}
	protected:
		virtual void mutations(ParamSet &C)override{
			FITGEN::mutations(C);
			std::uniform_real_distribution<double> Prob(0,1);
			if(Prob(generator)<P)
				for(int i=0;i<AbstractGenetic::ParamCount();i++)
					C.Set(i,C[i]*(1+distr[i](generator)));
		}
	};
	template<class FITGEN>
	class ExactCopying:public virtual FITGEN{
	private:
		double P;
		std::default_random_engine generator;
	public:
		ExactCopying():FITGEN(),P(0){}
		virtual ~ExactCopying(){}
		void SetExactCopyingProbability(double value){
			if((value<0)||(value>1))
				throw GeneticException("ExactCopying: probability value should fit the condition 0<=P<=1");
			P=value;
		}
		double ExactCopyingProbability(){
			return P;
		}
	protected:
		virtual void mutations(ParamSet &C)override{
			std::uniform_real_distribution<double> Prob(0,1);
			if(Prob(generator)>P)
				FITGEN::mutations(C);
		}
	};
}
#endif
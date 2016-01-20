// this file is distributed under 
// MIT license
#ifndef RGCFVQNCAQKWQXFU
#define RGCFVQNCAQKWQXFU
#include <random>
#include "../math_h/error.h"
#include "abstract.h"
namespace Genetic{
	using namespace std;
	using namespace MathTemplates;
	template<class FITGEN=AbstractGenetic>
	class DifferentialMutations: public virtual FITGEN{
	private:
		double M;
	public:
		DifferentialMutations():FITGEN(),M(0.5){}
		virtual ~DifferentialMutations(){}
		double MutationCoefficient(){
			return M;
		}
		void SetMutationCoefficient(double val){
			if(val<0)
				throw Exception<DifferentialMutations>("DifferentialMutations: mutation coefficient should be a positive value");
			M=val;
		}
	protected:
		virtual void mutations(ParamSet &C,RANDOM&R)override{
			FITGEN::mutations(C,R);
			std::uniform_int_distribution<int> randomelement(0,AbstractGenetic::PopulationSize()-1);
			auto A=AbstractGenetic::Parameters(randomelement(R));
			auto B=AbstractGenetic::Parameters(randomelement(R));
			for(int i=0; i<C.size();i++)
				C[i]+=M*(A[i]-B[i]);
		}
	};
	template<class FITGEN=AbstractGenetic>
	class Crossing:public virtual FITGEN{
	private:
		double P;
	public:
		Crossing():FITGEN(),P(0){}
		virtual ~Crossing(){}
		double CrossingProbability(){
			return P;
		}
		void SetCrossingProbability(double val){
			if((val<0)||(val>1))
				throw Exception<Crossing>("Crossing: probability value should fit the condition 0<=P<=1");
			P=val;
		}
	protected:
		virtual void mutations(ParamSet &C,RANDOM&R)override{
			FITGEN::mutations(C,R);
			std::uniform_real_distribution<double> Prob(0,1);
			if(Prob(R)<P){
				std::uniform_int_distribution<int> randomelement(0,AbstractGenetic::PopulationSize()-1);
				auto X=AbstractGenetic::Parameters(randomelement(R));
				FITGEN::mutations(X,R);
				for(int i=0; i<C.size();i++)
					if(Prob(R)<0.5)
						C[i]=X[i];
			}
		}
	};
	template<class FITGEN=AbstractGenetic>
	class AbsoluteMutations:public virtual FITGEN{
	private:
		double P;
		std::vector<std::normal_distribution<double>> distr;
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
					throw Exception<AbsoluteMutations>("AbsoluteMutations: mutation coefficient cannot be negative");
				distr.push_back(std::normal_distribution<double>(0,v));
			}
		}
		double AbsoluteMutationsProbability(){
			return P;
		}
		void SetAbsoluteMutationsProbability(double val){
			if((val<0)||(val>1))
				throw Exception<AbsoluteMutations>("AbsoluteMutations: probability value should fit the condition 0<=P<=1");
			P=val;
		}
	protected:
		virtual void mutations(ParamSet &C,RANDOM&R)override{
			FITGEN::mutations(C,R);
			std::uniform_real_distribution<double> Prob(0,1);
			if(Prob(R)<P)
				for(int i=0;i<AbstractGenetic::ParamCount();i++)
					C[i]+=distr[i](R);
		}
	};
	template<class FITGEN=AbstractGenetic>
	class RelativeMutations:public virtual FITGEN{
	private:
		double P;
		std::vector<std::normal_distribution<double>> distr;
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
					throw Exception<RelativeMutations>("AbsoluteMutations: mutation coefficient cannot be negative");
				distr.push_back(std::normal_distribution<double>(0,v));
			}
		}
		double RelativeMutationsProbability(){
			return P;
		}
		void SetRelativeMutationsProbability(double val){
			if((val<0)||(val>1))
				throw Exception<RelativeMutations>("RelativeMutations: probability value should fit the condition 0<=P<=1");
			P=val;
		}
	protected:
		virtual void mutations(ParamSet &C,RANDOM&R)override{
			FITGEN::mutations(C,R);
			std::uniform_real_distribution<double> Prob(0,1);
			if(Prob(R)<P)
				for(int i=0;i<AbstractGenetic::ParamCount();i++)
					C[i]*=(1+distr[i](R));
		}
	};
	template<class FITGEN>
	class ExactCopying:public virtual FITGEN{
	private:
		double P;
	public:
		ExactCopying():FITGEN(),P(0){}
		virtual ~ExactCopying(){}
		void SetExactCopyingProbability(double value){
			if((value<0)||(value>1))
				throw Exception<ExactCopying>("ExactCopying: probability value should fit the condition 0<=P<=1");
			P=value;
		}
		double ExactCopyingProbability(){
			return P;
		}
	protected:
		virtual void mutations(ParamSet &C,RANDOM&R)override{
			std::uniform_real_distribution<double> Prob(0,1);
			if(Prob(R)>P)
				FITGEN::mutations(C,R);
		}
	};
}
#endif
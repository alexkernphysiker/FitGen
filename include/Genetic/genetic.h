// this file is distributed under 
// MIT license
#ifndef RGCFVQNCAQKWQXFU
#define RGCFVQNCAQKWQXFU
#include <random>
#include <math_h/error.h>
#include "abstract.h"
namespace Genetic{
	template<class FITGEN=AbstractGenetic>
	class DifferentialMutations: public virtual FITGEN{
	private:
		double M;
	public:
		DifferentialMutations():FITGEN(),M(0.5){}
		virtual ~DifferentialMutations(){}
		const double&MutationCoefficient()const {return M;}
		void SetMutationCoefficient(const double val){
			if(val<0)
				throw MathTemplates::Exception<DifferentialMutations>("DifferentialMutations: mutation coefficient should be a positive value");
			M=val;
		}
	protected:
		virtual void mutations(ParamSet &C,RANDOM&R)const override{
			FITGEN::mutations(C,R);
			std::uniform_int_distribution<int> randomelement(0,AbstractGenetic::PopulationSize()-1);
			auto A=AbstractGenetic::Parameters(randomelement(R));
			auto B=AbstractGenetic::Parameters(randomelement(R));
			for(int i=0; i<C.size();i++)
				C(i)+=M*(A[i]-B[i]);
		}
	};
	template<class FITGEN=AbstractGenetic>
	class Crossing:public virtual FITGEN{
	private:
		double P;
	public:
		Crossing():FITGEN(),P(0){}
		virtual ~Crossing(){}
		const double&CrossingProbability()const {return P;}
		void SetCrossingProbability(const double val){
			if((val<0)||(val>1))
				throw MathTemplates::Exception<Crossing>("Crossing: probability value should fit the condition 0<=P<=1");
			P=val;
		}
	protected:
		virtual void mutations(ParamSet &C,RANDOM&R)const override{
			FITGEN::mutations(C,R);
			std::uniform_real_distribution<double> Prob(0,1);
			if(Prob(R)<P){
				std::uniform_int_distribution<int> randomelement(0,AbstractGenetic::PopulationSize()-1);
				auto X=AbstractGenetic::Parameters(randomelement(R));
				FITGEN::mutations(X,R);
				for(size_t i=0; i<C.size();i++)
					if(Prob(R)<0.5)
						C(i)=X[i];
			}
		}
	};
	template<class FITGEN=AbstractGenetic>
	class AbsoluteMutations:public virtual FITGEN{
	private:
		double P;
		ParamSet m_mutation;
	public:
		AbsoluteMutations():FITGEN(),P(0){}
		virtual ~AbsoluteMutations(){}
		const ParamSet&AbsoluteMutationCoefficients()const {return m_mutation;}
		void SetAbsoluteMutationCoefficients(const ParamSet&&p){
			m_mutation={};
			for(double v:p){
				if(v<0)
					throw MathTemplates::Exception<AbsoluteMutations>("AbsoluteMutations: mutation coefficient cannot be negative");
				m_mutation<<v;
			}
		}
		const double&AbsoluteMutationsProbability()const {return P;}
		void SetAbsoluteMutationsProbability(const double val){
			if((val<0)||(val>1))
				throw MathTemplates::Exception<AbsoluteMutations>("AbsoluteMutations: probability value should fit the condition 0<=P<=1");
			P=val;
		}
	protected:
		virtual void mutations(ParamSet &C,RANDOM&R)const override{
			FITGEN::mutations(C,R);
			std::uniform_real_distribution<double> Prob(0,1);
			if(Prob(R)<P)
				for(size_t i=0;i<AbstractGenetic::ParamCount();i++){
					std::normal_distribution<double>distr(0,m_mutation[i]);
					C(i)+=distr(R);
				}
		}
	};
	template<class FITGEN=AbstractGenetic>
	class RelativeMutations:public virtual FITGEN{
	private:
		double P;
		ParamSet m_mutation;
	public:
		RelativeMutations():FITGEN(),P(0){}
		virtual ~RelativeMutations(){}
		const ParamSet&RelativeMutationCoefficients()const {return m_mutation;}
		void SetRelativeMutationCoefficients(const ParamSet&&p){
			m_mutation={};
			for(double v:p){
				if(v<0)
					throw MathTemplates::Exception<RelativeMutations>("AbsoluteMutations: mutation coefficient cannot be negative");
				m_mutation<<v;
			}
		}
		const double&RelativeMutationsProbability()const {return P;}
		void SetRelativeMutationsProbability(const double val){
			if((val<0)||(val>1))
				throw MathTemplates::Exception<RelativeMutations>("RelativeMutations: probability value should fit the condition 0<=P<=1");
			P=val;
		}
	protected:
		virtual void mutations(ParamSet &C,RANDOM&R)const override{
			FITGEN::mutations(C,R);
			std::uniform_real_distribution<double> Prob(0,1);
			if(Prob(R)<P)
				for(size_t i=0;i<AbstractGenetic::ParamCount();i++){
					std::normal_distribution<double> distr(0,m_mutation[i]);
					C(i)*=(1+distr(R));
				}
		}
	};
	template<class FITGEN>
	class ExactCopying:public virtual FITGEN{
	private:
		double P;
	public:
		ExactCopying():FITGEN(),P(0){}
		virtual ~ExactCopying(){}
		void SetExactCopyingProbability(const double value){
			if((value<0)||(value>1))
				throw MathTemplates::Exception<ExactCopying>("ExactCopying: probability value should fit the condition 0<=P<=1");
			P=value;
		}
		const double&ExactCopyingProbability()const {return P;}
	protected:
		virtual void mutations(ParamSet &C,RANDOM&R)const override{
			std::uniform_real_distribution<double> Prob(0,1);
			if(Prob(R)>P)
				FITGEN::mutations(C,R);
		}
	};
}
#endif
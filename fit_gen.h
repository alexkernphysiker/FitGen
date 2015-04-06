#ifndef ____WrKDhKHP___
#define ____WrKDhKHP___
#include <mutex>
#include <vector>
#include <memory>
#include <utility>
#include <math.h>
#include "paramset.h"
namespace Fit{
	using namespace std;
	class IInitialConditions{
	public:
		virtual ~IInitialConditions(){}
		virtual ParamSet Generate()=0;
	};
	class IParamCheck{
	public:
		virtual ~IParamCheck(){}
		virtual bool CorrectParams(ParamSet&P)=0;
	};
	class IParamFunc:public IParamCheck{
	public:
		virtual ~IParamFunc(){}
		virtual double operator()(ParamSet&X,ParamSet&P)=0;
	};
	class IOptimalityFunction{
	public:
		virtual ~IOptimalityFunction(){}
		virtual double operator()(ParamSet&P,IParamFunc&F)=0;
	};
	class AbstractGenetic{
	protected:
		AbstractGenetic(shared_ptr<IParamFunc> function, 
			shared_ptr<IOptimalityFunction> optimality,
			unsigned int threads_count);
	public:
		virtual ~AbstractGenetic();
		
		shared_ptr<IParamFunc> Function();
		shared_ptr<IOptimalityFunction> OptimalityCalculator();
		void SetFilter(shared_ptr<IParamCheck> filter);
		void RemoveFilter();
		
		void Init(int population_size,shared_ptr<IInitialConditions> initial_conditions);
	protected:
		virtual void mutations(ParamSet&);
	public:
		void Iterate();
		
		unsigned int iteration_count();
		int PopulationSize();
		int ParamCount();
		double Optimality(int point_index=0);
		ParamSet& Parameters(int point_index=0);
		double operator[](int i);
		double operator()(ParamSet X);
		ParamSet ParamAverage();
		ParamSet ParamDispersion();
		ParamSet ParamMaxDeviation();
		
		bool ConcentratedInOnePoint();
		bool AbsoluteOptimalityExitCondition(double accuracy);
		bool RelativeOptimalityExitCondition(double accuracy);
		
		ParamSet GetParamParabolicError(ParamSet delta);
		
		typedef vector<double>::iterator iterator;
		typedef vector<double>::const_iterator const_iterator;
		iterator begin();
		const_iterator cbegin()const;
		iterator end();
		const_iterator cend() const;
	protected:
		mutex m_mutex;
	private:
		shared_ptr<IParamFunc> m_function;
		shared_ptr<IOptimalityFunction> m_optimality;
		shared_ptr<IParamCheck> m_filter;
		vector<pair<ParamSet,double>> m_population;
		ParamSet m_avr;
		ParamSet m_disp;
		ParamSet m_max_dev;
		unsigned int m_itercount;
		unsigned int threads;
	};
	template<class AdderOfPairs,class F, class S>
	inline shared_ptr<AdderOfPairs> operator<<(shared_ptr<AdderOfPairs> adder, pair<F,S> value){
		adder->Add(value.first,value.second);
		return adder;
	}
	#define use_num_type double
	#define use_indexer_type ParamSet&
	#include "math_h/wrap_func_indexer.h"
	#undef use_num_type
	#undef use_indexer_type
}
#endif

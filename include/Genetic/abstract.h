// this file is distributed under 
// MIT license
#ifndef ____WrKDhKHP___
#define ____WrKDhKHP___
#include <mutex>
#include <vector>
#include <memory>
#include <utility>
#include <functional>
#include <random>
#include <math.h>
#include <math_h/sigma.h>
#include <math_h/chains.h>
#include "paramset.h"
namespace Genetic{
	typedef std::mt19937 RANDOM;
	class IInitialConditions{
	public:
		virtual ~IInitialConditions(){}
		virtual ParamSet Generate(RANDOM&)const =0;
	};
	class IParamCheck{
	public:
		virtual ~IParamCheck(){}
		virtual bool operator()(const ParamSet&P)const=0;
	};
	class IOptimalityFunction{
	public:
		virtual ~IOptimalityFunction(){}
		virtual double operator()(const ParamSet&P)const=0;
	};
	class AbstractGenetic{
	protected:
		AbstractGenetic();
		AbstractGenetic(const std::shared_ptr<IOptimalityFunction> optimality);
	public:
		virtual ~AbstractGenetic();
		std::shared_ptr<IOptimalityFunction> OptimalityCalculator()const;
		
		AbstractGenetic&SetFilter(const std::shared_ptr<IParamCheck> filter);
		AbstractGenetic&SetFilter(const std::function<bool(const ParamSet&)>);
		AbstractGenetic&RemoveFilter();
		AbstractGenetic&SetThreadCount(const size_t threads_count);
		const size_t ThreadCount()const;
		
		AbstractGenetic&Init(const size_t population_size,const std::shared_ptr<IInitialConditions> initial_conditions,RANDOM&random);
		void Iterate(RANDOM&random);
		
		const unsigned long long int&iteration_count()const;
		const size_t PopulationSize()const;
		const size_t ParamCount()const;
		const double&Optimality(const size_t point_index=0)const;
		const ParamSet&Parameters(const size_t point_index=0)const;
		const std::vector<MathTemplates::value<double>>&ParametersStatistics()const;
		
		const bool ConcentratedInOnePoint()const;
		const bool AbsoluteOptimalityExitCondition(const double accuracy)const;
		const bool RelativeOptimalityExitCondition(const double accuracy)const;
		const bool ParametersDispersionExitCondition(const ParamSet&max_disp)const;
		const bool RelativeParametersDispersionExitCondition(const ParamSet&max_disp)const;
		const bool ParametersDispersionExitCondition(const ParamSet&&max_disp)const;
		const bool RelativeParametersDispersionExitCondition(const ParamSet&&max_disp)const;
	protected:
		//contains empty implementation for templates from genetic.h could work correctly
		virtual void mutations(ParamSet&,RANDOM&)const;
		virtual void HandleIteration();
	private:
		std::mutex m_mutex;
		std::shared_ptr<IOptimalityFunction> m_optimality;
		std::shared_ptr<IParamCheck> m_filter;
		MathTemplates::SortedChain<std::pair<ParamSet,double>> m_population;
		std::vector<MathTemplates::value<double>> m_stat;
		unsigned long long int m_itercount;
		size_t threads;
	};
	std::ostream&operator<<(std::ostream&str,const AbstractGenetic&P);
	inline void Find(AbstractGenetic&fit,RANDOM&engine){
		while(!fit.ConcentratedInOnePoint())
			fit.Iterate(engine);
	}

	class Filter:public IParamCheck{
	public:
		Filter(const std::function<bool(const ParamSet&)> c);
		virtual ~Filter();
		virtual bool operator()(const ParamSet&P)const override;
	private:
		std::function<bool(const ParamSet&)> condition;
	};
	class OptimalityFunction:public IOptimalityFunction{
	public:
		OptimalityFunction(const std::function<double(const ParamSet&)> f);
		virtual ~OptimalityFunction();
		virtual double operator()(const ParamSet&P)const override;
	private:
		std::function<double(const ParamSet&)> func;
	};
}
#endif

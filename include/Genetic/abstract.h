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
#include "paramset.h"
namespace Genetic{
	using namespace std;
	typedef mt19937 RANDOM;
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
	class Filter:public IParamCheck{
	public:
		Filter(const function<bool(const ParamSet&)> c);
		virtual ~Filter();
		virtual bool operator()(const ParamSet&P)const override;
	private:
		function<bool(const ParamSet&)> condition;
	};
	class IOptimalityFunction{
	public:
		virtual ~IOptimalityFunction(){}
		virtual double operator()(const ParamSet&P)const=0;
	};
	class OptimalityFunction:public IOptimalityFunction{
	public:
		OptimalityFunction(const function<double(const ParamSet&)> f);
		virtual ~OptimalityFunction();
		virtual double operator()(const ParamSet&P)const override;
	private:
		function<double(const ParamSet&)> func;
	};
	
	class AbstractGenetic{
	protected:
		AbstractGenetic();
		AbstractGenetic(const shared_ptr<IOptimalityFunction> optimality);
	public:
		virtual ~AbstractGenetic();
		shared_ptr<IOptimalityFunction> OptimalityCalculator()const;
		
		AbstractGenetic&SetFilter(const shared_ptr<IParamCheck> filter);
		AbstractGenetic&SetFilter(const function<bool(const ParamSet&)>);
		AbstractGenetic&RemoveFilter();
		
		AbstractGenetic&SetThreadCount(const size_t threads_count);
		const size_t ThreadCount()const;
		
		AbstractGenetic&Init(const size_t population_size,const shared_ptr<IInitialConditions> initial_conditions,RANDOM&random);
	protected:
		virtual void mutations(ParamSet&,RANDOM&)const;
	public:
		void Iterate(RANDOM&random);
		
		const unsigned long int iteration_count()const;
		const size_t PopulationSize()const;
		const size_t ParamCount()const;
		const double Optimality(const size_t point_index=0)const;
		const ParamSet&Parameters(const size_t point_index=0)const;
		const double operator[](const size_t i)const;
		const ParamSet&ParamAverage()const;
		const ParamSet&ParamDispersion()const;
		const ParamSet&ParamMaxDeviation()const;
		
		const bool ConcentratedInOnePoint()const;
		const bool AbsoluteOptimalityExitCondition(const double accuracy)const;
		const bool RelativeOptimalityExitCondition(const double accuracy)const;
		const bool ParametersDispersionExitCondition(const ParamSet&max_disp)const;
		const bool RelativeParametersDispersionExitCondition(const ParamSet&max_disp)const;
		const bool ParametersDispersionExitCondition(const ParamSet&&max_disp)const;
		const bool RelativeParametersDispersionExitCondition(const ParamSet&&max_disp)const;
		
		typedef vector<double>::const_iterator const_iterator;
		const_iterator begin()const;
		const_iterator cbegin()const;
		const_iterator end() const;
		const_iterator cend() const;
	protected:
		mutex m_mutex;
	private:
		shared_ptr<IOptimalityFunction> m_optimality;
		shared_ptr<IParamCheck> m_filter;
		vector<pair<ParamSet,double>> m_population;
		ParamSet m_avr;
		ParamSet m_disp;
		ParamSet m_max_dev;
		unsigned long int m_itercount;
		size_t threads;
	};
	ostream&operator<<(ostream&str,const AbstractGenetic&P);
	inline void Find(AbstractGenetic&fit,RANDOM&engine){
		while(!fit.ConcentratedInOnePoint())
			fit.Iterate(engine);
	}
}
#endif

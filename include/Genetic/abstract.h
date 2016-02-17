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
		virtual ParamSet Generate(RANDOM&)=0;
	};
	class IParamCheck{
	public:
		virtual ~IParamCheck(){}
		virtual bool operator()(const ParamSet&P)const=0;
	};
	class Filter:public IParamCheck{
	public:
		Filter(function<bool(const ParamSet&)> c);
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
		OptimalityFunction(function<double(const ParamSet&)> f);
		virtual ~OptimalityFunction();
		virtual double operator()(const ParamSet&P)const override;
	private:
		function<double(const ParamSet&)> func;
	};
	
	class AbstractGenetic{
	protected:
		AbstractGenetic();
		AbstractGenetic(shared_ptr<IOptimalityFunction> optimality);
	public:
		virtual ~AbstractGenetic();
		shared_ptr<IOptimalityFunction> OptimalityCalculator()const;
		
		AbstractGenetic&SetFilter(shared_ptr<IParamCheck> filter);
		AbstractGenetic&SetFilter(function<bool(const ParamSet&)>);
		AbstractGenetic&RemoveFilter();
		
		AbstractGenetic&SetThreadCount(size_t threads_count);
		size_t ThreadCount()const;
		
		AbstractGenetic&Init(size_t population_size,shared_ptr<IInitialConditions> initial_conditions,RANDOM&random);
	protected:
		virtual void mutations(ParamSet&,RANDOM&);
	public:
		void Iterate(RANDOM&random);
		
		unsigned long int iteration_count()const;
		size_t PopulationSize()const;
		size_t ParamCount()const;
		double Optimality(size_t point_index=0)const;
		const ParamSet&Parameters(size_t point_index=0)const;
		double operator[](size_t i)const;
		const ParamSet&ParamAverage()const;
		const ParamSet&ParamDispersion()const;
		const ParamSet&ParamMaxDeviation()const;
		
		bool ConcentratedInOnePoint()const;
		bool AbsoluteOptimalityExitCondition(double accuracy)const;
		bool RelativeOptimalityExitCondition(double accuracy)const;
		bool ParametersDispersionExitCondition(const ParamSet&max_disp)const;
		bool RelativeParametersDispersionExitCondition(const ParamSet&max_disp)const;
		bool ParametersDispersionExitCondition(ParamSet&&max_disp)const;
		bool RelativeParametersDispersionExitCondition(ParamSet&&max_disp)const;
		
		typedef vector<double>::iterator iterator;
		typedef vector<double>::const_iterator const_iterator;
		iterator begin();
		const_iterator begin()const;
		const_iterator cbegin()const;
		iterator end();
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

#ifndef ____WrKDhKHP___
#define ____WrKDhKHP___
#include <mutex>
#include <vector>
#include <memory>
#include <utility>
#include <functional>
#include <math.h>
#include "paramset.h"
namespace Genetic{
	using namespace std;
	class IInitialConditions{
	public:
		virtual ~IInitialConditions(){}
		virtual ParamSet Generate()=0;
	};
	class IParamCheck{
	public:
		virtual ~IParamCheck(){}
		virtual bool operator()(ParamSet&P)=0;
	};
	class Filter:public IParamCheck{
	public:
		Filter(function<bool(ParamSet&)> c);
		virtual ~Filter();
		virtual bool operator()(ParamSet&P)override;
	private:
		function<bool(ParamSet&)> condition;
	};
	class IOptimalityFunction{
	public:
		virtual ~IOptimalityFunction(){}
		virtual double operator()(ParamSet&P)=0;
	};
	class OptimalityFunction:public IOptimalityFunction{
	public:
		OptimalityFunction(function<double(ParamSet&)> f);
		virtual ~OptimalityFunction();
		virtual double operator()(ParamSet& P)override;
	private:
		function<double(ParamSet&)> func;
	};
	
	class AbstractGenetic{
	protected:
		AbstractGenetic();
		AbstractGenetic(shared_ptr<IOptimalityFunction> optimality);
	public:
		virtual ~AbstractGenetic();
		
		shared_ptr<IOptimalityFunction> OptimalityCalculator();
		void SetFilter(shared_ptr<IParamCheck> filter);
		void SetFilter(function<bool(ParamSet&)>);
		void RemoveFilter();
		
		void SetThreadCount(unsigned int threads_count);
		unsigned int ThreadCount();
		void Init(int population_size,shared_ptr<IInitialConditions> initial_conditions);
	protected:
		virtual void mutations(ParamSet&);
	public:
		void Iterate();
		
		unsigned long int iteration_count();
		int PopulationSize();
		int ParamCount();
		double Optimality(int point_index=0);
		ParamSet& Parameters(int point_index=0);
		double operator[](int i);
		ParamSet ParamAverage();
		ParamSet ParamDispersion();
		ParamSet ParamMaxDeviation();
		
		bool ConcentratedInOnePoint();
		bool AbsoluteOptimalityExitCondition(double accuracy);
		bool RelativeOptimalityExitCondition(double accuracy);
		
		typedef vector<double>::iterator iterator;
		typedef vector<double>::const_iterator const_iterator;
		iterator begin();
		const_iterator cbegin()const;
		iterator end();
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
		unsigned int threads;
	};
}
#endif

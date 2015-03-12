#ifndef ____WrKDhKHP___
#define ____WrKDhKHP___
#include <memory>
#include <mutex>
#include <vector>
#include <utility>
namespace Fit{
	using namespace std;
	class ParamSet{
	public:
		ParamSet();
		ParamSet(const ParamSet &source);
		ParamSet(double x);
		ParamSet(double x,double y);
		ParamSet(double x,double y,double z);
		ParamSet(double x,double y,double z,double zz);
		ParamSet(double x,double y,double z,double zz, double zzz);
		ParamSet(double x,double y,double z,double zz, double zzz, double zzzz);
		~ParamSet();
		ParamSet &operator=(const ParamSet &source);
		double operator[](int i);
		ParamSet &operator<<(double val);
		ParamSet &operator<<(ParamSet val);
		int Count();
		void Set(int i,double v);
		typedef vector<double>::iterator iterator;
		typedef vector<double>::const_iterator const_iterator;
		iterator begin();
		const_iterator cbegin()const;
		iterator end();
		const_iterator cend() const;
	protected:
		mutex m_mutex;
	private:
		vector<double> m_values;
	};
	template<class indexer> 
	ParamSet parFrom(int n,indexer x){
		ParamSet res;
		for(int i=0; i<n;i++)
			res<<(x[i]);
		return res;
	}
	ParamSet parEq(unsigned int cnt,double val);
	inline ParamSet parZeros(unsigned int cnt){
		return parEq(cnt,0);
	}
	inline ParamSet parOnes(unsigned int cnt){
		return parEq(cnt,1);
	}
	class IInitialConditions{
	public:
		virtual ~IInitialConditions(){}
		virtual ParamSet Generate()=0;
	};
	template<class AdderOfPairs>
	inline shared_ptr<AdderOfPairs> operator<<(shared_ptr<AdderOfPairs> adder, pair<double,double> value){
		adder->Add(value.first,value.second);
		return adder;
	}
	class IParamCheck{
	public:
		virtual ~IParamCheck(){}
		virtual bool CorrectParams(ParamSet params)=0;
	};
	class IParamFunc:public IParamCheck{
	public:
		virtual ~IParamFunc(){}
		virtual double operator()(ParamSet X, ParamSet P)=0;
	};
	class IOptimalityFunction{
	public:
		virtual ~IOptimalityFunction(){}
		virtual double operator()(ParamSet params, IParamFunc &func)=0;
	};
	class _gen{
	protected:
		_gen(shared_ptr<IParamFunc> function, shared_ptr<IOptimalityFunction> optimality);
	public:
		virtual ~_gen();
		void SetFilter(shared_ptr<IParamCheck> filter);
		void RemoveFilter();
		void Init(int population_size,shared_ptr<IInitialConditions> initial_conditions);
		void Iterate();
		int PopulationSize();
		ParamSet Parameters(int point_index=0);
		double Optimality(int point_index=0);
		int ParamCount();
		double operator[](int i);
		double operator()(ParamSet &X);
		unsigned int iteration_count();
		ParamSet ParamAverage();
		ParamSet ParamDispersion();
		ParamSet ParamMaxDeviation();
		ParamSet ParamParabolicError(ParamSet delta);
		shared_ptr<IParamFunc> Function();
		shared_ptr<IOptimalityFunction> OptimalityCalculator();
	protected:
		virtual ParamSet born(ParamSet)=0;
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
	};
	class FitGen: public _gen{
	public:
		FitGen(shared_ptr<IParamFunc> function, shared_ptr<IOptimalityFunction> optimality);
		virtual ~FitGen();
		double Mutation();
		void SetMutation(double val);
	protected:
		virtual ParamSet born(ParamSet C)override;
	private:
		double F;
	};
	class FitGenWithCrossing: public FitGen{
	public:
		FitGenWithCrossing(shared_ptr<IParamFunc> function, shared_ptr<IOptimalityFunction> optimality);
		virtual ~FitGenWithCrossing();
		double CrossingProbability();
		void SetCrossingProbability(double val);
	protected:
		virtual ParamSet born(ParamSet C)override;
	private:
		double P;
	};
}
#endif

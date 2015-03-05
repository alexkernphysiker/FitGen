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
		void operator=(const ParamSet &source);
		double operator[](int i);
		ParamSet &operator<<(double val);
		ParamSet &operator<<(ParamSet val);
		int Count();
		void Set(int i,double v);
	protected:
		mutex m_mutex;
	private:
		vector<double> m_values;
	};
	template<class indexer> 
	ParamSet CreateParamSet(int n,indexer x){
		ParamSet res;
		for(int i=0; i<n;i++)
			res<<(x[i]);
		return res;
	}
	ParamSet parEq(unsigned int cnt,double val);
	inline ParamSet parZeros(unsigned int cnt){return parEq(cnt,0);}
	inline ParamSet parOnes(unsigned int cnt){return parEq(cnt,1);}
	
	class IInitialConditions{
	public:
		virtual ~IInitialConditions(){}
		virtual ParamSet Generate()=0;
	};
	class IParamCheck{
	public:
		virtual ~IParamCheck(){}
		virtual bool CorrectParams(ParamSet &params)=0;
	};
	class IParamFunc:public IParamCheck{
	public:
		virtual ~IParamFunc(){}
		virtual double operator()(ParamSet &X, ParamSet &P)=0;
	};
	class IOptimalityFunction{
	public:
		virtual ~IOptimalityFunction(){}
		virtual double operator()(ParamSet &params, IParamFunc &func)=0;
	};
	
	class _gen{
	protected:
		_gen(std::shared_ptr<IParamFunc> function, std::shared_ptr<IOptimalityFunction> optimality);
	public:
		virtual ~_gen();
		void SetFilter(std::shared_ptr<IParamCheck> filter);
		void Init(int population_size,std::shared_ptr<IInitialConditions> initial_conditions);
		void Iterate();
		int PopulationSize();
		ParamSet GetParameters(int point_index=0);
		double GetOptimality(int point_index=0);
		int ParamCount();
		double operator[](int i);
		double operator()(ParamSet &X);
		unsigned int iteration_count();
		ParamSet ParamAverage();
		ParamSet ParamDispersion();
		ParamSet ParamMaxDeviation();
		ParamSet ParamParabolicError(ParamSet delta);
		shared_ptr<IParamFunc> GetFunction();
		shared_ptr<IOptimalityFunction> GetOptimalityCalculator();
	protected:
		virtual ParamSet born(ParamSet&)=0;
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
	
	enum MutationType{mutDifferential=0,mutRatio=1,mutAbsolute=2};
	class FitGenVeg: public _gen{
	public:
		FitGenVeg(std::shared_ptr<IParamFunc> function, std::shared_ptr<IOptimalityFunction> optimality);
		virtual ~FitGenVeg();
		ParamSet Mutation(MutationType index);
		void SetMutation(MutationType index,ParamSet val);
	protected:
		virtual ParamSet born(ParamSet &C)override;
	private:
		ParamSet m_Mut_Differential;
		ParamSet m_Mut_Ratio;
		ParamSet m_Mut_Absolute;
	};
	class FitGen: public FitGenVeg{
	public:
		FitGen(std::shared_ptr<IParamFunc> function, std::shared_ptr<IOptimalityFunction> optimality);
		virtual ~FitGen();
	protected:
		virtual ParamSet born(ParamSet &C)override;
	};
}
#endif

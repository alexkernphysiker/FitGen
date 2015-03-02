#ifndef FitGen_H
#define FitGen_H
#include <memory>
#include <mutex>
#include <vector>
namespace Fit{
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
	int Count();
	void Set(int i,double v);
protected:
	std::mutex m_mutex;
private:
	std::vector<double> m_values;
};
template<class indexer> 
ParamSet CreateParamSet(int n,indexer x){
	ParamSet res;
	for(int i=0; i<n;i++)
		res<<(x[i]);
	return res;
}
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
	_gen(std::shared_ptr<IParamFunc> function, std::shared_ptr<IOptimalityFunction> S);
public:
	virtual ~_gen();
	void SetFilter(std::shared_ptr<IParamCheck> filter);
	void Init(int N,std::shared_ptr<IInitialConditions> initial_conditions);
	void Iterate(unsigned char threads=1);
	int N();
	unsigned int iteration_count();
	double S(int point_index=0);
	ParamSet GetParameters(int point_index=0);
	int count();
	double operator[](int i);
	double operator()(ParamSet &X);
	ParamSet ParamDispersion();
	ParamSet ParamParabolicError(ParamSet delta);
	std::shared_ptr<IParamFunc> GetFunction();
	std::shared_ptr<IOptimalityFunction> GetOptimalityCalculator();
protected:
	ParamSet &Point(int point_index);
	virtual ParamSet born(ParamSet&)=0;
	std::mutex m_mutex;
private:
	std::shared_ptr<IParamFunc> m_function;
	std::shared_ptr<IOptimalityFunction> m_S;
	std::shared_ptr<IParamCheck> m_filter;
	std::vector<ParamSet> m_data;
	std::vector<double> S_cache;
	ParamSet m_disp;
	unsigned int m_itercount;
};
enum MutationType{mutDifferential,mutRatio,mutAbsolute};
class FitGenVeg: public _gen{
public:
	FitGenVeg(std::shared_ptr<IParamFunc> function, std::shared_ptr<IOptimalityFunction> S);
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
	FitGen(std::shared_ptr<IParamFunc> function, std::shared_ptr<IOptimalityFunction> S);
	virtual ~FitGen();
protected:
	virtual ParamSet born(ParamSet &C)override;
};
}
#endif // FitGen_H

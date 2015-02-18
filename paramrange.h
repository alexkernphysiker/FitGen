#ifndef PARAMRANGE_H
#define PARAMRANGE_H
#include "fit_gen.h"
namespace Fit{
class FilterRange:public IParamCheck{
public:
	virtual ~FilterRange(){}
	int Count();
	double Min(int i);
	double Max(int i);
	FilterRange &Add(double min,double max);
	virtual bool CorrectParams(ParamSet &params)override=0;
protected:
	std::vector<double> m_min;
	std::vector<double> m_max;
};
class FilterRangeIn:public FilterRange{
public:
	FilterRangeIn(){}
	FilterRangeIn(double min,double max){Add(min,max);}
	virtual ~FilterRangeIn(){}
	virtual bool CorrectParams(ParamSet &params)override;
};
class FilterRangeOut:public FilterRange{
public:
	FilterRangeOut(){}
	FilterRangeOut(double min,double max){Add(min,max);}
	virtual ~FilterRangeOut(){}
	virtual bool CorrectParams(ParamSet &params)override;
};

class FilterMulti:public IParamCheck{
public:
	virtual ~FilterMulti();
	int Count();
	IParamCheck &Get(int i);
	IParamCheck &operator ()(int i){return Get(i);}
	FilterMulti &Add(std::shared_ptr<IParamCheck> val);
	virtual bool CorrectParams(ParamSet &params)override=0;
protected:
	std::vector<std::shared_ptr<IParamCheck>> m_data;
};
class FilterAnd:public FilterMulti{
public:
	FilterAnd(){}
	FilterAnd(std::shared_ptr<IParamCheck> val){Add(val);}
	virtual ~FilterAnd(){}
	virtual bool CorrectParams(ParamSet &params)override;
};
class FilterOr:public FilterMulti{
public:
	FilterOr(){}
	FilterOr(std::shared_ptr<IParamCheck> val){Add(val);}
	virtual ~FilterOr(){}
	virtual bool CorrectParams(ParamSet &params)override;
};

class GenerateUniform:public IGenerator{
public:
	GenerateUniform();
	virtual ~GenerateUniform();
	int Count();
	double Min(int i);
	double Max(int i);
	GenerateUniform &Add(double min,double max);
	virtual ParamSet Generate()override;
private:
	std::vector<double> m_min;
	std::vector<double> m_max;
};

class GenerateByGauss:public IGenerator{
public:
	GenerateByGauss();
	virtual ~GenerateByGauss();
	int Count();
	double Mean(int i);
	double Sigma(int i);
	GenerateByGauss &Add(double mean,double sig);
	virtual ParamSet Generate()override;
private:
	std::vector<double> m_mean;
	std::vector<double> m_sig;
};
}
#endif // PARAMRANGE_H

#ifndef ____LRCfOgCQ___
#define ____LRCfOgCQ___
#include "fit_gen.h"
namespace Fit{
	using namespace std;
	
	template<class Filter>
	inline shared_ptr<Filter> operator<<(shared_ptr<Filter> filter, std::shared_ptr<IParamCheck> value){
		filter->Add(value);
		return filter;
	}
	
	
	class FilterAbove:public IParamCheck{
	public:
		FilterAbove(ParamSet v);
		virtual ~FilterAbove(){}
		virtual bool CorrectParams(ParamSet params)override;
	private:
		ParamSet m_data;
	};
	class FilterBelow:public IParamCheck{
	public:
		FilterBelow(ParamSet v);
		virtual ~FilterBelow(){}
		virtual bool CorrectParams(ParamSet params)override;
	private:
		ParamSet m_data;
	};
	
	class AbstractFilterRange:public IParamCheck{
	public:
		virtual ~AbstractFilterRange(){}
		int Count();
		double Min(int i);
		double Max(int i);
		AbstractFilterRange &Add(double min,double max);
		virtual bool CorrectParams(ParamSet params)override=0;
	protected:
		vector<double> m_min;
		vector<double> m_max;
	};
	class FilterRangeIn:public AbstractFilterRange{
	public:
		FilterRangeIn(){}
		virtual ~FilterRangeIn(){}
		virtual bool CorrectParams(ParamSet params)override;
	};
	class FilterRangeOut:public AbstractFilterRange{
	public:
		FilterRangeOut(){}
		virtual ~FilterRangeOut(){}
		virtual bool CorrectParams(ParamSet params)override;
	};
	
	class AbstractFilterMulti:public IParamCheck{
	public:
		virtual ~AbstractFilterMulti();
		int Count();
		IParamCheck &Get(int i);
		AbstractFilterMulti &Add(std::shared_ptr<IParamCheck> val);
		virtual bool CorrectParams(ParamSet params)override=0;
	protected:
		vector<shared_ptr<IParamCheck>> m_data;
	};
	class FilterAnd:public AbstractFilterMulti{
	public:
		FilterAnd(){}
		virtual ~FilterAnd(){}
		virtual bool CorrectParams(ParamSet params)override;
	};
	class FilterOr:public AbstractFilterMulti{
	public:
		FilterOr(){}
		virtual ~FilterOr(){}
		virtual bool CorrectParams(ParamSet params)override;
	};
}
#endif

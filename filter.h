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
	
	class Above:public IParamCheck{
	public:
		Above(ParamSet v);
		virtual ~Above(){}
		virtual bool CorrectParams(ParamSet params)override;
	private:
		ParamSet m_data;
	};
	class Below:public IParamCheck{
	public:
		Below(ParamSet v);
		virtual ~Below(){}
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
	class RangeIn:public AbstractFilterRange{
	public:
		RangeIn(){}
		virtual ~RangeIn(){}
		virtual bool CorrectParams(ParamSet params)override;
	};
	class RangeOut:public AbstractFilterRange{
	public:
		RangeOut(){}
		virtual ~RangeOut(){}
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
	class And:public AbstractFilterMulti{
	public:
		And(){}
		virtual ~And(){}
		virtual bool CorrectParams(ParamSet params)override;
	};
	class Or:public AbstractFilterMulti{
	public:
		Or(){}
		virtual ~Or(){}
		virtual bool CorrectParams(ParamSet params)override;
	};
}
#endif

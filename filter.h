#ifndef ____LRCfOgCQ___
#define ____LRCfOgCQ___
#include "fit_gen.h"
namespace Fit{
	using namespace std;
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
	class FilterRange:public IParamCheck{
	public:
		virtual ~FilterRange(){}
		int Count();
		double Min(int i);
		double Max(int i);
		FilterRange &Add(double min,double max);
		virtual bool CorrectParams(ParamSet params)override=0;
	protected:
		vector<double> m_min;
		vector<double> m_max;
	};
	class FilterRangeIn:public FilterRange{
	public:
		FilterRangeIn(){}
		FilterRangeIn(double min,double max){
			Add(min,max);
		}
		virtual ~FilterRangeIn(){}
		virtual bool CorrectParams(ParamSet params)override;
	};
	class FilterRangeOut:public FilterRange{
	public:
		FilterRangeOut(){}
		FilterRangeOut(double min,double max){
			Add(min,max);
		}
		virtual ~FilterRangeOut(){}
		virtual bool CorrectParams(ParamSet params)override;
	};
	
	class FilterMulti:public IParamCheck{
	public:
		virtual ~FilterMulti();
		int Count();
		IParamCheck &Get(int i);
		IParamCheck &operator ()(int i){
			return Get(i);
		}
		FilterMulti &Add(std::shared_ptr<IParamCheck> val);
		virtual bool CorrectParams(ParamSet params)override=0;
	protected:
		vector<shared_ptr<IParamCheck>> m_data;
	};
	class FilterAnd:public FilterMulti{
	public:
		FilterAnd(){}
		FilterAnd(shared_ptr<IParamCheck> val){
			Add(val);
		}
		virtual ~FilterAnd(){}
		virtual bool CorrectParams(ParamSet params)override;
	};
	class FilterOr:public FilterMulti{
	public:
		FilterOr(){}
		FilterOr(shared_ptr<IParamCheck> val){
			Add(val);
		}
		virtual ~FilterOr(){}
		virtual bool CorrectParams(ParamSet params)override;
	};
}
#endif

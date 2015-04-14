#ifndef CNTMUXAHBUIYZXIG
#define CNTMUXAHBUIYZXIG
#include "abstract.h"
#include "paramfunc.h"
namespace Genetic{
	using namespace std;
	class Above:public IParamCheck{
	public:
		Above();
		Above(ParamSet v);
		virtual ~Above(){}
		virtual bool operator()(ParamSet&P)override;
		Above &operator<<(double value);
	private:
		ParamSet m_data;
	};
	class Below:public IParamCheck{
	public:
		Below();
		Below(ParamSet v);
		virtual ~Below(){}
		virtual bool operator()(ParamSet&P)override;
		Below &operator<<(double value);
	private:
		ParamSet m_data;
	};
	template<class Filter>
	inline shared_ptr<Filter> operator<<(shared_ptr<Filter> filter, double value){
		filter->operator<<(value);
		return filter;
	}
	
	class AbstractFilterMulti:public IParamCheck{
	public:
		virtual ~AbstractFilterMulti();
		int Count();
		IParamCheck &Get(int i);
		AbstractFilterMulti &Add(shared_ptr<IParamCheck> val);
		AbstractFilterMulti &Add(function<bool(ParamSet&)> condition);
	protected:
		vector<shared_ptr<IParamCheck>> m_data;
	};
	shared_ptr<AbstractFilterMulti> operator<<(shared_ptr<AbstractFilterMulti> filter,shared_ptr<IParamCheck> value);
	shared_ptr<AbstractFilterMulti> operator<<(shared_ptr<AbstractFilterMulti> filter,function<bool(ParamSet&)> condition);
	class And:public AbstractFilterMulti{
	public:
		And(){}
		virtual ~And(){}
		virtual bool operator()(ParamSet&P)override;
	};
	class Or:public AbstractFilterMulti{
	public:
		Or(){}
		virtual ~Or(){}
		virtual bool operator()(ParamSet&P)override;
	};
}
#endif

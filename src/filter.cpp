// this file is distributed under 
// MIT license
#include <math_h/error.h>
#include <Genetic/filter.h>
namespace Genetic{
	using namespace std;
	using namespace MathTemplates;
	Above::Above(){}
	Above::Above(const ParamSet&&v):m_data(v){}
	bool Above::operator()(const ParamSet&P)const{
		bool res=true;
		int index=0;
		for(int i=0,n=m_data.size();i<n;i++){
			if(isfinite(m_data[i]))
				res&=(P[index]>=m_data[i]);
			index++;
		}
		return res;
	}
	Above& Above::operator<<(const double value){
		m_data<<value;
		return *this;
	}
	
	Below::Below(){}
	Below::Below(const ParamSet&&v):m_data(v){}
	bool Below::operator()(const ParamSet&P)const{
		bool res=true;
		int index=0;
		for(int i=0,n=m_data.size();i<n;i++){
			if(isfinite(m_data[i]))
				res&=(P[index]<=m_data[i]);
			index++;
		}
		return res;
	}
	Below& Below::operator<<(const double value){
		m_data<<value;
		return *this;
	}
	
	const size_t AbstractFilterMulti::Count()const{return m_data.size();}
	const IParamCheck &AbstractFilterMulti::Get(size_t i)const{
		if(i<m_data.size())
			return *m_data[i];
		else
			throw Exception<AbstractFilterMulti>("Range check error when getting an element from filter set");
	}
	AbstractFilterMulti &AbstractFilterMulti::Add(const std::shared_ptr<IParamCheck> val){
		m_data.push_back(val);
		return *this;
	}
	AbstractFilterMulti& AbstractFilterMulti::Add(const function<bool(const ParamSet&)> condition){
		m_data.push_back(make_shared<Filter>(condition));
		return *this;
	}
	shared_ptr<AbstractFilterMulti> operator<<(shared_ptr<AbstractFilterMulti> filter, const shared_ptr<IParamCheck> value){
		filter->Add(value);
		return filter;
	}
	shared_ptr<AbstractFilterMulti> operator<<(shared_ptr<AbstractFilterMulti> filter,const function<bool(const ParamSet&)> condition){
		filter->Add(condition);
		return filter;
	}
	AbstractFilterMulti::~AbstractFilterMulti(){}
	bool And::operator()(const ParamSet&P)const{
		for(auto f: m_data)
			if(!f->operator()(P))
				return false;
		return true;
	}
	bool Or::operator()(const ParamSet&P)const{
		for(auto f: m_data)
			if(f->operator()(P))
				return true;
		return false;
	}
}

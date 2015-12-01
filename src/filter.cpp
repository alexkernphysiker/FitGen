// this file is distributed under 
// MIT license
#include <exception_math_h.h>
#include <filter.h>
namespace Genetic{
	using namespace std;
	Above::Above(){}
	Above::Above(ParamSet&&v):m_data(v){}
	bool Above::operator()(const ParamSet&P)const{
		bool res=true;
		int index=0;
		for(int i=0,n=m_data.Count();i<n;i++){
			if(isfinite(m_data[i]))
				res&=(P[index]>=m_data[i]);
			index++;
		}
		return res;
	}
	Above& Above::operator<<(double value){
		m_data<<value;
		return *this;
	}
	
	Below::Below(){}
	Below::Below(ParamSet&&v):m_data(v){}
	bool Below::operator()(const ParamSet&P)const{
		bool res=true;
		int index=0;
		for(int i=0,n=m_data.Count();i<n;i++){
			if(isfinite(m_data[i]))
				res&=(P[index]<=m_data[i]);
			index++;
		}
		return res;
	}
	Below& Below::operator<<(double value){
		m_data<<value;
		return *this;
	}
	
	size_t AbstractFilterMulti::Count()const{return m_data.size();}
	IParamCheck &AbstractFilterMulti::Get(size_t i)const{
		if(i<m_data.size())
			return *m_data[i];
		else
			throw math_h_error<AbstractFilterMulti>("Range check error when getting an element from filter set");
	}
	AbstractFilterMulti &AbstractFilterMulti::Add(std::shared_ptr<IParamCheck> val){
		m_data.push_back(val);
		return *this;
	}
	AbstractFilterMulti& AbstractFilterMulti::Add(function<bool(const ParamSet&)> condition){
		m_data.push_back(make_shared<Filter>(condition));
		return *this;
	}
	shared_ptr<AbstractFilterMulti> operator<<(shared_ptr<AbstractFilterMulti> filter, shared_ptr<IParamCheck> value){
		filter->Add(value);
		return filter;
	}
	shared_ptr<AbstractFilterMulti> operator<<(shared_ptr<AbstractFilterMulti> filter,function<bool(const ParamSet&)> condition){
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

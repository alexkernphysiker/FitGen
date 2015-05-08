#include "filter.h"
#include "genetic_exception.h"
namespace Genetic{
	using namespace std;
	Above::Above(){}
	Above::Above(ParamSet v):m_data(v){}
	bool Above::operator()(ParamSet&P){
		bool res=true;
		int index=0;
		for(auto value:m_data){
			if(isfinite(value))
				res&=(P[index]>=value);
			index++;
		}
		return res;
	}
	Above& Above::operator<<(double value){
		m_data<<value;
		return *this;
	}
	
	Below::Below(){}
	Below::Below(ParamSet v):m_data(v){}
	bool Below::operator()(ParamSet&P){
		bool res=true;
		int index=0;
		for(auto value:m_data){
			if(isfinite(value))
				res&=(P[index]<=value);
			index++;
		}
		return res;
	}
	Below& Below::operator<<(double value){
		m_data<<value;
		return *this;
	}
	
	int AbstractFilterMulti::Count(){
		return m_data.size();
	}
	IParamCheck &AbstractFilterMulti::Get(int i){
		if((i>=0)&(i<m_data.size()))
			return *m_data[i];
		else
			throw GeneticException("Range check error when getting an element from filter set");
	}
	AbstractFilterMulti &AbstractFilterMulti::Add(std::shared_ptr<IParamCheck> val){
		m_data.push_back(val);
		return *this;
	}
	AbstractFilterMulti& AbstractFilterMulti::Add(function<bool(ParamSet&)> condition){
		m_data.push_back(make_shared<Filter>(condition));
		return *this;
	}
	shared_ptr<AbstractFilterMulti> operator<<(shared_ptr<AbstractFilterMulti> filter, shared_ptr<IParamCheck> value){
		filter->Add(value);
		return filter;
	}
	shared_ptr<AbstractFilterMulti> operator<<(shared_ptr<AbstractFilterMulti> filter,function<bool(ParamSet&)> condition){
		filter->Add(condition);
		return filter;
	}
	AbstractFilterMulti::~AbstractFilterMulti(){}
	bool And::operator()(ParamSet&P){
		for(auto f: m_data)
			if(!f->operator()(P))
				return false;
		return true;
	}
	bool Or::operator()(ParamSet&P){
		for(auto f: m_data)
			if(f->operator()(P))
				return true;
		return false;
	}
}

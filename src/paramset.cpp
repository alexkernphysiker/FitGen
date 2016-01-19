// this file is distributed under 
// MIT license
#include <math_h/exception_math_h.h>
#include <Genetic/paramset.h>
namespace Genetic{
	using namespace std;
	typedef lock_guard<mutex> Lock;
	ParamSet::ParamSet(){}
	ParamSet::ParamSet(const initializer_list< double >& source){
		for(auto value:source)m_values.push_back(value);
	}
	ParamSet::ParamSet(initializer_list< double >&& source):ParamSet(source){}
	ParamSet::ParamSet(const vector< double >& source){
		for(auto value:source)m_values.push_back(value);
	}
	ParamSet::ParamSet(vector<double>&&source):ParamSet(source){}
	ParamSet::ParamSet(const ParamSet&source):ParamSet(source.m_values){}
	ParamSet::ParamSet(ParamSet&&source):ParamSet(source.m_values){}
	ParamSet::~ParamSet(){}

	size_t ParamSet::size()const{return m_values.size();}
	double ParamSet::operator[](size_t i)const{
		if(i>=m_values.size())
			throw math_h_error<ParamSet>("Range check error when accessing ParamSet's element");
		return m_values[i];
	}
	double& ParamSet::operator[](size_t i){
		if(i>=m_values.size())
			throw math_h_error<ParamSet>("Range check error when accessing ParamSet's element");
		return m_values[i];
	}

	ParamSet& ParamSet::operator<<(double p){
		Lock lock(m_mutex);
		m_values.push_back(p);
		return *this;
	}
	ParamSet& ParamSet::operator<<(const initializer_list< double >& source){
		Lock lock(m_mutex);
		for(double p:source)m_values.push_back(p);
		return *this;
	}
	ParamSet& ParamSet::operator<<(initializer_list<double>&&source){return operator<<(source);}
	ParamSet& ParamSet::operator<<(const vector<double>&V){
		Lock lock(m_mutex);
		for(double p:V)m_values.push_back(p);
		return *this;
	}
	ParamSet& ParamSet::operator<<(vector<double>&&V){return operator<<(V);}
	ParamSet& ParamSet::operator<<(const ParamSet&P){return operator<<(P.m_values);}
	ParamSet& ParamSet::operator<<(ParamSet&&P){return operator<<(P);}
	
	ParamSet& ParamSet::operator>>(double& p){
		Lock lock(m_mutex);
		if(m_values.size()==0)
			throw math_h_error<ParamSet>("Attempt to take a value from empty paramset");
		p=m_values[m_values.size()-1];
		m_values.pop_back();
		return *this;
	}
	
	ParamSet& ParamSet::operator=(const initializer_list< double >& source){
		Lock lock(m_mutex);
		m_values.clear();
		for(double p:source)m_values.push_back(p);
		return *this;
	}
	ParamSet& ParamSet::operator=(initializer_list<double>&&source){return operator=(source);}
	ParamSet& ParamSet::operator=(const vector<double>&V){
		Lock lock(m_mutex);
		m_values.clear();
		for(double p:V)m_values.push_back(p);
		return *this;
	}
	ParamSet& ParamSet::operator=(vector< double >&& V){return operator=(V);}
	ParamSet& ParamSet::operator=(const ParamSet& P){return operator=(P.m_values);}
	ParamSet& ParamSet::operator=(ParamSet&& P){return operator=(P);}

	ParamSet::iterator ParamSet::begin(){return m_values.begin();}
	ParamSet::const_iterator ParamSet::begin()const{return m_values.begin();}
	ParamSet::const_iterator ParamSet::cbegin()const{return m_values.cbegin();}
	ParamSet::iterator ParamSet::end(){return m_values.end();}
	ParamSet::const_iterator ParamSet::end() const{return m_values.end();}
	ParamSet::const_iterator ParamSet::cend() const{return m_values.cend();}
	
	ostream& operator<<(ostream& str, const ParamSet& P){
		for(double x:P)str<<x<<"\t";
		return str;
	}
	istream& operator>>(istream& str, ParamSet& P){
		double x;
		while(str>>x)P<<x;
		return str;
	}
	
	ParamSet parEq(size_t cnt,double val){
		ParamSet res;
		for(size_t i=0;i<cnt;i++)
			res<<val;
		return res;
	}
	
	BinningParam::BinningParam(size_t ind, const pair< double, double >& R, size_t c){
		index=ind;cnt=c;
		range=R;
		CheckCorrectness();
	}
	BinningParam::BinningParam(size_t ind, pair< double, double >&& R, size_t cnt):BinningParam(ind,R,cnt){}
	BinningParam::BinningParam(const BinningParam& source):BinningParam(source.index,source.range,source.cnt){}
	pair<double,double>& BinningParam::limits()const{return const_cast<pair<double,double>&>(range);}
	size_t BinningParam::param_index() const{return index;}
	size_t BinningParam::count() const{return cnt;}
	void BinningParam::CheckCorrectness() const{
		if(range.second<=range.first)
			throw math_h_error<BinningParam>("wrong binning ranges");
		if(0==cnt)
			throw math_h_error<BinningParam>("there cannot be zero bins");
	}
	double BinningParam::bin_width() const{
		if(0==count())throw math_h_error<BinningParam>("count==0");
		return (range.second-range.first)/double(count());
	}
	double BinningParam::bin_center(size_t i) const{
		if(i>=count())throw math_h_error<BinningParam>("bin range check error");
		return range.first+bin_width()*(double(i)+0.5);
	}
	bool BinningParam::FindBinIndex(const ParamSet& P, size_t& res) const{
		double x=P[index],delta=bin_width()/2.0;
		if(delta<=0)throw math_h_error<BinningParam>("delta<=0");
		for(size_t i=0,n=count();i<n;i++){
			double pos=bin_center(i);
			if((x>=(pos-delta))&&(x<(pos+delta))){
				res=i;
				return true;
			}
		}
		return false;
	}
}
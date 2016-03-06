// this file is distributed under 
// MIT license
#ifndef PUWJZCZDMORMZODA
#define PUWJZCZDMORMZODA
#include <list>
#include <vector>
#include <mutex>
#include <iostream>
namespace Genetic{
	using namespace std;
	class ParamSet{
	public:
		ParamSet();
		ParamSet(const initializer_list<double>&source);
		ParamSet(const initializer_list<double>&&source);
		ParamSet(const vector<double> &source);
		ParamSet(const vector<double>&&source);
		ParamSet(const ParamSet&source);
		ParamSet(const ParamSet&&source);
		~ParamSet();
		
		ParamSet&operator<<(const double p);
		ParamSet&operator<<(const initializer_list<double>&source);
		ParamSet&operator<<(const initializer_list<double>&&source);
		ParamSet&operator<<(const vector<double>&V);
		ParamSet&operator<<(const vector<double>&&V);
		ParamSet&operator<<(const ParamSet&P);
		ParamSet&operator<<(const ParamSet&&P);
		ParamSet&operator>>(double&p);
		
		ParamSet&operator=(const initializer_list<double>&source);
		ParamSet&operator=(const initializer_list<double>&&source);
		ParamSet&operator=(const vector<double>&V);
		ParamSet&operator=(const vector<double>&&V);
		ParamSet&operator=(const ParamSet&P);
		ParamSet&operator=(const ParamSet&&P);
		
		const size_t size()const;
		double operator[](const size_t i)const;
		double&operator()(const size_t i);

		typedef vector<double>::iterator iterator;
		typedef vector<double>::const_iterator const_iterator;
		iterator begin();
		const_iterator begin()const;
		const_iterator cbegin()const;
		iterator end();
		const_iterator end() const;
		const_iterator cend() const;
	protected:
		mutex m_mutex;
	private:
		vector<double> m_values;
	};
	template<class indexer> 
	ParamSet parFrom(const size_t n,const indexer x){
		ParamSet res;
		for(size_t i=0; i<n;i++)
			res<<(x[i]);
		return res;
	}
	ostream&operator<<(ostream&str,const ParamSet&P);
	istream&operator>>(istream&str,ParamSet&P);
	ParamSet parEq(const size_t cnt,const double val);
	inline ParamSet parZeros(const size_t cnt){return parEq(cnt,0);}
	inline ParamSet parOnes(const size_t cnt){return parEq(cnt,1);}
};

#endif
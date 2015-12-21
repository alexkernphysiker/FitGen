// this file is distributed under 
// MIT license
#ifndef PUWJZCZDMORMZODA
#define PUWJZCZDMORMZODA
#include <vector>
#include <mutex>
#include <iostream>
namespace Genetic{
	using namespace std;
	class ParamSet{
	public:
		ParamSet();
		ParamSet(const ParamSet &source);
		ParamSet(double x);
		ParamSet(double x,double y);
		ParamSet(double x,double y,double z);
		ParamSet(double x,double y,double z,double zz);
		ParamSet(double x,double y,double z,double zz, double zzz);
		ParamSet(double x,double y,double z,double zz, double zzz, double zzzz);
		~ParamSet();

		ParamSet &operator=(const ParamSet &source);
		double operator[](size_t i)const;
		double&operator[](size_t i);
		ParamSet &operator<<(double val);
		ParamSet &operator<<(ParamSet val);
		size_t Count()const;

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
	ParamSet parFrom(size_t n,indexer x){
		ParamSet res;
		for(size_t i=0; i<n;i++)
			res<<(x[i]);
		return res;
	}
	ostream&operator<<(ostream&str,const ParamSet&P);
	istream&operator>>(istream&str,ParamSet&P);
	ParamSet parEq(size_t cnt,double val);
	inline ParamSet parZeros(size_t cnt){return parEq(cnt,0);}
	inline ParamSet parOnes(size_t cnt){return parEq(cnt,1);}
};

#endif
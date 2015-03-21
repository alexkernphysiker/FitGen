#ifndef PUWJZCZDMORMZODA
#define PUWJZCZDMORMZODA
#include <vector>
#include <mutex>
namespace Fit{
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
		double operator[](int i);
		ParamSet &operator<<(double val);
		ParamSet &operator<<(ParamSet val);
		int Count();
		void Set(int i,double v);
		typedef vector<double>::iterator iterator;
		typedef vector<double>::const_iterator const_iterator;
		iterator begin();
		const_iterator cbegin()const;
		iterator end();
		const_iterator cend() const;
	protected:
		mutex m_mutex;
	private:
		vector<double> m_values;
	};
	template<class indexer> 
	ParamSet parFrom(int n,indexer x){
		ParamSet res;
		for(int i=0; i<n;i++)
			res<<(x[i]);
		return res;
	}
	ParamSet parEq(unsigned int cnt,double val);
	inline ParamSet parZeros(unsigned int cnt){
		return parEq(cnt,0);
	}
	inline ParamSet parOnes(unsigned int cnt){
		return parEq(cnt,1);
	}
};

#endif
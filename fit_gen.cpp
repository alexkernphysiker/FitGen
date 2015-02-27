#include <thread>
#include "fit_gen.h"
#include "fitexception.h"
#include "math_h/interpolate.h"
#include "math_h/randomfunc.h"
#include "math_h/sigma.h"
using namespace std;
typedef lock_guard<mutex> Lock;
namespace Fit{
class EmptyFilter:public IParamCheck{
public:
	EmptyFilter(){}
	virtual ~EmptyFilter(){}
	virtual bool CorrectParams(ParamSet &)override{return true;}
};
_gen::_gen(std::shared_ptr<IParamFunc> function, std::shared_ptr<IOptimalityFunction> S){
	m_function=function;
	m_S=S;
	m_itercount=0;
	m_filter=std::make_shared<EmptyFilter>();
}
void _gen::Init(int N, std::shared_ptr<IInitialConditions> initial_conditions){
	if(m_data.size()>0)throw new FitException("Fitting algorithm cannot be inited twice");
	if(N<=0)throw new FitException("Fitting algorithm got incorrect parameters");
	m_itercount=0;
	for(int i=0;i<N;i++){
		{
			Lock lock(m_mutex);
			bool oncemore=true;
			while(oncemore){
				ParamSet tmp=initial_conditions->Generate();
				if(m_function->CorrectParams(tmp)&&m_filter->CorrectParams(tmp)){
					m_data.push_back(tmp);
					S_cache.push_back((*m_S)(tmp,*m_function));
					oncemore=false;
				}
			}
		}
	}
}
_gen::~_gen(){}
void _gen::SetFilter(std::shared_ptr<IParamCheck> filter){
	m_filter=filter;
}
int _gen::N(){
	Lock lock(m_mutex);
	return m_data.size();
}
uint _gen::iteration_count(){
	Lock lock(m_mutex);
	return m_itercount;
}
double _gen::S(int point_index){
	if((point_index<0)|(point_index>=S_cache.size()))
		throw new FitException("Point index out of range or no results were calculated");
	Lock lock(m_mutex);
	return S_cache[point_index];
}
ParamSet &_gen::Point(int point_index){
	if((point_index<0)|(point_index>=m_data.size()))
		throw new FitException("Point index out of range");
	Lock lock(m_mutex);
	return m_data[point_index];
}
ParamSet _gen::GetParameters(int point_index){return Point(point_index);}
int _gen::count(){
	if(m_data.size()==0)
		throw new FitException("Attempt to obtain unexisting results");
	Lock lock(m_mutex);
	return m_data[0].Count();
}
double _gen::operator [](int i){
	if((i<0)|(i>=count()))
		throw new FitException("Parameter index out of range");
	Lock lock(m_mutex);
	return m_data[0][i];
}
ParamSet _gen::GetParametersDispersion(){return m_disp;}
double _gen::operator ()(ParamSet &X){
	if(m_data.size()==0)
		throw new FitException("Attempt to obtain unexisting results");
	Lock lock(m_mutex);
	return (*m_function)(X,m_data[0]);
}
void _gen::Iterate(unsigned char threads){
	if(0==threads)
		throw new FitException("Cannot run fitting on zero threads");
	if(m_data.size()==0)
		throw new FitException("Fitting algorithm cannot work with zero size of population");
	S_tmp_cache.clear();m_tmp_data.clear();int n=N();
	auto func=[this](int a, int b){
		for(int i=a; i<=b; i++){
			double old_s;ParamSet oldC;
			{Lock lock(m_mutex);
				oldC=m_data[i];
				old_s=S_cache[i];
			}
			ParamSet newC=born(oldC);
			double new_s=+INFINITY;
			if(m_function->CorrectParams(newC)&& m_filter->CorrectParams(newC))
				new_s=m_S->operator ()(newC,*m_function);
			{Lock lock(m_mutex);
				int ind = WhereToInsert(0,S_tmp_cache.size()-1,S_tmp_cache,old_s);
				S_tmp_cache.insert(S_tmp_cache.begin()+ind,old_s);m_tmp_data.insert(m_tmp_data.begin()+ind,oldC);
				ind = WhereToInsert(0,S_tmp_cache.size()-1,S_tmp_cache,new_s);
				S_tmp_cache.insert(S_tmp_cache.begin()+ind,new_s);m_tmp_data.insert(m_tmp_data.begin()+ind,newC);
			}
		}
	};
	list<thread> pool;
	int partsize=n/threads;{
		int tf=threads-1;
		for(int t=0; t<tf;t++){
			int index=t*partsize;
			pool.push_back(thread(func,index,index+partsize-1));
		}
	}
	func((threads-1)*partsize,n-1);
	for(auto& thr:pool)
		thr.join();
	{
		Lock locker(m_mutex);
		m_data.clear();S_cache.clear();
		int cnt=m_tmp_data[0].Count();
		Sigma<double> disp[cnt];
		for(int i=0; i<n;i++){
			m_data.push_back(m_tmp_data[i]);
			S_cache.push_back(S_tmp_cache[i]);
			for(int j=0;j<cnt;j++)
				disp[j].AddValue(m_tmp_data[i][j]);
		}
		m_disp=ParamSet();
		for(int j=0;j<cnt;j++)
			m_disp<<disp[j].getSigma();
		m_tmp_data.clear();S_tmp_cache.clear();
		m_itercount++;
	}
}
std::shared_ptr<IParamFunc> _gen::GetFunction(){return m_function;}
std::shared_ptr<IOptimalityFunction> _gen::GetOptimalityCalculator(){return m_S;}
FitGenVeg::FitGenVeg(std::shared_ptr<IParamFunc> function, std::shared_ptr<IOptimalityFunction> S):
	_gen(function,S){}
FitGenVeg::~FitGenVeg(){}
ParamSet FitGenVeg::Mutation(int index){
	if(index==0)
		return m_f;
	if(index==1)
		return m_f_2;
	if(index==2)
		return m_f_3;
	throw new FitException("Index out of range");
}
void FitGenVeg::SetMutation(ParamSet val, int index){
	if(index==0){
		Lock locker(m_mutex);
		m_f=ParamSet();
		for(int i=0;i<val.Count();i++)
			if((val[i]<0)|(val[i]>1))
				throw new FitException("Attempt to set invalid M_0: should be in the range [0:1]");
			else
				m_f<<val[i];
	}
	if(index==1){
		Lock locker(m_mutex);
		for(int i=0;i<val.Count();i++)
			if(val[i]<0)
				throw new FitException("Attempt to set invalid M_1: should be a positive value");
			else
				m_f_2<<val[i];
	}
	if(index==2){
		Lock locker(m_mutex);
		for(int i=0;i<val.Count();i++)
			if(val[i]<0)
				throw new FitException("Attempt to set invalid M_2: should be a positive value");
			else
				m_f_3<<val[i];
	}
}
ParamSet FitGenVeg::born(ParamSet &C){
	ParamSet res;
	int a=rand()%N();int b=rand()%N();
	ParamSet A=Point(a), B=Point(b);
	for(int j=0; j<count();j++)
		res<<(C[j]+m_f[j]*(A[j]-B[j])+RandomGauss(m_f_2[j]*C[j])+RandomGauss(m_f_3[j]));
	return res;
}
FitGen::FitGen(std::shared_ptr<IParamFunc> function, std::shared_ptr<IOptimalityFunction> S):
	FitGenVeg(function,S){}
FitGen::~FitGen(){}
ParamSet FitGen::born(ParamSet &C){
	ParamSet res;
	ParamSet veg=FitGenVeg::born(C);
	int gen=rand()%N();
	for(int j=0; j<count();j++){
		if(rand()%2==0)
			res<<Point(gen)[j];
		else
			res<<veg[j];
	}
	return res;
}
}

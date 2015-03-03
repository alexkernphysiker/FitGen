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
_gen::_gen(std::shared_ptr<IParamFunc> function, std::shared_ptr<IOptimalityFunction> optimality){
	m_function=function;
	m_optimality=optimality;
	m_itercount=0;
	m_filter=std::make_shared<EmptyFilter>();
}
void _gen::Init(int population_size, std::shared_ptr<IInitialConditions> initial_conditions){
	if(m_data.size()>0)throw new FitException("Fitting algorithm cannot be inited twice");
	if(population_size<=0)throw new FitException("Fitting algorithm got incorrect parameters");
	Lock lock(m_mutex);
	m_itercount=0;
	for(int i=0;i<population_size;i++){
		bool oncemore=true;
		while(oncemore){
			ParamSet tmp=initial_conditions->Generate();
			if(m_function->CorrectParams(tmp)&&
				m_filter->CorrectParams(tmp)
			){
				m_data.push_back(tmp);
				S_cache.push_back((*m_optimality)(tmp,*m_function));
				oncemore=false;
			}
		}
	}
}
_gen::~_gen(){}
void _gen::SetFilter(std::shared_ptr<IParamCheck> filter){
	m_filter=filter;
}
int _gen::PopulationSize(){
	Lock lock(m_mutex);
	return m_data.size();
}
uint _gen::iteration_count(){
	Lock lock(m_mutex);
	return m_itercount;
}
double _gen::GetOptimality(int point_index){
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
int _gen::ParamCount(){
	if(m_data.size()==0)
		throw new FitException("Attempt to obtain unexisting results");
	Lock lock(m_mutex);
	return m_data[0].Count();
}
double _gen::operator [](int i){
	if((i<0)|(i>=ParamCount()))
		throw new FitException("Parameter index out of range");
	Lock lock(m_mutex);
	return m_data[0][i];
}
ParamSet _gen::ParamDispersion(){return m_disp;}
ParamSet _gen::ParamParabolicError(ParamSet delta){
	if(m_data.size()==0)
		throw new FitException("Attempt to calculate parabolic error with no results");
	if(delta.Count()!=m_data[0].Count())
		throw new FitException("Error in parabolic error calculation: incorrect delta size");
	ParamSet res;
	for(int i=0;i<delta.Count();i++){
		if(delta[i]<=0)
			throw new FitException("Error in parabolic error calculation: delta cannot be zero or negative");
		double s=GetOptimality();
		ParamSet ab=m_data[0];ab.Set(i,ab[i]+delta[i]);
		ParamSet be=m_data[0];be.Set(i,be[i]-delta[i]);
		double sa=m_optimality->operator()(ab,*m_function);
		double sb=m_optimality->operator()(be,*m_function);
		double da=(sa-s)/delta[i];double db=(s-sb)/delta[i];
		double dd=(da-db)/delta[i];
		if(dd<0){
			res<<INFINITY;
		}else{
			res<<sqrt(2.0/dd);
		}
	}
	return res;
}

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
	int n=PopulationSize();
	vector<ParamSet> tmp_data;
	vector<double> tmp_S;
	auto func=[this,&tmp_data,&tmp_S](int a, int b){
		for(int i=a; i<=b; i++){
			ParamSet oldC=m_data[i];
			double old_s=S_cache[i];
			ParamSet newC=born(oldC);
			double new_s=+INFINITY;
			if(
				m_function->CorrectParams(newC)&& 
				m_filter->CorrectParams(newC)
			)
				new_s=m_optimality->operator ()(newC,*m_function);
			{Lock lock(m_mutex);
				int ind = WhereToInsert(0,tmp_S.size()-1,tmp_S,old_s);
				tmp_S.insert(tmp_S.begin()+ind,old_s);
				tmp_data.insert(tmp_data.begin()+ind,oldC);
				ind = WhereToInsert(0,tmp_S.size()-1,tmp_S,new_s);
				tmp_S.insert(tmp_S.begin()+ind,new_s);
				tmp_data.insert(tmp_data.begin()+ind,newC);
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
	{Lock locker(m_mutex);
		m_data.clear();
		S_cache.clear();
		int cnt=tmp_data[0].Count();
		Sigma<double> disp[cnt];
		for(int i=0; i<n;i++){
			m_data.push_back(tmp_data[i]);
			S_cache.push_back(tmp_S[i]);
			for(int j=0;j<cnt;j++)
				disp[j].AddValue(tmp_data[i][j]);
		}
		m_disp=ParamSet();
		for(int j=0;j<cnt;j++)
			m_disp<<disp[j].getSigma();
		m_itercount++;
	}
}
std::shared_ptr<IParamFunc> _gen::GetFunction(){return m_function;}
std::shared_ptr<IOptimalityFunction> _gen::GetOptimalityCalculator(){return m_optimality;}
FitGenVeg::FitGenVeg(std::shared_ptr<IParamFunc> function, std::shared_ptr<IOptimalityFunction> optimality):
	_gen(function,optimality){}
FitGenVeg::~FitGenVeg(){}
ParamSet FitGenVeg::Mutation(MutationType index){
	switch(index){
		case mutDifferential:
			return m_Mut_Differential;
		case mutRatio:
			return m_Mut_Ratio;
		case mutAbsolute:
			return m_Mut_Absolute;
		default:
			throw new FitException("Mutation: Index out of range");
	};
}
void FitGenVeg::SetMutation(MutationType index,ParamSet val){
	switch(index){
		case mutDifferential:{
			Lock locker(m_mutex);
			m_Mut_Differential=ParamSet();
			for(int i=0;i<val.Count();i++)
				if((val[i]<0)|(val[i]>1))
					throw new FitException("Attempt to set invalid M_0: should be in the range [0:1]");
				else
					m_Mut_Differential<<val[i];
		}break;
		case mutRatio:{
			Lock locker(m_mutex);
			m_Mut_Ratio=ParamSet();
			for(int i=0;i<val.Count();i++)
				if(val[i]<0)
					throw new FitException("Attempt to set invalid M_1: should be a positive value");
				else
					m_Mut_Ratio<<val[i];
		}break;
		case mutAbsolute:{
			Lock locker(m_mutex);
			m_Mut_Absolute=ParamSet();
			for(int i=0;i<val.Count();i++)
				if(val[i]<0)
					throw new FitException("Attempt to set invalid M_2: should be a positive value");
				else
					m_Mut_Absolute<<val[i];
		}break;
		default:
			throw new FitException("Mutation: Index out of range");
	};
}
ParamSet FitGenVeg::born(ParamSet &C){
	ParamSet res;
	int a=rand()%PopulationSize();int b=rand()%PopulationSize();
	ParamSet A=Point(a), B=Point(b);
	for(int j=0; j<ParamCount();j++)
		res<<(C[j]+m_Mut_Differential[j]*(A[j]-B[j])
			+RandomGauss(m_Mut_Ratio[j]*C[j])
			+RandomGauss(m_Mut_Absolute[j]));
	return res;
}
FitGen::FitGen(std::shared_ptr<IParamFunc> function, std::shared_ptr<IOptimalityFunction> optimality):
	FitGenVeg(function,optimality){}
FitGen::~FitGen(){}
ParamSet FitGen::born(ParamSet &C){
	ParamSet res;
	ParamSet veg=FitGenVeg::born(C);
	int gen=rand()%PopulationSize();
	for(int j=0; j<ParamCount();j++){
		if(rand()%2==0)
			res<<Point(gen)[j];
		else
			res<<veg[j];
	}
	return res;
}
}

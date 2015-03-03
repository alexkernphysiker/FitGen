#include <thread>
#include "fit_gen.h"
#include "fitexception.h"
typedef std::pair<Fit::ParamSet,double> Point;
bool operator>(Point a,Point b){return a.second>b.second;}
bool operator<(Point a,Point b){return a.second<b.second;}
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
	if(m_population.size()>0)throw new FitException("Fitting algorithm cannot be inited twice");
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
				m_population.push_back(make_pair(tmp,m_optimality->operator()(tmp,*m_function)));
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
	return m_population.size();
}
uint _gen::iteration_count(){
	Lock lock(m_mutex);
	return m_itercount;
}
double _gen::GetOptimality(int point_index){
	if((point_index<0)|(point_index>=m_population.size()))
		throw new FitException("Point index out of range or no results were calculated");
	Lock lock(m_mutex);
	return m_population[point_index].second;
}
ParamSet _gen::GetParameters(int point_index){
	if((point_index<0)|(point_index>=m_population.size()))
		throw new FitException("Point index out of range or no results were calculated");
	Lock lock(m_mutex);
	return m_population[point_index].first;
}
int _gen::ParamCount(){
	if(m_population.size()==0)
		throw new FitException("Attempt to obtain unexisting results");
	Lock lock(m_mutex);
	return m_population[0].first.Count();
}
double _gen::operator [](int i){
	if((i<0)|(i>=ParamCount()))
		throw new FitException("Parameter index out of range");
	Lock lock(m_mutex);
	return m_population[0].first[i];
}
ParamSet _gen::ParamDispersion(){return m_disp;}
ParamSet _gen::ParamParabolicError(ParamSet delta){
	if(m_population.size()==0)
		throw new FitException("Attempt to calculate parabolic error with no results");
	if(delta.Count()!=m_population[0].first.Count())
		throw new FitException("Error in parabolic error calculation: incorrect delta size");
	ParamSet res;
	for(int i=0;i<delta.Count();i++){
		if(delta[i]<=0)
			throw new FitException("Error in parabolic error calculation: delta cannot be zero or negative");
		double s=GetOptimality();
		ParamSet ab=GetParameters();
		ParamSet be=ab;
		ab.Set(i,ab[i]+delta[i]);
		be.Set(i,be[i]-delta[i]);
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
	ParamSet P=GetParameters();
	return m_function->operator()(X,P);
}
void _gen::Iterate(){
	int n=m_population.size();
	if(n==0)
		throw new FitException("Fitting algorithm cannot work with zero size of population");
	vector<pair<ParamSet,double>> tmp_population;
	for(auto point:m_population){
		bool once_more=true;
		while(once_more){
			ParamSet new_member=born(point.first);
			if(
				m_function->CorrectParams(new_member)&& 
				m_filter->CorrectParams(new_member)
			){
				once_more=false;
				auto new_point=make_pair(new_member,m_optimality->operator()(new_member,*m_function));
				{Lock lock(m_mutex);
					tmp_population.insert(tmp_population.begin()+WhereToInsert(0,tmp_population.size()-1,tmp_population,point),point);
					tmp_population.insert(tmp_population.begin()+WhereToInsert(0,tmp_population.size()-1,tmp_population,new_point),new_point);
				}
			}
		}
	}
	int par_cnt=ParamCount();
	{Lock locker(m_mutex);
		Sigma<double> disp[par_cnt];
		m_population.clear();
		for(int i=0; i<n;i++){
			m_population.push_back(tmp_population[i]);
			for(int j=0;j<par_cnt;j++)
				disp[j].AddValue(tmp_population[i].first[j]);
		}
		m_disp=ParamSet();
		for(int j=0;j<par_cnt;j++)
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
	ParamSet A=GetParameters(a), 
		B=GetParameters(b);
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
			res<<GetParameters(gen)[j];
		else
			res<<veg[j];
	}
	return res;
}
}

#include <thread>
#include "fit_gen.h"
#include "fitexception.h"
#include "math_h/interpolate.h"
#include "math_h/randomfunc.h"
#include "math_h/sigma.h"
using namespace std;
namespace Fit{
	class EmptyFilter:public IParamCheck{
	public:
		EmptyFilter(){}
		virtual ~EmptyFilter(){}
		virtual bool CorrectParams(ParamSet &)override{return true;}
	};

	typedef lock_guard<mutex> Lock;
	typedef pair<ParamSet,double> Point;
	bool operator>(Point a,Point b){return a.second>b.second;}
	bool operator<(Point a,Point b){return a.second<b.second;}
	template<class Create,class Condition>
	inline ParamSet CreateNew(Create create,Condition condition){
		while(true){
			ParamSet res=create();
			if(condition(res))return res;
		}
	}
	
	_gen::_gen(shared_ptr<IParamFunc> function, shared_ptr<IOptimalityFunction> optimality){
		m_function=function;
		m_optimality=optimality;
		m_itercount=0;
		m_filter=std::make_shared<EmptyFilter>();
	}
	_gen::~_gen(){}
	void _gen::Init(int population_size, shared_ptr<IInitialConditions> initial_conditions){
		if(m_population.size()>0)throw new FitException("Fitting algorithm cannot be inited twice");
		if(population_size<=0)throw new FitException("Fitting algorithm got incorrect parameters");
		Lock lock(m_mutex);
		m_itercount=0;
		for(int i=0;i<population_size;i++){
			ParamSet new_member=CreateNew(
				[initial_conditions](){return initial_conditions->Generate();},
				[this](ParamSet p){return m_function->CorrectParams(p) && m_filter->CorrectParams(p);}
			);
			m_population.push_back(make_pair(new_member,m_optimality->operator()(new_member,*m_function)));
		}
	}
	void _gen::Iterate(){
		int n=PopulationSize();
		int par_cnt=ParamCount();
		if(n==0)
			throw new FitException("Fitting algorithm cannot work with zero size of population");
		vector<Point> tmp_population;
		for(auto point:m_population){
			ParamSet new_member=CreateNew(
				[this,&point](){return born(point.first);},
				[this](ParamSet p){return m_function->CorrectParams(p) && m_filter->CorrectParams(p);}
			);
			auto new_point=make_pair(new_member,m_optimality->operator()(new_member,*m_function));
			{Lock lock(m_mutex);
				InsertSorted(point,tmp_population,std_size(tmp_population),std_insert(tmp_population,Point));
				InsertSorted(new_point,tmp_population,std_size(tmp_population),std_insert(tmp_population,Point));
			}
		}
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
	shared_ptr<IParamFunc> _gen::GetFunction(){return m_function;}
	shared_ptr<IOptimalityFunction> _gen::GetOptimalityCalculator(){return m_optimality;}
	void _gen::SetFilter(shared_ptr<IParamCheck> filter){
		Lock lock(m_mutex);
		m_filter=filter;
	}
	int _gen::PopulationSize(){
		Lock lock(m_mutex);
		return m_population.size();
	}
	int _gen::ParamCount(){
		if(m_population.size()==0)
			throw new FitException("Attempt to obtain unexisting results");
		Lock lock(m_mutex);
		return m_population[0].first.Count();
	}
	unsigned int _gen::iteration_count(){
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
	double _gen::operator ()(ParamSet &X){
		ParamSet P=GetParameters();
		return m_function->operator()(X,P);
	}
	double _gen::operator [](int i){
		if((i<0)|(i>=ParamCount()))
			throw new FitException("Parameter index out of range");
		Lock lock(m_mutex);
		return m_population[0].first[i];
	}
	ParamSet _gen::ParamDispersion(){return m_disp;}
	ParamSet _gen::ParamParabolicError(ParamSet delta){
		if(PopulationSize()==0)
			throw new FitException("Attempt to calculate parabolic error with no results");
		if(delta.Count()!=ParamCount())
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
			double da=(sa-s)/delta[i];
			double db=(s-sb)/delta[i];
			double dd=(da-db)/delta[i];
			if(dd<=0)
				res<<INFINITY;
			else
				res<<sqrt(2.0/dd);
		}
		return res;
	}
	
	FitGenVeg::FitGenVeg(shared_ptr<IParamFunc> function, shared_ptr<IOptimalityFunction> optimality):_gen(function,optimality){}
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
	
	FitGen::FitGen(shared_ptr<IParamFunc> function, shared_ptr<IOptimalityFunction> optimality):FitGenVeg(function,optimality){}
	FitGen::~FitGen(){}
	ParamSet FitGen::born(ParamSet &C){
		ParamSet res;
		ParamSet veg=FitGenVeg::born(C);
		auto gen=GetParameters(rand()%PopulationSize());
		for(int j=0; j<ParamCount();j++){
			if(rand()%2==0)
				res<<gen[j];
			else
				res<<veg[j];
		}
		return res;
	}
}

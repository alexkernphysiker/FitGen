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
		virtual bool CorrectParams(ParamSet)override{return true;}
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
			ParamSet new_param=CreateNew(
				[initial_conditions](){return initial_conditions->Generate();},
				[this](ParamSet p){return m_function->CorrectParams(p) && m_filter->CorrectParams(p);}
			);
			auto new_point=make_pair(new_param,m_optimality->operator()(new_param,*m_function));
			InsertSorted(new_point,m_population,field_size(m_population),field_insert(m_population,Point));
		}
	}
	void _gen::Iterate(){
		int n=PopulationSize();
		int par_cnt=ParamCount();
		if(n==0)
			throw new FitException("Fitting algorithm cannot work with zero size of population");
		vector<Point> tmp_population;
		for(auto point:m_population){
			ParamSet new_param=CreateNew(
				[this,&point](){return born(point.first);},
				[this](ParamSet p){return m_function->CorrectParams(p) && m_filter->CorrectParams(p);}
			);
			auto new_point=make_pair(new_param,m_optimality->operator()(new_param,*m_function));
			{Lock lock(m_mutex);
				InsertSorted(point,tmp_population,std_size(tmp_population),std_insert(tmp_population,Point));
				InsertSorted(new_point,tmp_population,std_size(tmp_population),std_insert(tmp_population,Point));
			}
		}
		{Lock locker(m_mutex);
			Sigma<double> disp[par_cnt];
			m_max_dev=parZeros(par_cnt);
			m_population.clear();
			for(int i=0; i<n;i++){
				m_population.push_back(tmp_population[i]);
				for(int j=0;j<par_cnt;j++){
					disp[j].AddValue(tmp_population[i].first[j]);
					double dev=abs(tmp_population[i].first[j]-tmp_population[0].first[j]);
					if(dev>m_max_dev[j])
						m_max_dev.Set(j,dev);
				}
			}
			m_avr=ParamSet();
			m_disp=ParamSet();
			for(int j=0;j<par_cnt;j++){
				m_avr<<disp[j].getAverage();
				m_disp<<disp[j].getSigma();
			}
			m_itercount++;
		}
	}
	shared_ptr<IParamFunc> _gen::Function(){return m_function;}
	shared_ptr<IOptimalityFunction> _gen::OptimalityCalculator(){return m_optimality;}
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
	double _gen::Optimality(int point_index){
		if((point_index<0)|(point_index>=m_population.size()))
			throw new FitException("Point index out of range or no results were calculated");
		Lock lock(m_mutex);
		return m_population[point_index].second;
	}
	ParamSet _gen::Parameters(int point_index){
		if((point_index<0)|(point_index>=m_population.size()))
			throw new FitException("Point index out of range or no results were calculated");
		Lock lock(m_mutex);
		return m_population[point_index].first;
	}
	double _gen::operator ()(ParamSet &X){
		ParamSet P=Parameters();
		return m_function->operator()(X,P);
	}
	double _gen::operator [](int i){
		if((i<0)|(i>=ParamCount()))
			throw new FitException("Parameter index out of range");
		Lock lock(m_mutex);
		return m_population[0].first[i];
	}
	ParamSet _gen::ParamAverage(){return m_avr;}
	ParamSet _gen::ParamDispersion(){return m_disp;}
	ParamSet _gen::ParamMaxDeviation(){return m_max_dev;}
	ParamSet _gen::ParamParabolicError(ParamSet delta){
		if(PopulationSize()==0)
			throw new FitException("Attempt to calculate parabolic error with no results");
		if(delta.Count()!=ParamCount())
			throw new FitException("Error in parabolic error calculation: incorrect delta size");
		ParamSet res;
		for(int i=0;i<delta.Count();i++){
			if(delta[i]<=0)
				throw new FitException("Error in parabolic error calculation: delta cannot be zero or negative");
			double s=Optimality();
			ParamSet ab=Parameters();
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
	
	FitGen::FitGen(shared_ptr<IParamFunc> function, shared_ptr<IOptimalityFunction> optimality):_gen(function,optimality){
		F=0.5;
	}
	FitGen::~FitGen(){}
	double FitGen::Mutation(){return F;}
	void FitGen::SetMutation(double val){
		if(val<=0)
			throw new FitException("Invalid mutation coefficient value");
		F=val;
	}
	ParamSet FitGen::born(ParamSet C){
		ParamSet res;
		int a=rand()%PopulationSize();int b=rand()%PopulationSize();
		ParamSet A=Parameters(a), 
		B=Parameters(b);
		for(int j=0; j<ParamCount();j++)
			res<<(C[j]+F*(A[j]-B[j]));
		return res;
	}
	
	FitGenWithCrossing::FitGenWithCrossing(shared_ptr<IParamFunc> function, shared_ptr<IOptimalityFunction> optimality):
		FitGen(function,optimality){
		P=0.5;
	}
	FitGenWithCrossing::~FitGenWithCrossing(){}
	double FitGenWithCrossing::CrossingProbability(){return P;}
	void FitGenWithCrossing::SetCrossingProbability(double val){
		if((val<0)||(val>1))
			throw new FitException("Invalid crossing probability value");
		P=val;
	}
	ParamSet FitGenWithCrossing::born(ParamSet parent){
		auto X=FitGen::born(parent);
		if(P>0){
			auto C=FitGen::born(Parameters(rand()%PopulationSize()));
			for(int i=0; i<ParamCount();i++)
				if(RandomUniformly(0.0,1.0)<=P)
					X.Set(i,C[i]);
		}
		return X;
	}
}

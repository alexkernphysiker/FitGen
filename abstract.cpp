#include <thread>
#include <math.h>
#include "abstract.h"
#include "genetic_exception.h"
#include "math_h/interpolate.h"
#include "math_h/sigma.h"
using namespace std;
namespace Genetic{
	Filter::Filter(function<bool(ParamSet&)> c){
		condition=c;
	}
	Filter::~Filter(){}
	bool Filter::operator()(ParamSet& P){
		return condition(P);
	}
	OptimalityFunction::OptimalityFunction(function<double(ParamSet&)> f){
		func=f;
	}
	OptimalityFunction::~OptimalityFunction(){}
	double OptimalityFunction::operator()(ParamSet& P){
		return func(P);
	}
	
	typedef lock_guard<mutex> Lock;
	typedef pair<ParamSet,double> Point;
	bool operator>(Point a,Point b){return a.second>b.second;}
	bool operator<(Point a,Point b){return a.second<b.second;}
	template<class Create,class Condition>
	inline ParamSet CreateNew(Create create,Condition condition){
		while(true){//ToDo: this may be a problem source
			ParamSet res=create();
			if(condition(res))return res;
		}
	}
	AbstractGenetic::AbstractGenetic(){}
	AbstractGenetic::AbstractGenetic(
		shared_ptr<IOptimalityFunction> optimality
	):AbstractGenetic(){
		threads=thread::hardware_concurrency();
		if(threads==0)
			threads=1;
		m_optimality=optimality;
		m_itercount=0;
		m_filter=make_shared<Filter>([](ParamSet&){return true;});
	}
	AbstractGenetic::~AbstractGenetic(){}
	shared_ptr<IOptimalityFunction> AbstractGenetic::OptimalityCalculator(){
		Lock lock(m_mutex);
		return m_optimality;
	}
	void AbstractGenetic::SetFilter(shared_ptr<IParamCheck> filter){
		Lock lock(m_mutex);
		m_filter=filter;
	}
	void AbstractGenetic::SetFilter(function<bool(ParamSet&)> f){
		Lock lock(m_mutex);
		m_filter=make_shared<Filter>(f);
	}
	void AbstractGenetic::RemoveFilter(){
		Lock lock(m_mutex);
		m_filter=make_shared<Filter>([](ParamSet&){return true;});
	}
	
	void AbstractGenetic::SetThreadCount(unsigned int threads_count){
		Lock lock(m_mutex);
		if(threads_count==0)
			throw GeneticException("Cannot run genetic algorithm with zero threads");
		threads=threads_count;
	}
	unsigned int AbstractGenetic::ThreadCount(){
		Lock lock(m_mutex);
		return threads;
	}
	void AbstractGenetic::Init(int population_size, shared_ptr<IInitialConditions> initial_conditions){
		if(m_population.size()>0)
			throw GeneticException("Fitting algorithm cannot be inited twice");
		if(population_size<=0)
			throw GeneticException("Fitting algorithm got incorrect parameters");
		auto add_to_population=[this,initial_conditions](int count){
			for(int i=0;i<count;i++){
				double s=INFINITY;
				ParamSet new_param=CreateNew(
					[initial_conditions](){return initial_conditions->Generate();},
					[this,&s](ParamSet&p){
						if(!(m_filter->operator()(p)))return false;
						s=m_optimality->operator()(p);
						return isfinite(s)!=0;
					}
				);
				auto new_point=make_pair(new_param,s);
				{Lock lock(m_mutex);
					InsertSorted(new_point,m_population,field_size(m_population),field_insert(m_population,Point));
				}
			}
		};
		{
			int piece_size=population_size/threads;
			int rest=population_size%threads;
			vector<shared_ptr<thread>> thread_vector;
			for(int i=1;i<threads;i++)
				thread_vector.push_back(make_shared<thread>(add_to_population,piece_size));
			add_to_population(piece_size+rest);
			for(auto thr:thread_vector)
				thr->join();
		}
		{Lock lock(m_mutex);
			m_itercount=0;
		}
	}
	void AbstractGenetic::mutations(ParamSet&){}
	void AbstractGenetic::Iterate(){
		int n=PopulationSize();
		int par_cnt=ParamCount();
		if(n==0)
			throw GeneticException("Fitting algorithm cannot work with zero size of population");
		vector<Point> tmp_population;
		auto process_elements=[this,&tmp_population](int from,int to){
			for(int i=from;i<=to;i++){
				Point point;
				{Lock lock(m_mutex);
					point=m_population[i];
				}
				double s=INFINITY;
				ParamSet new_param=CreateNew(
					[this,&point](){
						ParamSet p=point.first;
						mutations(p);
						return p;
					},
					[this,&s](ParamSet&p){
						if(!(m_filter->operator()(p)))return false;
						s=m_optimality->operator()(p);
						return isfinite(s)!=0;
					}
				);
				auto new_point=make_pair(new_param,s);
				{Lock lock(m_mutex);
					InsertSorted(point,tmp_population,std_size(tmp_population),std_insert(tmp_population,Point));
					InsertSorted(new_point,tmp_population,std_size(tmp_population),std_insert(tmp_population,Point));
				}
			}
		};
		{
			int piece_size=n/threads;
			vector<shared_ptr<thread>> thread_vector;
			for(int i=1;i<threads;i++)
				thread_vector.push_back(make_shared<thread>(process_elements,(i-1)*piece_size,(i*piece_size)-1));
			process_elements((threads-1)*piece_size,n-1);
			for(auto thr:thread_vector)
				thr->join();
		}
		{Lock locker(m_mutex);
			Sigma<double> disp[par_cnt];
			m_max_dev=parZeros(par_cnt);
			m_population.clear();
			for(int i=0; i<n;i++){
				m_population.push_back(tmp_population[i]);
				for(int j=0;j<par_cnt;j++){
					disp[j].AddValue(tmp_population[i].first[j]);
					double dev=tmp_population[i].first[j]-tmp_population[0].first[j];
					if(dev<0)dev=-dev;
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
	unsigned long int AbstractGenetic::iteration_count(){
		Lock lock(m_mutex);
		return m_itercount;
	}
	int AbstractGenetic::PopulationSize(){
		Lock lock(m_mutex);
		return m_population.size();
	}
	int AbstractGenetic::ParamCount(){
		if(m_population.size()==0)
			throw GeneticException("Attempt to obtain unexisting results");
		Lock lock(m_mutex);
		return m_population[0].first.Count();
	}
	double AbstractGenetic::Optimality(int point_index){
		if(m_population.size()==0)
			throw GeneticException("Attempt to obtain unexisting results");
		if((point_index<0)|(point_index>=m_population.size()))
			throw GeneticException("Point index out of range or no results were calculated");
		Lock lock(m_mutex);
		return m_population[point_index].second;
	}
	ParamSet& AbstractGenetic::Parameters(int point_index){
		if(m_population.size()==0)
			throw GeneticException("Attempt to obtain unexisting results");
		if((point_index<0)|(point_index>=m_population.size()))
			throw GeneticException("Point index out of range or no results were calculated");
		Lock lock(m_mutex);
		return m_population[point_index].first;
	}
	double AbstractGenetic::operator [](int i){
		if((i<0)|(i>=ParamCount()))
			throw GeneticException("Parameter index out of range");
		Lock lock(m_mutex);
		return m_population[0].first[i];
	}
	ParamSet AbstractGenetic::ParamAverage(){
		Lock lock(m_mutex);
		return m_avr;
	}
	ParamSet AbstractGenetic::ParamDispersion(){
		Lock lock(m_mutex);
		return m_disp;
	}
	ParamSet AbstractGenetic::ParamMaxDeviation(){
		Lock lock(m_mutex);
		return m_max_dev;
	}
	
	bool AbstractGenetic::ConcentratedInOnePoint(){
		if(m_population.size()==0)
			throw GeneticException("Attempt to obtain unexisting results");
		if(m_itercount==0)
			return false;
		if(Optimality(PopulationSize()-1)>Optimality())
			return false;
		bool res=true;
		for(auto v:m_max_dev)
			res&=(v==0);
		return res;
	}
	bool AbstractGenetic::AbsoluteOptimalityExitCondition(double accuracy){
		if(m_population.size()==0)
			throw GeneticException("Attempt to obtain unexisting results");
		if(accuracy<0)
			throw GeneticException("Wrong optimality exit condition.");
		if(m_itercount==0)
			return false;
		if(accuracy==0)
			return Optimality(PopulationSize()-1)==Optimality();
		return (Optimality(PopulationSize()-1)-Optimality())<=accuracy;
	}
	bool AbstractGenetic::RelativeOptimalityExitCondition(double accuracy){
		if(m_population.size()==0)
			throw GeneticException("Attempt to obtain unexisting results");
		if(accuracy<0)
			throw GeneticException("Wrong optimality exit condition.");
		if(m_itercount==0)
			return false;
		if(accuracy==0)
			return Optimality(PopulationSize()-1)==Optimality();
		return Optimality(PopulationSize()-1)<=(Optimality()*(1.0+accuracy));
	}
	
	AbstractGenetic::iterator AbstractGenetic::begin(){
		if(m_population.size()==0)
			throw GeneticException("Population size is zero. Attempt to get unexisting results");
		return m_population[0].first.begin();
	}
	AbstractGenetic::const_iterator AbstractGenetic::cbegin()const{
		if(m_population.size()==0)
			throw GeneticException("Population size is zero. Attempt to get unexisting results");
		return m_population[0].first.cbegin();
	}
	AbstractGenetic::iterator AbstractGenetic::end(){
		if(m_population.size()==0)
			throw GeneticException("Population size is zero. Attempt to get unexisting results");
		return m_population[0].first.end();
	}
	AbstractGenetic::const_iterator AbstractGenetic::cend() const{
		if(m_population.size()==0)
			throw GeneticException("Population size is zero. Attempt to get unexisting results");
		return m_population[0].first.cend();
	}
}

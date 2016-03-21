// this file is distributed under 
// MIT license
#include <thread>
#include <math.h>
#include <math_h/error.h>
#include <math_h/interpolate.h>
#include <math_h/sigma.h>
#include <Genetic/abstract.h>
namespace Genetic{
	using namespace std;
	using namespace MathTemplates;
	Filter::Filter(function<bool(const ParamSet&)> c){
		condition=c;
	}
	Filter::~Filter(){}
	bool Filter::operator()(const ParamSet&P)const{return condition(P);}
	OptimalityFunction::OptimalityFunction(const function<double(const ParamSet&)> f){func=f;}
	OptimalityFunction::~OptimalityFunction(){}
	double OptimalityFunction::operator()(const ParamSet&P)const{return func(P);}
	
	typedef lock_guard<mutex> Lock;
	typedef pair<ParamSet,double> Point;
	bool operator>(const Point&a,const Point&b){return a.second>b.second;}
	bool operator<(const Point&a,const Point&b){return a.second<b.second;}
	template<class Create,class Condition>
	inline ParamSet CreateNew(const Create create,const Condition condition){
		//ToDo: this may be a problem source
		while(true){
			ParamSet res=create();
			if(condition(res))return res;
		}
	}
	AbstractGenetic::AbstractGenetic(){
		threads=1;
		m_itercount=0;
		m_filter=make_shared<Filter>([](const ParamSet&){return true;});
	}
	AbstractGenetic::AbstractGenetic(
		const shared_ptr<IOptimalityFunction> optimality
	):AbstractGenetic(){
		m_optimality=optimality;
	}
	AbstractGenetic::~AbstractGenetic(){}
	shared_ptr<IOptimalityFunction> AbstractGenetic::OptimalityCalculator()const{return m_optimality;}
	AbstractGenetic&AbstractGenetic::SetFilter(const shared_ptr<IParamCheck> filter){
		Lock lock(m_mutex);
		m_filter=filter;
		return *this;
	}
	AbstractGenetic&AbstractGenetic::SetFilter(const function<bool(const ParamSet&)> f){
		Lock lock(m_mutex);
		m_filter=make_shared<Filter>(f);
		return *this;
	}
	AbstractGenetic&AbstractGenetic::RemoveFilter(){
		Lock lock(m_mutex);
		m_filter=make_shared<Filter>([](const ParamSet&){return true;});
		return *this;
	}
	
	AbstractGenetic&AbstractGenetic::SetThreadCount(const size_t threads_count){
		Lock lock(m_mutex);
		if(threads_count==0)
			throw Exception<AbstractGenetic>("Thread count cannot be zero");
		threads=threads_count;
		return *this;
	}
	const size_t AbstractGenetic::ThreadCount()const{return threads;}
	AbstractGenetic&AbstractGenetic::Init(const size_t population_size, const shared_ptr<IInitialConditions> initial_conditions,RANDOM&random){
		if(m_population.size()>0)
			throw Exception<AbstractGenetic>("Genetic algorithm cannot be inited twice");
		if(population_size<=0)
			throw Exception<AbstractGenetic>("Polulation size must be a positive number");
		if(ThreadCount()>population_size)
			SetThreadCount(population_size);
		auto add_to_population=[this,initial_conditions,&random](size_t count){
			for(size_t i=0;i<count;i++){
				double s=INFINITY;
				ParamSet new_param=CreateNew(
					[this,initial_conditions,&random](){
						Lock lock(m_mutex);
						return initial_conditions->Generate(random);
					},
					[this,&s](const ParamSet&p){
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
			size_t piece_size=population_size/threads;
			size_t rest=population_size%threads;
			vector<shared_ptr<thread>> thread_vector;
			for(size_t i=1;i<threads;i++)
				thread_vector.push_back(make_shared<thread>(add_to_population,piece_size));
			add_to_population(piece_size+rest);
			for(auto thr:thread_vector)
				thr->join();
		}
		{Lock lock(m_mutex);
			m_itercount=0;
		}
		return *this;
	}
	void AbstractGenetic::mutations(ParamSet&,RANDOM&)const{}
	void AbstractGenetic::Iterate(RANDOM&random){
		size_t n=PopulationSize();
		size_t par_cnt=ParamCount();
		if(n==0)
			throw Exception<AbstractGenetic>("Cannot perform the calculation when population size is zero");
		if(ThreadCount()>n)
			SetThreadCount(n);
		vector<Point> tmp_population;
		auto process_elements=[this,&tmp_population,&random](size_t from,size_t to){
			for(size_t i=from;i<=to;i++){
				Point point;
				{Lock lock(m_mutex);
					point=m_population[i];
				}
				double s=INFINITY;
				ParamSet new_param=CreateNew(
					[this,&point,&random](){
						ParamSet p=point.first;
						mutations(p,random);
						return p;
					},
					[this,&s](const ParamSet&p){
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
			size_t piece_size=n/threads;
			vector<shared_ptr<thread>> thread_vector;
			for(size_t i=1;i<threads;i++)
				thread_vector.push_back(make_shared<thread>(process_elements,(i-1)*piece_size,(i*piece_size)-1));
			process_elements((threads-1)*piece_size,n-1);
			for(auto thr:thread_vector)
				thr->join();
		}
		{Lock locker(m_mutex);
			Sigma<double> disp[par_cnt];
			m_max_dev=parZeros(par_cnt);
			m_population.clear();
			for(size_t i=0; i<n;i++){
				m_population.push_back(tmp_population[i]);
				for(size_t j=0;j<par_cnt;j++){
					disp[j]<<tmp_population[i].first[j];
					double dev=tmp_population[i].first[j]-tmp_population[0].first[j];
					if(dev<0)dev=-dev;
					if(dev>m_max_dev[j])
						m_max_dev(j)=j;
				}
			}
			m_avr=ParamSet();
			m_disp=ParamSet();
			for(size_t j=0;j<par_cnt;j++){
				m_avr<<disp[j]().val();
				m_disp<<disp[j]().delta();
			}
			m_itercount++;
		}
	}
	const unsigned long long int AbstractGenetic::iteration_count()const{
		return m_itercount;
	}
	const size_t AbstractGenetic::PopulationSize()const{
		return m_population.size();
	}
	const size_t AbstractGenetic::ParamCount()const{
		if(m_population.size()==0)
			throw Exception<AbstractGenetic>("Cannot obtain any parameters when population size is zero");
		return m_population[0].first.size();
	}
	const double AbstractGenetic::Optimality(const size_t point_index)const{
		if(m_population.size()==0)
			throw Exception<AbstractGenetic>("Cannot obtain any parameters when population size is zero");
		if(point_index>=m_population.size())
			throw Exception<AbstractGenetic>("Range check error when accessing an element in the population");
		return m_population[point_index].second;
	}
	const ParamSet&AbstractGenetic::Parameters(const size_t point_index)const{
		if(m_population.size()==0)
			throw Exception<AbstractGenetic>("Cannot obtain any parameters when population size is zero");
		if(point_index>=m_population.size())
			throw Exception<AbstractGenetic>("Range check error when accessing an element in the population");
		return m_population[point_index].first;
	}
	const double AbstractGenetic::operator [](const size_t i)const{
		if(i>=ParamCount())
			throw Exception<AbstractGenetic>("Parameter index out of range");
		return m_population[0].first[i];
	}
	const ParamSet&AbstractGenetic::ParamAverage()const{
		return m_avr;
	}
	const ParamSet&AbstractGenetic::ParamDispersion()const{
		return m_disp;
	}
	const ParamSet&AbstractGenetic::ParamMaxDeviation()const{
		return m_max_dev;
	}
	
	const bool AbstractGenetic::ConcentratedInOnePoint()const{
		if(m_population.size()==0)
			throw Exception<AbstractGenetic>("Cannot obtain any parameters when population size is zero");
		if(m_itercount==0)
			return false;
		if(Optimality(PopulationSize()-1)>Optimality())
			return false;
		bool res=true;
		for(size_t i=0,n=m_max_dev.size();i<n;i++)
			res&=(m_max_dev[i]==0);
		return res;
	}
	const bool AbstractGenetic::AbsoluteOptimalityExitCondition(const double accuracy)const{
		if(m_population.size()==0)
			throw Exception<AbstractGenetic>("Cannot obtain any parameters when population size is zero");
		if(accuracy<0)
			throw Exception<AbstractGenetic>("Cannot use negative accuracy parameter");
		if(m_itercount==0)
			return false;
		if(accuracy==0)
			return Optimality(PopulationSize()-1)==Optimality();
		return (Optimality(PopulationSize()-1)-Optimality())<=accuracy;
	}
	const bool AbstractGenetic::RelativeOptimalityExitCondition(const double accuracy)const{
		if(m_population.size()==0)
			throw Exception<AbstractGenetic>("Cannot obtain any parameters when population size is zero");
		if(accuracy<0)
			throw Exception<AbstractGenetic>("Cannot use negative accuracy parameter");
		if(m_itercount==0)
			return false;
		if(accuracy==0)
			return Optimality(PopulationSize()-1)==Optimality();
		return Optimality(PopulationSize()-1)<=(Optimality()*(1.0+accuracy));
	}
	const bool AbstractGenetic::ParametersDispersionExitCondition(const ParamSet&max_disp)const{
		if(m_population.size()==0)
			throw Exception<AbstractGenetic>("Cannot obtain any parameters when population size is zero");
		if(m_itercount==0)
			return false;
		if(max_disp.size()>m_disp.size())
			throw Exception<AbstractGenetic>("There number of parameters of GA is less than number of maximum dispersion exit conditions");
		for(size_t i=0;i<max_disp.size();i++){
			double m=max_disp[i];
			if(isfinite(m)){
				if(m<0)
					throw Exception<AbstractGenetic>("Dispersion value cannot be negative");
				if(m<m_disp[i])
					return false;
			}
		}
		return true;
	}
	const bool AbstractGenetic::ParametersDispersionExitCondition(const ParamSet&& max_disp) const{return ParametersDispersionExitCondition(max_disp);}
	const bool AbstractGenetic::RelativeParametersDispersionExitCondition(const ParamSet&max_disp)const{
		if(m_population.size()==0)
			throw Exception<AbstractGenetic>("Cannot obtain any parameters when population size is zero");
		if(m_itercount==0)
			return false;
		if(max_disp.size()>m_disp.size())
			throw Exception<AbstractGenetic>("There number of parameters of GA is less than number of maximum dispersion exit conditions");
		for(size_t i=0;i<max_disp.size();i++){
			double m=max_disp[i];
			if(isfinite(m)){
				if(m<0)
					throw Exception<AbstractGenetic>("Dispersion value cannot be negative");
				m*=m;
				if(m<pow(m_disp[i]/m_population[0].first[i],2))
					return false;
			}
		}
		return true;
	}
	const bool AbstractGenetic::RelativeParametersDispersionExitCondition(const ParamSet&& max_disp) const{return RelativeParametersDispersionExitCondition(max_disp);}
	
	AbstractGenetic::const_iterator AbstractGenetic::begin()const{
		if(m_population.size()==0)
			throw Exception<AbstractGenetic>("Cannot obtain any parameters when population size is zero");
		return m_population[0].first.begin();
	}
	AbstractGenetic::const_iterator AbstractGenetic::cbegin()const{
		if(m_population.size()==0)
			throw Exception<AbstractGenetic>("Cannot obtain any parameters when population size is zero");
		return m_population[0].first.cbegin();
	}
	AbstractGenetic::const_iterator AbstractGenetic::end() const{
		if(m_population.size()==0)
			throw Exception<AbstractGenetic>("Cannot obtain any parameters when population size is zero");
		return m_population[0].first.end();
	}
	AbstractGenetic::const_iterator AbstractGenetic::cend() const{
		if(m_population.size()==0)
			throw Exception<AbstractGenetic>("Cannot obtain any parameters when population size is zero");
		return m_population[0].first.cend();
	}
	ostream& operator<<(ostream& str, const AbstractGenetic& P){
		return str<<P.Parameters();
	}
}

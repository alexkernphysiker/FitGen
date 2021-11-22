// this file is distributed under
// LGPLv3 license
#ifdef using_multithread
#include <thread>
#endif
#include <math.h>
#include <math_h/error.h>
#include <math_h/interpolate.h>
#include <Genetic/abstract.h>
namespace Genetic
{
using namespace std;
using namespace MathTemplates;
#ifdef using_multithread
typedef lock_guard<mutex> Lock;
#endif

AbstractGenetic::AbstractGenetic()
{
#ifdef using_multithread
    threads = 1;
#endif
    m_itercount = 0;
    m_filter = make_shared<Filter>([](const ParamSet &) {
        return true;
    });
}
AbstractGenetic::AbstractGenetic(AbstractGenetic&&source)
    :m_population(std::move(source.m_population))
{
#ifdef using_multithread
    threads = source.threads;
#endif
    m_itercount = source.m_itercount;
    m_filter = source.m_filter;
    m_optimality=source.m_optimality;
    m_stat=std::move(source.m_stat);
    source.m_stat={};
}
AbstractGenetic::AbstractGenetic(
    const shared_ptr<IOptimalityFunction> optimality
): AbstractGenetic()
{
    m_optimality = optimality;
}
AbstractGenetic::~AbstractGenetic() {}
shared_ptr<IOptimalityFunction> AbstractGenetic::OptimalityCalculator()const
{
    return m_optimality;
}
AbstractGenetic &AbstractGenetic::SetFilter(const shared_ptr<IParamCheck> filter)
{
#ifdef using_multithread
    Lock lock(m_mutex);
#endif
    m_filter = filter;
    return *this;
}
AbstractGenetic &AbstractGenetic::SetFilter(const function<bool(const ParamSet &)> f)
{
#ifdef using_multithread
    Lock lock(m_mutex);
#endif
    m_filter = make_shared<Filter>(f);
    return *this;
}
AbstractGenetic &AbstractGenetic::RemoveFilter()
{
#ifdef using_multithread
    Lock lock(m_mutex);
#endif
    m_filter = make_shared<Filter>([](const ParamSet &) {
        return true;
    });
    return *this;
}
#ifdef using_multithread
AbstractGenetic &AbstractGenetic::SetThreadCount(const size_t threads_count)
{
    Lock lock(m_mutex);
    if (threads_count == 0)
        throw Exception<AbstractGenetic>("Thread count cannot be zero");
    threads = threads_count;
    return *this;
}
const size_t AbstractGenetic::ThreadCount()const
{
    return threads;
}
#endif
AbstractGenetic &AbstractGenetic::Init(const size_t population_size, const shared_ptr<IInitialConditions> initial_conditions)
{
    if (m_population.size() > 0)
        throw Exception<AbstractGenetic>("Genetic algorithm cannot be inited twice");
    if (population_size <= 0)
        throw Exception<AbstractGenetic>("Polulation size must be a positive number");
#ifdef using_multithread
    if (ThreadCount() > population_size)
        SetThreadCount(population_size);
#endif
    auto add_to_population = [this, initial_conditions](size_t count) {
        for (size_t i = 0; i < count; i++) {
            double s = INFINITY;
        ParamSet new_param;
        while (true) {
            new_param = initial_conditions->Generate();
            if (!(m_filter->operator()(new_param)))continue;
            s = m_optimality->operator()(new_param);
            if(isfinite(s))break;
        }
            auto new_point = make_point(s,new_param);
            {
#ifdef using_multithread
                Lock lock(m_mutex);
#endif
                m_population << new_point;
            }
        }
    };
    {
#ifdef using_multithread
        size_t piece_size = population_size / threads;
        size_t rest = population_size % threads;
        vector<shared_ptr<thread>> thread_vector;
        for (size_t i = 1; i < threads; i++)
            thread_vector.push_back(make_shared<thread>(add_to_population, piece_size));
        add_to_population(piece_size + rest);
        for (auto thr : thread_vector)
            thr->join();
#else
        add_to_population(population_size);
#endif
    }
    {
#ifdef using_multithread
        Lock lock(m_mutex);
#endif
        m_itercount = 0;
    }
    return *this;
}
void AbstractGenetic::Iterate()
{
    size_t n = PopulationSize();
    size_t par_cnt = ParamCount();
#ifdef using_multithread
    if (ThreadCount() > n)
        SetThreadCount(n);
#endif
    SortedPoints<double,ParamSet> tmp_population;
    auto process_elements = [this, &tmp_population](size_t from, size_t to) {
        for (size_t i = from; i <= to; i++) {
	    ParamSet new_param;
            double s = INFINITY;
	    while (true) {
		new_param = m_population[i].Y();
		mutations(new_param);
		if (!(m_filter->operator()(new_param)))continue;
		s = m_optimality->operator()(new_param);
		if(isfinite(s))break;
	    }
            auto new_point = make_point(s,new_param);
            {
#ifdef using_multithread
                Lock lock(m_mutex);
#endif
                tmp_population << m_population[i] << new_point;
            }
        }
    };
    {
#ifdef using_multithread
        size_t piece_size = n / threads;
        vector<shared_ptr<thread>> thread_vector;
        for (size_t i = 1; i < threads; i++)
            thread_vector.push_back(make_shared<thread>(process_elements, (i - 1)*piece_size, (i * piece_size) - 1));
        process_elements((threads - 1)*piece_size, n - 1);
        for (auto thr : thread_vector)
            thr->join();
#else
        process_elements(0,n - 1);
#endif
    }
    {
#ifdef using_multithread
        Lock locker(m_mutex);
#endif
        m_population = tmp_population.IndexRange(0, n);
        StandardDeviation<double> STAT[par_cnt];
        for (size_t i = 0; i < n; i++)
            for (size_t j = 0; j < par_cnt; j++)
                STAT[j] << tmp_population[i].Y()[j];
        m_stat.clear();
        for (size_t j = 0; j < par_cnt; j++)
            m_stat.push_back(STAT[j]);
        m_itercount++;
    }
    {
#ifdef using_multithread
        Lock locker(m_mutex);
#endif
        HandleIteration();
    }
}
const unsigned long long int &AbstractGenetic::iteration_count()const
{
    return m_itercount;
}
const size_t AbstractGenetic::PopulationSize()const
{
    return m_population.size();
}
const size_t AbstractGenetic::ParamCount()const
{
    if (m_population.size() == 0)
        throw Exception<AbstractGenetic>("Cannot obtain any parameters when population size is zero");
    return m_population[0].Y().size();
}
const double &AbstractGenetic::Optimality(const size_t point_index)const
{
    if (m_population.size() == 0)
        throw Exception<AbstractGenetic>("Cannot obtain any parameters when population size is zero");
    if (point_index >= m_population.size())
        throw Exception<AbstractGenetic>("Range check error when accessing an element in the population");
    return m_population[point_index].X();
}
const ParamSet &AbstractGenetic::Parameters(const size_t point_index)const
{
    if (m_population.size() == 0)
        throw Exception<AbstractGenetic>("Cannot obtain any parameters when population size is zero");
    if (point_index >= m_population.size())
        throw Exception<AbstractGenetic>("Range check error when accessing an element in the population");
    return m_population[point_index].Y();
}
const vector< value< double > > &AbstractGenetic::ParametersStatistics() const
{
    return m_stat;
}

const bool AbstractGenetic::ConcentratedInOnePoint()const
{
    if (m_population.size() == 0)
        throw Exception<AbstractGenetic>("Cannot obtain any parameters when population size is zero");
    if (m_itercount == 0)
        return false;
    if (Optimality(PopulationSize() - 1) > Optimality())
        return false;
    bool res = true;
    for (size_t i = 0, n = m_stat.size(); i < n; i++)
        res &= (m_stat[i].uncertainty() == 0);
    return res;
}
const bool AbstractGenetic::AbsoluteOptimalityExitCondition(const double &accuracy)const
{
    if (m_population.size() == 0)
        throw Exception<AbstractGenetic>("Cannot obtain any parameters when population size is zero");
    if (accuracy < 0)
        throw Exception<AbstractGenetic>("Cannot use negative accuracy parameter");
    if (m_itercount == 0)
        return false;
    if (accuracy == 0)
        return Optimality(PopulationSize() - 1) == Optimality();
    return (Optimality(PopulationSize() - 1) - Optimality()) <= accuracy;
}
const bool AbstractGenetic::RelativeOptimalityExitCondition(const double &accuracy) const
{
    if (m_population.size() == 0)
        throw Exception<AbstractGenetic>("Cannot obtain any parameters when population size is zero");
    if (accuracy < 0)
        throw Exception<AbstractGenetic>("Cannot use negative accuracy parameter");
    if (m_itercount == 0)
        return false;
    if (accuracy == 0)
        return Optimality(PopulationSize() - 1) == Optimality();
    return Optimality(PopulationSize() - 1) <= (Optimality() * (1.0 + accuracy));
}
const bool AbstractGenetic::ParametersDispersionExitCondition(const ParamSet &max_disp)const
{
    if (m_population.size() == 0)
        throw Exception<AbstractGenetic>("Cannot obtain any parameters when population size is zero");
    if (m_itercount == 0)
        return false;
    if (max_disp.size() > m_stat.size())
        throw Exception<AbstractGenetic>("There number of parameters of GA is less than number of maximum dispersion exit conditions");
    for (size_t i = 0; i < max_disp.size(); i++) {
        double m = max_disp[i];
        if (isfinite(m)) {
            if (m < 0)
                throw Exception<AbstractGenetic>("Dispersion value cannot be negative");
            if (m < m_stat[i].uncertainty())
                return false;
        }
    }
    return true;
}
const bool AbstractGenetic::RelativeParametersDispersionExitCondition(const ParamSet &max_disp)const
{
    if (m_population.size() == 0)
        throw Exception<AbstractGenetic>("Cannot obtain any parameters when population size is zero");
    if (m_itercount == 0)
        return false;
    if (max_disp.size() > m_stat.size())
        throw Exception<AbstractGenetic>("There number of parameters of GA is less than number of maximum dispersion exit conditions");
    for (size_t i = 0; i < max_disp.size(); i++) {
        double m = max_disp[i];
        if (isfinite(m)) {
            if (m < 0)
                throw Exception<AbstractGenetic>("Dispersion value cannot be negative");
            if (m < m_stat[i].epsilon())
                return false;
        }
    }
    return true;
}
void AbstractGenetic::mutations(ParamSet &)const {}
void AbstractGenetic::HandleIteration() {}

Filter::Filter(function<bool(const ParamSet &)> c)
{
    condition = c;
}
Filter::~Filter() {}
bool Filter::operator()(const ParamSet &P)const
{
    return condition(P);
}
}

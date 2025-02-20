// this file is distributed under
// LGPLv3 license
#include <math_h/error.h>
#include <Genetic/paramset.h>
namespace Genetic
{
    using namespace std;
    using namespace MathTemplates;
#ifdef using_multithread
    typedef lock_guard<mutex> Lock;
#endif
    ParamSet::ParamSet() {}
    ParamSet::ParamSet(const initializer_list< double >& source)
    {
        for (auto value : source)m_values.push_back(value);
    }
    ParamSet::ParamSet(const vector< double >& source)
    {
        for (auto value : source)m_values.push_back(value);
    }
    ParamSet::ParamSet(const ParamSet& source) : ParamSet(source.m_values) {}

    ParamSet::~ParamSet() {}

    const size_t ParamSet::size()const
    {
        return m_values.size();
    }
    const double& ParamSet::operator[](const size_t i)const
    {
        if (i >= m_values.size())
            throw Exception<ParamSet>("Range check error when accessing ParamSet's element");
        return m_values[i];
    }
    double& ParamSet::operator()(const size_t i)
    {
        if (i >= m_values.size())
            throw Exception<ParamSet>("Range check error when accessing ParamSet's element");
        return m_values[i];
    }

    ParamSet& ParamSet::push_back(const double& p)
    {
#ifdef using_multithread
        Lock lock(m_mutex);
#endif
        m_values.push_back(p);
        return *this;
    }

    ParamSet& ParamSet::operator=(const ParamSet& P)
    {
#ifdef using_multithread
        Lock lock(m_mutex);
#endif
        m_values.clear();
        for (double p : P)m_values.push_back(p);
        return *this;
    }

    ParamSet::iterator ParamSet::begin()
    {
        return m_values.begin();
    }
    ParamSet::const_iterator ParamSet::begin()const
    {
        return m_values.begin();
    }
    ParamSet::iterator ParamSet::end()
    {
        return m_values.end();
    }
    ParamSet::const_iterator ParamSet::end() const
    {
        return m_values.end();
    }

    ostream& operator<<(ostream& str, const ParamSet& P)
    {
        for (double x : P)str << x << "\t";
        return str;
    }
    istream& operator>>(istream& str, ParamSet& P)
    {
        double x;
        while (str >> x)P.push_back(x);
        return str;
    }

    ParamSet parEq(const size_t cnt, const double val)
    {
        ParamSet res;
        for (size_t i = 0; i < cnt; i++)
            res.push_back(val);
        return res;
    }
}

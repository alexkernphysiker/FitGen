// this file is distributed under
// LGPL license
#ifndef PUWJZCZDMORMZODA
#define PUWJZCZDMORMZODA
#include <list>
#include <vector>
#include <mutex>
#include <iostream>
namespace Genetic
{
class ParamSet
{
public:
    ParamSet();
    ParamSet(const std::initializer_list<double> &source);
    ParamSet(const std::vector<double> &source);
    ParamSet(const ParamSet &source);
    ~ParamSet();

    ParamSet &operator<<(const double p);
    ParamSet &operator<<(const std::vector<double> &V);
    ParamSet &operator<<(const ParamSet &P);

    ParamSet &operator=(const std::initializer_list<double> &source);
    ParamSet &operator=(const std::vector<double> &V);
    ParamSet &operator=(const ParamSet &P);

    const size_t size()const;
    const double &operator[](const size_t i)const;
    double &operator()(const size_t i);

    typedef std::vector<double>::iterator iterator;
    typedef std::vector<double>::const_iterator const_iterator;
    iterator begin();
    const_iterator begin()const;
    iterator end();
    const_iterator end() const;
protected:
#ifdef using_multithread
    std::mutex m_mutex;
#endif
private:
    std::vector<double> m_values;
};
template<class indexer>
ParamSet parFrom(const size_t n, const indexer x)
{
    ParamSet res;
    for (size_t i = 0; i < n; i++)
        res << (x[i]);
    return res;
}
std::ostream &operator<<(std::ostream &str, const ParamSet &P);
std::istream &operator>>(std::istream &str, ParamSet &P);
ParamSet parEq(const size_t cnt, const double val);
inline const ParamSet parZeros(const size_t cnt)
{
    return parEq(cnt, 0);
}
inline const ParamSet parOnes(const size_t cnt)
{
    return parEq(cnt, 1);
}
};

#endif

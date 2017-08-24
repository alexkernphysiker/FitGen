// this file is distributed under
// LGPL license
#ifndef CNTMUXAHBUIYZXIG
#define CNTMUXAHBUIYZXIG
#include "abstract.h"
#include "paramfunc.h"
namespace Genetic
{
class Above: public IParamCheck
{
public:
    Above();
    Above(const ParamSet &v);
    virtual ~Above() {}
    virtual bool operator()(const ParamSet &P)const override;
    Above &operator<<(const double &value);
private:
    ParamSet m_data;
};
class Below: public IParamCheck
{
public:
    Below();
    Below(const ParamSet &v);
    virtual ~Below() {}
    virtual bool operator()(const ParamSet &P)const override;
    Below &operator<<(const double &value);
private:
    ParamSet m_data;
};
template<class Filter>
inline std::shared_ptr<Filter> operator<<(std::shared_ptr<Filter> filter, const double value)
{
    filter->operator<<(value);
    return filter;
}

class AbstractFilterMulti: public IParamCheck
{
public:
    virtual ~AbstractFilterMulti();
    const size_t Count()const;
    const IParamCheck &Get(size_t i)const;
    AbstractFilterMulti &Add(const std::shared_ptr<IParamCheck> val);
    AbstractFilterMulti &Add(const std::function<bool(const ParamSet &)> condition);
protected:
    std::vector<std::shared_ptr<IParamCheck>> m_data;
};
std::shared_ptr<AbstractFilterMulti> operator<<(std::shared_ptr<AbstractFilterMulti>filter, const std::shared_ptr<IParamCheck> value);
std::shared_ptr<AbstractFilterMulti> operator<<(std::shared_ptr<AbstractFilterMulti>filter, const std::function<bool(const ParamSet &)> condition);
class And: public AbstractFilterMulti
{
public:
    And() {}
    virtual ~And() {}
    virtual bool operator()(const ParamSet &P)const override;
};
class Or: public AbstractFilterMulti
{
public:
    Or() {}
    virtual ~Or() {}
    virtual bool operator()(const ParamSet &P)const override;
};
}
#endif

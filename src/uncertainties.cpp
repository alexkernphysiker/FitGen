// this file is distributed under
// LGPLv3 license
#include <math.h>
#include <climits>
#include <math_h/error.h>
#include <Genetic/uncertainties.h>
namespace Genetic
{
using namespace std;
using namespace MathTemplates;
UncertaintiesEstimation::UncertaintiesEstimation()
{
    m_uncertainty_cache = make_shared<vector<value_numeric_distr<>>>();
}
UncertaintiesEstimation::UncertaintiesEstimation(UncertaintiesEstimation&&source)
    :AbstractGenetic(std::move(source))
{
    m_uncertainty_cache = source.m_uncertainty_cache;
    m_delta=source.m_delta;
}
UncertaintiesEstimation::~UncertaintiesEstimation() {}
double UncertaintiesEstimation::GetParamParabolicError(const double &delta, const size_t i)const
{
    if (delta <= 0)
        throw Exception<UncertaintiesEstimation>("Exception in parabolic error calculation: delta cannot be zero or negative");

    double s = Optimality();
    const auto& P = Parameters();
    const auto  S = [this](const ParamSet&p) { return OptimalityCalculator()->operator()(p); };
    ParamSet Pminus = P;
    ParamSet Pplus = P;
    Pplus(i) += delta;
    Pminus(i) -= delta;
    double second_deriv = ( S(Pminus) - 2.0 * S(P) + S(Pplus) ) / pow(delta, 2);
    if (second_deriv <= 0) {
        return INFINITY;
    } else {
        return sqrt( 2.0 / second_deriv );
    }
}
void UncertaintiesEstimation::HandleIteration()
{
    Genetic::AbstractGenetic::HandleIteration();
    m_uncertainty_cache->clear();
}
UncertaintiesEstimation &UncertaintiesEstimation::SetUncertaintyCalcDeltas(const ParamSet &P)
{
    m_delta = P;
    HandleIteration();
    return *this;
}

const vector<value_numeric_distr<double>> &UncertaintiesEstimation::ParametersWithUncertainties()const
{
    if (m_uncertainty_cache->size() == 0) {
        for (size_t i = 0; (i < ParamCount()) && (i < m_delta.size()); i++) {
            m_uncertainty_cache->push_back(value_numeric_distr<>(Parameters()[i], GetParamParabolicError(m_delta[i], i)));
        }
    }
    return *m_uncertainty_cache;
}

FunctionUncertaintiesEstimation::FunctionUncertaintiesEstimation():FunctionContainer(nullptr),UncertaintiesEstimation(){}
FunctionUncertaintiesEstimation::FunctionUncertaintiesEstimation(FunctionUncertaintiesEstimation&&source):FunctionContainer(source),UncertaintiesEstimation(std::move(source)){}
FunctionUncertaintiesEstimation::~FunctionUncertaintiesEstimation(){}
value<> FunctionUncertaintiesEstimation::FuncWithUncertainties(const ParamSet&X)const{
    const double val=FunctionContainer::func(X,AbstractGenetic::Parameters());
    double unc_sqr=0;
    const auto&PU=UncertaintiesEstimation::ParametersWithUncertainties();
    for(size_t i=0;i<AbstractGenetic::ParamCount();i++){
        ParamSet P=AbstractGenetic::Parameters();
        P(i)=PU[i].max();
        const double valu=FunctionContainer::func(X,P);
        P(i)=PU[i].min();
        const double vald=FunctionContainer::func(X,P);
        unc_sqr+=pow((valu-vald)/2.0,2);
    }
    return value<>(val,sqrt(unc_sqr));
}
}

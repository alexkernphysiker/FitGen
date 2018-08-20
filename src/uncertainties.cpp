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
UncertaintiesEstimation::UncertaintiesEstimation(const UncertaintiesEstimation&source)
    :AbstractGenetic(source)
{
    m_delta=source.m_delta;
    m_uncertainty_cache=source.m_uncertainty_cache;
}
UncertaintiesEstimation::~UncertaintiesEstimation() {}
double UncertaintiesEstimation::GetParamParabolicError(const double &delta, const size_t i)const
{
    if (delta <= 0)
        throw Exception<UncertaintiesEstimation>("Exception in parabolic error calculation: delta cannot be zero or negative");
    double s = Optimality();
    ParamSet ab = Parameters();
    ParamSet be = ab;
    ab(i) += delta;
    be(i) -= delta;
    double sa = OptimalityCalculator()->operator()(ab);
    double sb = OptimalityCalculator()->operator()(be);
    double dd = (sa - 2.0 * s + sb) / pow(delta, 2);
    if (dd <= 0)
        return INFINITY;
    else
        return sqrt(2.0 / dd);
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

FunctionUncertaintiesEstimation::FunctionUncertaintiesEstimation():UncertaintiesEstimation(),FunctionContainer(nullptr){}
FunctionUncertaintiesEstimation::FunctionUncertaintiesEstimation(const FunctionUncertaintiesEstimation&source):UncertaintiesEstimation(source),FunctionContainer(source){}
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

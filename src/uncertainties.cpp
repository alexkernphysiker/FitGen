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
ParabolicErrorEstimationFromChisq::ParabolicErrorEstimationFromChisq()
{
    m_uncertainty_cache = make_shared<vector<value_numeric_distr<>>>();
}
ParabolicErrorEstimationFromChisq::~ParabolicErrorEstimationFromChisq() {}
double ParabolicErrorEstimationFromChisq::GetParamParabolicError(const double &delta, const size_t i)const
{
    if (delta <= 0)
        throw Exception<ParabolicErrorEstimationFromChisq>("Exception in parabolic error calculation: delta cannot be zero or negative");
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
void ParabolicErrorEstimationFromChisq::HandleIteration()
{
    Genetic::AbstractGenetic::HandleIteration();
    m_uncertainty_cache->clear();
}
ParabolicErrorEstimationFromChisq &ParabolicErrorEstimationFromChisq::SetUncertaintyCalcDeltas(const ParamSet &P)
{
    m_delta = P;
    HandleIteration();
    return *this;
}

const vector<value_numeric_distr<double>> &ParabolicErrorEstimationFromChisq::ParametersWithUncertainties()const
{
    if (m_uncertainty_cache->size() == 0) {
        for (size_t i = 0; (i < ParamCount()) && (i < m_delta.size()); i++) {
            m_uncertainty_cache->push_back(value_numeric_distr<>(Parameters()[i], GetParamParabolicError(m_delta[i], i)));
        }
    }
    return *m_uncertainty_cache;
}

}
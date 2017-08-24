// this file is distributed under
// LGPL license
#ifndef ____PARABOLIC_H______
#	define ____PARABOLIC_H______
#include "abstract.h"
#include "genetic.h"
namespace Genetic
{
class ParabolicErrorEstimationFromChisq: public virtual AbstractGenetic
{
protected:
    ParabolicErrorEstimationFromChisq();
public:
    virtual ~ParabolicErrorEstimationFromChisq();
    ParabolicErrorEstimationFromChisq &SetUncertaintyCalcDeltas(const ParamSet &P);
    const std::vector<MathTemplates::value<double>> &ParametersWithUncertainties()const;
protected:
    virtual void HandleIteration()override;
private:
    double GetParamParabolicError(const double &delta, const size_t i)const;
    ParamSet m_delta;
    std::shared_ptr<std::vector<MathTemplates::value<double>>> m_uncertainty_cache;
};
typedef ParabolicErrorEstimationFromChisq Uncertainty;
}
#endif

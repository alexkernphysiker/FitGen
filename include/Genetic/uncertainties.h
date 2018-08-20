// this file is distributed under
// LGPLv3 license
#ifndef ____PARABOLIC_H______
#	define ____PARABOLIC_H______
#if __cplusplus<201100L
#error c++>=11 is needed for using math_h headers
#endif
#include <math_h/sigma2.h>
#include "abstract.h"
#include "genetic.h"
namespace Genetic
{
class UncertaintiesEstimation: public virtual AbstractGenetic
{
protected:
    UncertaintiesEstimation();
    UncertaintiesEstimation(const UncertaintiesEstimation&source);
public:
    virtual ~UncertaintiesEstimation();
    UncertaintiesEstimation &SetUncertaintyCalcDeltas(const ParamSet &P);
    const std::vector<MathTemplates::value_numeric_distr<>> &ParametersWithUncertainties()const;
protected:
    virtual void HandleIteration()override;
private:
    double GetParamParabolicError(const double &delta, const size_t i)const;
    ParamSet m_delta;
    std::shared_ptr<std::vector<MathTemplates::value_numeric_distr<>>> m_uncertainty_cache;
};
}
#endif

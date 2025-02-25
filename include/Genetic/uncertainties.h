// this file is distributed under
// LGPLv3 license
#ifndef ____PARABOLIC_H______
#	define ____PARABOLIC_H______
#if __cplusplus<201100L
#error c++>=11 is needed for using math_h headers
#endif
#include <math_h/sigma2.h>
#include "abstract.h"
#include "fit.h"
namespace Genetic
{
    class UncertaintiesEstimation : public virtual AbstractGenetic
    {
    protected:
        UncertaintiesEstimation();
        UncertaintiesEstimation(UncertaintiesEstimation&& source);
    public:
        virtual ~UncertaintiesEstimation();
        UncertaintiesEstimation& SetUncertaintyCalcDeltas(const ParamSet& P);
        const std::vector<MathTemplates::value_numeric_distr<>>& ParametersWithUncertainties()const;
    protected:
        virtual void HandleIteration()override;
    private:
        double GetParamParabolicError(const double& delta, const size_t i)const;
        ParamSet m_delta;
        std::shared_ptr<std::vector<MathTemplates::value_numeric_distr<>>> m_uncertainty_cache;
    };
    class FunctionUncertaintiesEstimation :public UncertaintiesEstimation, public virtual FunctionContainer {
    protected:
        FunctionUncertaintiesEstimation();
        FunctionUncertaintiesEstimation(FunctionUncertaintiesEstimation&& source);
    public:
        virtual ~FunctionUncertaintiesEstimation();
        MathTemplates::value<> FuncWithUncertainties(const ParamSet& X)const;
    };
}
#endif

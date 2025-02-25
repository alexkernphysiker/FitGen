// this file is distributed under
// LGPLv3 license
#ifndef RGCFVQNCAQKWQXFU
#define RGCFVQNCAQKWQXFU
#if __cplusplus<201100L
#error c++>=11 is needed for using math_h headers
#endif
#include <random>
#include <type_traits>
#include <math_h/error.h>
#include "abstract.h"
namespace Genetic
{
    class EmptyMutation :public virtual AbstractGenetic {
    public:
        EmptyMutation() {}
        EmptyMutation(EmptyMutation&& source) :AbstractGenetic(std::move(source)) {}
        virtual ~EmptyMutation() {}
    };
    template<class FITGEN = EmptyMutation>
    class DifferentialMutations : public FITGEN, public virtual AbstractGenetic
    {
        static_assert(std::is_base_of<EmptyMutation, FITGEN>::value, "Mutation algorithm must be a class derived from EmptyMutation");
    private:
        double M;
    public:
        DifferentialMutations() : FITGEN(), M(0.5) {}
        DifferentialMutations(DifferentialMutations&& source) : FITGEN(std::move(source)), M(std::move(source.M)), AbstractGenetic(std::move(source)) {}
        virtual ~DifferentialMutations() {}
        const double& MutationCoefficient()const
        {
            return M;
        }
        DifferentialMutations& SetMutationCoefficient(const double& val)
        {
            if (val < 0)
                throw MathTemplates::Exception<DifferentialMutations>("DifferentialMutations: mutation coefficient should be a positive value");
            M = val;
            return *this;
        }
    protected:
        virtual void mutations(ParamSet& C)const override
        {
            FITGEN::mutations(C);
            std::uniform_int_distribution<int> randomelement(0, AbstractGenetic::PopulationSize() - 1);
            const auto& A = AbstractGenetic::Parameters(randomelement(MathTemplates::RandomEngine<>::Instance()));
            const auto& B = AbstractGenetic::Parameters(randomelement(MathTemplates::RandomEngine<>::Instance()));
            for (size_t i = 0; i < C.size(); i++)
                C(i) += M * (A[i] - B[i]);
        }
    };
    template<class FITGEN = EmptyMutation>
    class Crossing : public FITGEN, public virtual AbstractGenetic
    {
        static_assert(std::is_base_of<EmptyMutation, FITGEN>::value, "Mutation algorithm must be a class derived from EmptyMutation");
    private:
        double P;
    public:
        Crossing() : FITGEN(), P(0) {}
        Crossing(Crossing&& source) : FITGEN(std::move(source)), P(std::move(source.P)), AbstractGenetic(std::move(source)) {}
        virtual ~Crossing() {}
        const double& CrossingProbability()const
        {
            return P;
        }
        Crossing& SetCrossingProbability(const double& val)
        {
            if ((val < 0) || (val > 1))
                throw MathTemplates::Exception<Crossing>("Crossing: probability value should fit the condition 0<=P<=1");
            P = val;
            return *this;
        }
    protected:
        virtual void mutations(ParamSet& C)const override
        {
            FITGEN::mutations(C);
            MathTemplates::RandomUniform<> Prob(0, 1);
            if (Prob() < P) {
                std::uniform_int_distribution<int> randomelement(0, AbstractGenetic::PopulationSize() - 1);
                auto X = AbstractGenetic::Parameters(randomelement(MathTemplates::RandomEngine<>::Instance()));
                FITGEN::mutations(X);
                for (size_t i = 0; i < C.size(); i++) {
                    if (Prob() < 0.5)
                        C(i) = X[i];
                }
            }
        }
    };
    template<class FITGEN = EmptyMutation>
    class AbsoluteMutations : public FITGEN, public virtual AbstractGenetic
    {
        static_assert(std::is_base_of<EmptyMutation, FITGEN>::value, "Mutation algorithm must be a class derived from EmptyMutation");
    private:
        double P;
        ParamSet m_mutation;
    public:
        AbsoluteMutations() : FITGEN(), P(0) {}
        AbsoluteMutations(AbsoluteMutations&& source) :
            FITGEN(std::move(source)), P(std::move(source.P)),
            m_mutation(std::move(source.m_mutation)),
            AbstractGenetic(std::move(source)) {
        }
        virtual ~AbsoluteMutations() {}
        const ParamSet& AbsoluteMutationCoefficients()const
        {
            return m_mutation;
        }
        AbsoluteMutations& SetAbsoluteMutationCoefficients(const ParamSet& p)
        {
            m_mutation = {};
            for (double v : p) {
                if (v < 0)
                    throw MathTemplates::Exception<AbsoluteMutations>("AbsoluteMutations: mutation coefficient cannot be negative");
                m_mutation.push_back(v);
            }
            return *this;
        }
        const double& AbsoluteMutationsProbability()const
        {
            return P;
        }
        AbsoluteMutations& SetAbsoluteMutationsProbability(const double& val)
        {
            if ((val < 0) || (val > 1))
                throw MathTemplates::Exception<AbsoluteMutations>("AbsoluteMutations: probability value should fit the condition 0<=P<=1");
            P = val;
            return *this;
        }
    protected:
        virtual void mutations(ParamSet& C)const override
        {
            FITGEN::mutations(C);
            MathTemplates::RandomUniform<> Prob(0, 1);
            if (Prob() < P)
                for (size_t i = 0; i < AbstractGenetic::ParamCount(); i++) {
                    MathTemplates::RandomGauss<> distr(0, m_mutation[i]);
                    C(i) += distr();
                }
        }
    };
    template<class FITGEN = EmptyMutation>
    class RelativeMutations : public FITGEN, public virtual AbstractGenetic
    {
        static_assert(std::is_base_of<EmptyMutation, FITGEN>::value, "Mutation algorithm must be a class derived from EmptyMutation");
    private:
        double P;
        ParamSet m_mutation;
    public:
        RelativeMutations() : FITGEN(), P(0) {}
        RelativeMutations(RelativeMutations&& source) :
            FITGEN(std::move(source)), P(std::move(source.P)),
            m_mutation(std::move(source.m_mutation)), AbstractGenetic(std::move(source)) {
        }
        virtual ~RelativeMutations() {}
        const ParamSet& RelativeMutationCoefficients()const
        {
            return m_mutation;
        }
        RelativeMutations& SetRelativeMutationCoefficients(const ParamSet& p)
        {
            m_mutation = {};
            for (double v : p) {
                if (v < 0)
                    throw MathTemplates::Exception<RelativeMutations>("AbsoluteMutations: mutation coefficient cannot be negative");
                m_mutation.push_back(v);
            }
            return *this;
        }
        const double& RelativeMutationsProbability()const
        {
            return P;
        }
        RelativeMutations& SetRelativeMutationsProbability(const double& val)
        {
            if ((val < 0) || (val > 1))
                throw MathTemplates::Exception<RelativeMutations>("RelativeMutations: probability value should fit the condition 0<=P<=1");
            P = val;
            return *this;
        }
    protected:
        virtual void mutations(ParamSet& C)const override
        {
            FITGEN::mutations(C);
            MathTemplates::RandomUniform<> Prob(0, 1);
            if (Prob() < P)
                for (size_t i = 0; i < AbstractGenetic::ParamCount(); i++) {
                    MathTemplates::RandomGauss<> distr(0, m_mutation[i]);
                    C(i) *= (1 + distr());
                }
        }
    };
    template<class FITGEN>
    class ExactCopying : public FITGEN, public virtual AbstractGenetic
    {
        static_assert(std::is_base_of<EmptyMutation, FITGEN>::value, "Mutation algorithm must be a class derived from EmptyMutation");
    private:
        double P;
    public:
        ExactCopying() : FITGEN(), P(0) {}
        ExactCopying(const ExactCopying& source) : FITGEN(source), P(source.P), AbstractGenetic(source) {}
        virtual ~ExactCopying() {}
        ExactCopying& SetExactCopyingProbability(const double& value)
        {
            if ((value < 0) || (value > 1))
                throw MathTemplates::Exception<ExactCopying>("ExactCopying: probability value should fit the condition 0<=P<=1");
            P = value;
            return *this;
        }
        const double& ExactCopyingProbability()const
        {
            return P;
        }
    protected:
        virtual void mutations(ParamSet& C)const override
        {
            MathTemplates::RandomUniform<> Prob(0, 1);
            if (Prob() > P)
                FITGEN::mutations(C);
        }
    };
}
#endif

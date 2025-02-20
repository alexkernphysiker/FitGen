// this file is distributed under
// LGPLv3 license
#include <math_h/error.h>
#include <Genetic/equation2.h>
namespace Genetic
{
    using namespace std;
    using namespace MathTemplates;
    InexactEquationSystem::InexactEquationSystem(const initializer_list<InexactEquation>& source)
    {
        for (const auto& item : source)f_data.push_back(item);
    }
    InexactEquationSystem::InexactEquationSystem(const vector<InexactEquation>& source)
    {
        for (const auto& item : source)f_data.push_back(item);
    }
    InexactEquationSystem::~InexactEquationSystem() {}
    double InexactEquationSystem::operator()(const ParamSet& P) const
    {
        double result = 0;
        for (const auto& item : f_data)
            result += item.right.NumCompare(item.left(P));
        return result;
    }
    const std::vector<InexactEquation>& Genetic::InexactEquationSystem::equations()const
    {
        return f_data;
    }
}

#ifndef ___DnoQpUEN
#define ___DnoQpUEN
#include "fit_gen.h"
namespace Fit{
	using namespace std;
	class GenerateUniform:public IInitialConditions{
	public:
		GenerateUniform();
		virtual ~GenerateUniform();
		int Count();
		double Min(int i);
		double Max(int i);
		GenerateUniform &Add(double min,double max);
		virtual ParamSet Generate()override;
	private:
		vector<double> m_min;
		vector<double> m_max;
	};
	class GenerateByGauss:public IInitialConditions{
	public:
		GenerateByGauss();
		virtual ~GenerateByGauss();
		int Count();
		double Mean(int i);
		double Sigma(int i);
		GenerateByGauss &Add(double mean,double sig);
		virtual ParamSet Generate()override;
	private:
		vector<double> m_mean;
		vector<double> m_sig;
	};
}
#endif
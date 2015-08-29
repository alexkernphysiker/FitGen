#include <iostream>
#include <fstream>
#include <random>
#include <fit.h>
#include <initialconditions.h>
#include <math_h/functions.h>
using namespace std;
using namespace Genetic;
int main(int argcnt, char **arg){
	double left=0;
	double right=10;
	unsigned int bins=2;
	int count=5000;
	auto distribution=make_shared<Distribution1D>(left,right,int(right-left)*bins);
	printf("Filling...\n");
	RANDOM engine;
	normal_distribution<double> gauss((right+left)/2.0,(right-left)/10.0);
	for(int i=0;i<count;i++)
		distribution->Fill(gauss(engine));
	printf("Prepare fitting...\n");
	Fit<DifferentialMutations<>,ChiSquareWithXError>
		fit(distribution,[](ParamSet&&X,ParamSet&&P){return Gaussian(X[0],P[0],P[1])*P[2];});
	fit.SetFilter([](ParamSet&&P){return (P[1]>0)&&(P[2]>0);});
	fit.Init(30,make_shared<GenerateUniform>()<<make_pair(left,right)<<make_pair(0,right-left)<<make_pair(0,2.0*count/bins),engine);
	printf("Parameter count: %i\n",fit.ParamCount());
	printf("Population size: %i\n",fit.PopulationSize());
	printf("Fitting...\n");
	while(!fit.AbsoluteOptimalityExitCondition(0.00000001)){
		fit.Iterate(engine);
		printf("%i iterations;   %f<=chi^2<=%f         \r",fit.iteration_count(),fit.Optimality(),fit.Optimality(fit.PopulationSize()-1));
	}
	printf("\nParameters:\n");
	for(double param:fit)
		printf("\t%f",param);
	printf("\nErrors:\n");
	for(double p:fit.GetParamParabolicErrors(parEq(fit.ParamCount(),0.01)))
		printf("\t%f",p);
	printf("\n");
	Plotter::Instance().SetOutput(".","distribution");
	PlotFit1D<decltype(fit)>().Points("Generated distribution",distribution).Fit("Fit distribution",fit)<<"set xlabel 'parameter value'\nset ylabel 'counts'";
	printf("Plot saved.\n");
	return 0;
}

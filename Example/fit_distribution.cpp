// this file is distributed under 
// MIT license
#include <iostream>
#include <fstream>
#include <random>
#include <gnuplot_wrap.h>
#include <math_h/functions.h>
#include <Genetic/fit.h>
#include <Genetic/initialconditions.h>
using namespace std;
using namespace Genetic;
using namespace GnuplotWrap;
using namespace MathTemplates;
int main(){
	RANDOM engine;
	double left=0;
	double right=10;
	int count=500;
	
	Distribution1D<double> distribution(BinsByStep(left,1.0,right));
	normal_distribution<double> gauss((right+left)/2.0,(right-left)/10.0);
	for(int i=0;i<count;i++)distribution.Fill(gauss(engine));
	
	Fit<DifferentialMutations<>,ChiSquareWithXError> fit(
		make_shared<FitPoints>(distribution),
		[](const ParamSet&X,const ParamSet&P){
			return Gaussian(X[0],P[0],P[1])*P[2];
		}
	);
	fit.SetFilter([](const ParamSet&P){return (P[1]>0)&&(P[2]>0);});
	fit.Init(30,
		make_shared<GenerateUniform>()
			<<make_pair(left,right)
			<<make_pair(0,right-left)
			<<make_pair(0,5.0*count)
		,engine
	);
	while(!fit.AbsoluteOptimalityExitCondition(0.00000001)){
		fit.Iterate(engine);
		cout<<fit.iteration_count()<<" iterations; "
			<<fit.Optimality()<<"<S<"
			<<fit.Optimality(fit.PopulationSize()-1)
			<<"        \r";
	}
	cout<<endl;
	cout<<"Fit parameters:"<<endl<<fit.Parameters()<<endl;
	cout<<"Fit parameters dispersion:"<<endl<<fit.ParamDispersion()<<endl;

	Plotter::Instance().SetOutput(".","distribution");
	Plot<double>().Hist(distribution)
	.Line(SortedPoints<double>(
		[&fit](double x)->double{return fit({x});},
		ChainWithStep(0.0,0.01,10.0)
	));

	return 0;
}

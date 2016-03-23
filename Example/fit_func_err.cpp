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
	
	FitFunctionWithError<DifferentialMutations<>,ChiSquareWithXError> fit(
		make_shared<FitPoints>(distribution),
		[](const ParamSet&X,const ParamSet&P){
			return value<double>::std_error(Gaussian(X[0],P[0],P[1])*P[2]);
		}
	);
	fit.SetFilter([](const ParamSet&P){return (P[1]>0)&&(P[2]>0);});
	fit.Init(30,make_shared<GenerateUniform>()<<make_pair(left,right)<<make_pair(0,right-left)<<make_pair(0,5.0*count),engine);
	while(!fit.AbsoluteOptimalityExitCondition(0.00000001))fit.Iterate(engine);

	cout<<"Fit parameters (most optimal values):"<<endl<<fit.Parameters()<<endl;
	fit.SetUncertaintyCalcDeltas({0.01,0.01,0.01});
	cout<<"Fit parameters (uncertainty)"<<endl;
	for(const auto&P:fit.ParametersWithUncertainties())cout<<P<<endl;
	
	Plotter::Instance().SetOutput(".","func_with_error");
	SortedPoints<value<double>> func_draw;
	for(double x:ChainWithStep(0.0,0.1,10.0))func_draw<<point<value<double>>(value<double>(x,0),fit({x}));
	Plot<double>().Hist(distribution).Hist(func_draw,"fit");
	
	return 0;
}

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

    //Generating the distribution for fitting
    //we use class from math_h submodule
    double left=0,right=10,sigma=1.0;size_t count=500;
    Distribution1D<double> distribution(BinsByStep(left,1.0,right));
    normal_distribution<double> gauss((right+left)/2.0,sigma);
    for(size_t i=0;i<count;i++)
	distribution.Fill(gauss(engine));

    //Fitting generated distribution
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
	     <<make_pair(0,5.0*count),
	 engine
    );
    while(!fit.AbsoluteOptimalityExitCondition(0.000001))
	fit.Iterate(engine);

    //Output results
    cout<<"Chi^2 = "<<fit.Optimality()<<endl;
    cout<<"Chi^2 divided by degrees of freedom = "
    <<fit.Optimality()/(fit.Points()->size()-fit.ParamCount())<<endl;
    cout<<endl;
    fit.SetUncertaintyCalcDeltas({0.01,0.01,0.01});//for estimating the uncertainties
    cout<<"Fit parameters"<<endl;
    for(const auto&P:fit.ParametersWithUncertainties())
	cout<<P<<endl;
    Plotter::Instance().SetOutput(".","gauss-fit");
    Plot<double>().Hist(distribution).Line(
	SortedPoints<double>(
	    [&fit](double x)->double{return fit({x});},
	     ChainWithStep(0.0,0.01,10.0)
	)
    );
    return 0;
}

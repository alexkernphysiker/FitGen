// this file is distributed under 
// MIT license
#include <iostream>
#include <fstream>
#include <gnuplot_wrap.h>
#include <Genetic/fit.h>
#include <Genetic/parabolic.h>
#include <Genetic/filter.h>
#include <Genetic/initialconditions.h>
using namespace std;
using namespace Genetic;
using namespace MathTemplates;
using namespace GnuplotWrap;
int main(){
	//points
	auto points_to_fit=make_shared<FitPoints>()
		<<Point({{-67.5,2.5}},{179.4,12.5})
		<<Point({{-62.5,2.5}},{213.1,13.0})
		<<Point({{-57.5,2.5}},{221.6,12.0})
		<<Point({{-52.5,2.5}},{220.1,11.5})
		<<Point({{-47.5,2.5}},{214.4,10.5})
		<<Point({{-42.5,2.5}},{241.0,10.5})
		<<Point({{-37.5,2.5}},{237.3,10.5})
		<<Point({{-32.5,2.5}},{240.1,10.0})
		<<Point({{-27.5,2.5}},{268.7,10.0})
		<<Point({{-22.5,2.5}},{316.3,10.0})
		<<Point({{-17.5,2.5}},{329.5,11.0})
		<<Point({{-12.5,2.5}},{317.2,10.5})
		<<Point({{-07.5,2.5}},{360.4,12.0})
		<<Point({{-02.5,2.5}},{371.7,11.5})
		<<Point({{ 02.5,2.5}},{406.6,12.0})
		<<Point({{ 07.5,2.5}},{413.2,12.0})
		<<Point({{ 12.5,2.5}},{455.1,13.0})
		<<Point({{ 17.5,2.5}},{495.3,14.0})
		<<Point({{ 22.5,2.5}},{497.3,14.5})
		<<Point({{ 27.5,2.5}},{511.4,15.0});
		
	//Foreground, background and total sum for fitting
	typedef Mul<Func3<Gaussian,Arg<0>,Par<2>,Par<1>>,Par<0>> Foreground;
	const int background_polynom_power=5;
	typedef PolynomFunc<0,Foreground::ParamCount,background_polynom_power> Background;
	typedef Add<Foreground,Background> TotalFunc;

	//Fitting
	RANDOM random_engine;
	FitFunction<DifferentialMutations<Uncertainty>,TotalFunc> fit(points_to_fit);
	fit.SetFilter(
	    make_shared<And>()
		<<(make_shared<Above>()<<0<<0)
		<<[](const ParamSet&P){
		    static Foreground F;
		    return F({P[2]},P)<P[1]*5.0;
		}
	);
	auto initial=make_shared<InitialDistributions>()
	    <<make_shared<DistribGauss>(100.,100.)
	    <<make_shared<DistribUniform>(0.,20.)
	    <<make_shared<FixParam>(-20.)
	    <<make_shared<DistribGauss>(400.,100.)
	    <<make_shared<DistribGauss>(5.,2.)
	    <<make_shared<DistribGauss>(0.,0.5);
	while(initial->Count()<TotalFunc::ParamCount)
	    initial<<make_shared<DistribGauss>(0.,0.01);
	fit.Init(TotalFunc::ParamCount*15,initial,random_engine);
	
	while(!fit.AbsoluteOptimalityExitCondition(0.0001)){
	    fit.Iterate(random_engine);
	    cout<<fit.iteration_count()<<" iterations; "
		<<fit.Optimality()<<" < Chi^2 < "
		<<fit.Optimality(fit.PopulationSize()-1)
		<<"           \r";
	}
	cout<<endl;
	
	//Output results
	cout<<"Chi^2 divided by degrees of freedom = "
	<<fit.Optimality()/(fit.Points()->size()-fit.ParamCount())
	<<endl;
	cout<<"Fit parameters with uncertainties"<<endl;
	fit.SetUncertaintyCalcDeltas(parEq(fit.ParamCount(),0.01));
	for(const auto&P:fit.ParametersWithUncertainties())
	    cout<<P<<endl;
	
	//Plotting total fit and background
	Plotter::Instance().SetOutput(".","foreground-background-fit");
	const auto&P=fit.Parameters();
	const auto chain=ChainWithStep(-70.0,0.1,30.0);
	const SortedPoints<double>
	    totalfit([&fit](double x)->double{return fit({x});},chain),
	    background([&P](double x)->double{return Background()({x},P);},chain);
	Plot<double>().Hist(points_to_fit->Hist1(0),"points").Line(totalfit,"fit")
	.Line(background,"background")<<"set key on";
	return 0;
}

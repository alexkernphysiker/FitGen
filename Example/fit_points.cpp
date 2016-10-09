// this file is distributed under 
// MIT license
#include <iostream>
#include <fstream>
#include <gnuplot_wrap.h>
#include <Genetic/fit.h>
#include <Genetic/filter.h>
#include <Genetic/initialconditions.h>
const int background_polynom_power=4;
using namespace std;
using namespace Genetic;
using namespace MathTemplates;
using namespace GnuplotWrap;
int main(){
	RANDOM engine;
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
		
	typedef Mul<Func3<Gaussian,Arg<0>,Par<2>,Par<1>>,Par<0>> Foreground;
	typedef PolynomFunc<0,Foreground::ParamCount,background_polynom_power> Background;
	typedef Add<Foreground,Background> TotalFunc;
	FitFunction<DifferentialMutations<>,TotalFunc,ChiSquareWithXError> fit(points_to_fit);
	fit.SetFilter(make_shared<And>()
		<<(make_shared<Above>()<<0<<0)
		<<[](const ParamSet&P){
			static Foreground F;
			return F({P[2]},P)<P[1]*5.0;
		}
	);
	auto initial=make_shared<GenerateByGauss>()
		<<make_pair(100.,100.)<<make_pair(30.,30.)<<make_pair(-20.,10.)
		<<make_pair(400.,100.)<<make_pair(5.,1.)<<make_pair(0.,0.5);
	while(initial->Count()<TotalFunc::ParamCount)
		initial<<make_pair(0.,0.01);
	fit.Init(TotalFunc::ParamCount*15,initial,engine);
	
	while(!fit.AbsoluteOptimalityExitCondition(0.0001)){
		fit.Iterate(engine);
		cout<<fit.iteration_count()<<" iterations; "
			<<fit.Optimality()<<" < Chi^2 < "
			<<fit.Optimality(fit.PopulationSize()-1)
			<<"           \r";
	}
	cout<<endl;
	cout<<"Chi^2 divided by degrees of freedom = "<<fit.Optimality()/(fit.Points()->size()-fit.ParamCount())<<endl;
	cout<<endl;
	cout<<"Fit parameters with uncertainties"<<endl;
	fit.SetUncertaintyCalcDeltas(parEq(fit.ParamCount(),0.01));
	for(const auto&P:fit.ParametersWithUncertainties())cout<<P<<endl;
	Plotter::Instance().SetOutput(".","points");
	SortedPoints<double>
		totalfit([&fit](double x)->double{return fit({x});},ChainWithStep(-70.0,0.1,30.0)),
		background([&fit](double x)->double{return Background()({x},fit.Parameters());},ChainWithStep(-70.0,0.1,30.0));
	Plot<double>().Hist(points_to_fit->Hist1(0)).Line(totalfit,"Fit").Line(background,"Background");
	return 0;
}

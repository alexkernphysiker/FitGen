// this file is distributed under 
// MIT license
#include <iostream>
#include <fstream>
#include <Genetic/fit.h>
#include <Genetic/filter.h>
#include <Genetic/initialconditions.h>
const int background_polynom_power=5;
using namespace std;
using namespace Genetic;
typedef Mul<Func3<BreitWigner,Arg<0>,Par<2>,Par<1>>,Par<0>> Foreground;
typedef PolynomFunc<0,Foreground::ParamCount,background_polynom_power> Background;
typedef Add<Foreground,Background> TotalFunc;
int main(){
	RANDOM engine;
	auto points_to_fit=make_shared<FitPoints>()
		<<Point({-67.5},{2.5},179.4,12.5)
		<<Point({-62.5},{2.5},223.1,13.0)
		<<Point({-57.5},{2.5},221.6,12.0)
		<<Point({-52.5},{2.5},220.1,11.5)
		<<Point({-47.5},{2.5},214.4,10.5)
		<<Point({-42.5},{2.5},241.0,10.5)
		<<Point({-37.5},{2.5},237.3,10.5)
		<<Point({-32.5},{2.5},230.1,10.0)
		<<Point({-27.5},{2.5},258.7,10.0)
		<<Point({-22.5},{2.5},290.3,10.5)
		<<Point({-17.5},{2.5},307.5,11.0)
		<<Point({-12.5},{2.5},317.2,10.5)
		<<Point({-07.5},{2.5},365.4,11.5)
		<<Point({-02.5},{2.5},371.7,11.5)
		<<Point({ 02.5},{2.5},406.6,12.0)
		<<Point({ 07.5},{2.5},413.2,12.0)
		<<Point({ 12.5},{2.5},455.1,13.0)
		<<Point({ 17.5},{2.5},495.3,14.0)
		<<Point({ 22.5},{2.5},497.3,14.5)
		<<Point({ 27.5},{2.5},511.4,15.0);
	FitFunction<DifferentialMutations<>,TotalFunc,ChiSquareWithXError> fit(points_to_fit);
	fit.SetFilter(make_shared<And>()
		<<(make_shared<Above>()<<0<<0)
		<<(make_shared<Below>()<<INFINITY<<5)
		<<[](const ParamSet&P){Foreground F;return F({P[2]},P)<P[1]*5.0;}
	);
	auto initial=make_shared<GenerateByGauss>()
		<<make_pair(100,100)<<make_pair(20,20)<<make_pair(-20,0)
		<<make_pair(400,100)<<make_pair(5,1)<<make_pair(0,0.5);
	while(initial->Count()<TotalFunc::ParamCount)
		initial<<make_pair(0,0.01);
	fit.Init(TotalFunc::ParamCount*20,initial,engine);
	cout<<"Population:"<<fit.PopulationSize()<<endl;
	cout<<"Parameters:"<<fit.ParamCount()<<endl;
	while(!(
		fit.AbsoluteOptimalityExitCondition(0.000001)&&
		fit.RelativeParametersDispersionExitCondition(parEq(TotalFunc::ParamCount,0.01))
	)){
		fit.Iterate(engine);
		cout<<fit.iteration_count()<<" iterations; "<<fit.Optimality()<<"<S<"<<fit.Optimality(fit.PopulationSize()-1)<<"        \r";
	}
	cout<<endl;
	cout<<"Fit parameters:"<<endl<<fit<<endl;
	cout<<"Fit parameters errors:"<<endl<<fit.GetParamParabolicErrors(parEq(fit.ParamCount(),0.001))<<endl;
	cout<<"Fit parameters averages:"<<endl<<fit.ParamAverage()<<endl;
	cout<<"Fit parameters dispersions:"<<endl<<fit.ParamDispersion()<<endl;

	Plotter::Instance().SetOutput(".","points");
	PlotFit1D<decltype(fit)>().Points("Generated distribution",points_to_fit).Fit("Fit distribution",fit,0.5)
		.ParamFunc("Foreground",Foreground(),fit,0.5).ParamFunc("Background",Background(),fit,0.5);
	return 0;
}

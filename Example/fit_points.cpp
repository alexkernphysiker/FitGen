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
		<<make_pair(make_pair(-67.5,2.5),make_pair(179.4,12.5))
		<<make_pair(make_pair(-62.5,2.5),make_pair(223.1,13.))
		<<make_pair(make_pair(-57.5,2.5),make_pair(221.6,12.))
		<<make_pair(make_pair(-52.5,2.5),make_pair(220.1,11.5))
		<<make_pair(make_pair(-47.5,2.5),make_pair(214.4,10.5))
		<<make_pair(make_pair(-42.5,2.5),make_pair(241.0,10.5))
		<<make_pair(make_pair(-37.5,2.5),make_pair(237.3,10.5))
		<<make_pair(make_pair(-32.5,2.5),make_pair(230.1,10.))
		<<make_pair(make_pair(-27.5,2.5),make_pair(258.7,10.))
		<<make_pair(make_pair(-22.5,2.5),make_pair(290.3,10.5))
		<<make_pair(make_pair(-17.5,2.5),make_pair(307.5,11.))
		<<make_pair(make_pair(-12.5,2.5),make_pair(317.2,10.5))
		<<make_pair(make_pair(-07.5,2.5),make_pair(365.4,11.5))
		<<make_pair(make_pair(-02.5,2.5),make_pair(371.7,11.5))
		<<make_pair(make_pair( 02.5,2.5),make_pair(406.6,12.))
		<<make_pair(make_pair( 07.5,2.5),make_pair(413.2,12.))
		<<make_pair(make_pair( 12.5,2.5),make_pair(455.1,13.))
		<<make_pair(make_pair( 17.5,2.5),make_pair(495.3,14.))
		<<make_pair(make_pair( 22.5,2.5),make_pair(497.3,14.5))
		<<make_pair(make_pair( 27.5,2.5),make_pair(511.4,15.));
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

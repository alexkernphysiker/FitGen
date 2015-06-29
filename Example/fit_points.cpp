#include <iostream>
#include <fstream>
#include <fit.h>
#include <paramfunc.h>
#include <filter.h>
#include <initialconditions.h>
const int background_polynom_power=6;
using namespace std;
using namespace Genetic;
typedef Mul<Func3<BreitWigner,Arg<0>,Par<2>,Par<1>>,Par<0>> Foreground;
typedef PolynomFunc<0,3,background_polynom_power> Background;
typedef Add<Foreground,Background> TotalFunc;

double X[]={-67.5,-62.5,-57.5,-52.5,-47.5,-42.5,-37.5,-32.5,-27.5,
	-22.5,-17.5,-12.5,-7.5,-2.5,2.5,7.5,12.5,17.5,22.5,27.5};
double dX[]={2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5,
	2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5};
double Y[]={179.355, 223.078, 211.644, 220.101, 214.429, 241.046, 237.26, 230.127, 258.658, 290.287, 307.54, 317.218,
	365.392, 371.688, 406.558, 413.209, 455.075, 495.286, 497.253, 511.347};
double dY[]={12.4159, 13.178, 11.8098, 11.4024, 10.555, 10.7758, 10.3217,  9.81692, 10.1093, 10.5735,
	10.7944, 10.679, 11.4741, 11.4566, 12.0048, 12.1155, 13.117, 14.1408, 14.3717, 15.0141};

int main(int argcnt, char **arg){
	auto points_to_fit=FitPointsXdXYdY(0,19,X,dX,Y,dY);
	FitFunction<DifferentialMutations<>,TotalFunc,ChiSquareWithXError> fit(points_to_fit);
	fit.SetFilter(make_shared<And>()
		<<(make_shared<Above>()<<0<<0)
		<<(make_shared<Below>()<<INFINITY<<5)
		<<[](ParamSet&&P){Foreground F;return F.F(ParamSet(P[2]),P)<P[1]*5.0;}
	);
	auto initial=make_shared<GenerateByGauss>()
		<<make_pair(100,100)<<make_pair(20,20)<<make_pair(-20,0)
		<<make_pair(400,100)<<make_pair(5,1)<<make_pair(0,0.5);
	while(initial->Count()<TotalFunc::ParamCount)
		initial<<make_pair(0,0.01);
	fit.Init(TotalFunc::ParamCount*20,initial);
	printf("Parameter count: %i\n",fit.ParamCount());
	printf("Population size: %i\n",fit.PopulationSize());
	while(!(
		fit.AbsoluteOptimalityExitCondition(0.000001)&&
		fit.RelativeParametersDispersionExitCondition(parEq(TotalFunc::ParamCount,0.01))
	)){
		fit.Iterate();
		printf("Iteration count: %i; %f <= chi^2 <= %f         \r",fit.iteration_count(),fit.Optimality(),fit.Optimality(fit.PopulationSize()-1));
	}
	printf("\nParameters:\n");
	for(double p:fit)
		printf("\t%f",p);
	printf("\nErrors:\n");
	for(double p:fit.GetParamParabolicErrors(parEq(fit.ParamCount(),0.01)))
		printf("\t%f",p);
	printf("\nAverage:\n");
	for(double p:fit.ParamAverage())
		printf("\t%f",p);
	printf("\nDispersion:\n");
	for(double p:fit.ParamDispersion())
		printf("\t%f",p);
	printf("\n");
	{//plot calculation results
		ofstream data;
		data.open("output.data.txt");
		if(data.is_open()){
			for(auto p:(*points_to_fit))
				data<<p.X[0]<<" "<<p.y<<" "<<p.WX[0]<<" "<<p.wy<<"\n";
			data.close();
		}
		ofstream out;
		ofstream outbg;
		ofstream outfg;
		out.open("output.txt");
		outbg.open("output.bg.txt");
		outfg.open("output.fg.txt");
		if(out.is_open()&&outbg.is_open()&&outfg.is_open()){
			Background bg_func;
			Foreground fg_func;
			for(double x=-70; x<=30; x+=0.5){
				out<<x<<" "<<fit(ParamSet(x))<<"\n";
				outbg<<x<<" "<<bg_func(ParamSet(x),fit.Parameters())<<"\n";
				outfg<<x<<" "<<fg_func(ParamSet(x),fit.Parameters())<<"\n";
			}
			out.close();
			outbg.close();
			outfg.close();
		}
		ofstream script;
		script.open(".plotscript.gp");
		if(script.is_open()){
			script << "plot ";
			script << "\"output.data.txt\" using 1:2:($1-$3):($1+$3):($2-$4):($2+$4) with xyerrorbars title \"data\"";
			script << ",\\\n";
			script << "\"output.txt\" w l title \"total fit\"";
			script << ",\\\n";
			script << "\"output.fg.txt\" w l title \"foreground\"";
			script << ",\\\n";
			script << "\"output.bg.txt\" w l title \"background\"";
			script << "\n";
			script << "\npause -1";
			script.close();
		}
		system("gnuplot .plotscript.gp");
	}
	return 0;
}

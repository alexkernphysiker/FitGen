#include <iostream>
#include <fstream>

#include <fit_gen.h>
#include <fitpoints.h>
#include <paramfunc.h>
#include <filter.h>
#include <initialconditions.h>
const int background_polynom_power=5;
using namespace std;
using namespace Fit;
typedef Func4<BreitWigner,Arg<0>,Par<0>,Par<1>,Par<2>> Foreground;
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
	auto points_to_fit=FitPointsXdXYdY<ChiSquareWithXError>(0,19,X,dX,Y,dY);
	FitGen fit(make_shared<TotalFunc>(),points_to_fit,THREADS_COUNT);
	auto initial_cond=make_shared<GenerateByGauss>()
		<<make_pair(1.0,20.0)<<make_pair(20.0,20.0)<<make_pair(-20.0,0)
		<<make_pair(300.0,300.0)<<make_pair(4.0,4.0);
	while(initial_cond->Count()<TotalFunc::ParamCount)
		initial_cond<<make_pair(0.0,0.01);
	fit.SetFilter(make_shared<FilterAnd>()
		<<make_shared<FilterAbove>(ParamSet(0,7,INFINITY,0,0))
		<<make_shared<FilterBelow>(ParamSet(INFINITY,40))
	);
	fit.Init(TotalFunc::ParamCount*10,initial_cond);
	printf("Population size: %i\n",fit.PopulationSize());
	while(!fit.AbsoluteOptimalityExitCondition(0.0000001)){
		fit.Iterate();
		printf("%f <= chi^2 <= %f     \r",fit.Optimality(),fit.Optimality(fit.PopulationSize()-1));
	}
	printf("Iteration count: %i           \nchi^2 = %f\n",fit.iteration_count(),fit.Optimality());
	printf("par\t\toptimal\t\tParabolicErr\t\tmax_dev\t\taverage\t\tdisp\n");
	ParamSet err=fit.ParamParabolicError(parEq(fit.ParamCount(),0.001));
	for(int i=0; i<fit.ParamCount();i++)
		printf("par%i \t\t%f \t\t%f \t\t%f \t\t%f \t\t%f \n",i,
			   fit[i],err[i],fit.ParamMaxDeviation()[i],fit.ParamAverage()[i],fit.ParamDispersion()[i]);
	{//plot calculation results
		ofstream data;
		data.open("output.data.txt");
		if(data.is_open()){
			for(int i=0;i<points_to_fit->Count();i++)
				data<<points_to_fit->X(i)[0]<<" "
					<<points_to_fit->Y(i)<<" "
					<<points_to_fit->X_w(i)[0]<<" "
					<<points_to_fit->W(i)<<"\n";
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
				ParamSet X(x);
				out<<x<<" "<<fit(X)<<"\n";
				outbg<<x<<" "<<bg_func(X,fit.Parameters())<<"\n";
				outfg<<x<<" "<<fg_func(X,fit.Parameters())<<"\n";
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

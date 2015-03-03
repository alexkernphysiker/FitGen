#include <iostream>
#include <fstream>

#include <fit_gen.h>
#include <fitpoints.h>
#include <paramfunc.h>
#include <filter.h>
#include <initialconditions.h>
using namespace std;
using namespace Fit;
const unsigned char threads=1;
typedef Func4<BreitWigner,Arg<0>,Par<0>,Par<1>,Par<2>> Foreground;
typedef PolynomFunc<0,3,4> Background;

int main(int argcnt, char **arg){
	auto points_to_fit=make_shared<chi_2_wx>();
	points_to_fit->
		 Add(ParamSet(-67.5),ParamSet(2,5), 179.355, 12.4159)
		.Add(ParamSet(-62.5),ParamSet(2,5), 223.078, 13.178)
		.Add(ParamSet(-57.5),ParamSet(2,5), 211.644, 11.8098)
		.Add(ParamSet(-52.5),ParamSet(2,5), 220.101, 11.4024)
		.Add(ParamSet(-47.5),ParamSet(2,5), 214.429, 10.555)
		.Add(ParamSet(-42.5),ParamSet(2,5), 241.046, 10.7758)
		.Add(ParamSet(-37.5),ParamSet(2,5), 237.26 , 10.3217)
		.Add(ParamSet(-32.5),ParamSet(2,5), 230.127,  9.81692)
		.Add(ParamSet(-27.5),ParamSet(2,5), 258.658, 10.1093)
		.Add(ParamSet(-22.5),ParamSet(2,5), 290.287, 10.5735)
		.Add(ParamSet(-17.5),ParamSet(2,5), 307.54 , 10.7944)
		.Add(ParamSet(-12.5),ParamSet(2,5), 317.218, 10.679)
		.Add(ParamSet(- 7.5),ParamSet(2,5), 365.392, 11.4741)
		.Add(ParamSet(- 2.5),ParamSet(2,5), 371.688, 11.4566)
		.Add(ParamSet(  2.5),ParamSet(2,5), 406.558, 12.0048)
		.Add(ParamSet(  7.5),ParamSet(2,5), 413.209, 12.1155)
		.Add(ParamSet( 12.5),ParamSet(2,5), 455.075, 13.117)
		.Add(ParamSet( 17.5),ParamSet(2,5), 495.286, 14.1408)
		.Add(ParamSet( 22.5),ParamSet(2,5), 497.253, 14.3717)
		.Add(ParamSet( 27.5),ParamSet(2,5), 511.347, 15.0141);
	FitGen fit(make_shared<Add<Foreground,Background>>(),points_to_fit);

	auto initial_cond=make_shared<GenerateByGauss>();
	initial_cond->Add(1,20).Add(20,20).Add(-20,0).Add(300,300).Add(4,4).Add(0,0.01).Add(0,0.01).Add(0,0.01);
	fit.Init(1000,initial_cond);

	auto filter=make_shared<FilterRangeIn>();
	filter->Add(0,30).Add(5,50).Add(-100,0).Add(0,1000).Add(0,10).Add(-1,1).Add(-0.1,0.1).Add(-0.1,0.1);
	fit.SetFilter(filter);
	fit.SetMutation(Fit::mutDifferential,parEq(7,0.9)<<1);
	fit.SetMutation(Fit::mutRatio,ParamSet(0.01,0,0,0.1)<<parEq(2,0.05)<<0.1<<0.2);
	fit.SetMutation(Fit::mutAbsolute,parZeros(8));

	do{
		fit.Iterate(threads);
		printf("%f <= chi^2 <= %f     \r",fit.GetOptimality(),fit.GetOptimality(fit.PopulationSize()-1));
	}while (fit.GetOptimality(fit.PopulationSize()-1)>fit.GetOptimality());
	printf("Iteration count: %i           \nchi^2 = %f\n",fit.iteration_count(),fit.GetOptimality());

	ParamSet delta;
	for(int i=0; i<fit.ParamCount();i++)
		delta<<((fit.ParamDispersion()[i]>0.1)?fit.ParamDispersion()[i]:0.1);
	ParamSet er=fit.ParamParabolicError(delta);
	for(int i=0; i<fit.ParamCount();i++)
		printf("par%i = %f +/- %f\n",i,fit[i],er[i]);

	{
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
		if(
			out.is_open()&&
			outbg.is_open()&&
			outfg.is_open()
		){
			Background bg_func;
			Foreground fg_func;
			ParamSet P=fit.GetParameters();
			for(double x=-70; x<=30; x+=1){
				ParamSet X(x);
				out<<x<<" "<<fit(X)<<"\n";
				outbg<<x<<" "<<bg_func(X,P)<<"\n";
				outfg<<x<<" "<<fg_func(X,P)<<"\n";
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

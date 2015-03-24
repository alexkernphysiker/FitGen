#include <iostream>
#include <fstream>
#include <paramfunc.h>
#include <fitpoints.h>
#include <genetic.h>
#include <filter.h>
#include <initialconditions.h>
using namespace std;
using namespace Fit;
int main(int argcnt, char **arg){
	auto points_to_fit=make_shared<Distribution1D<ChiSquareWithXError>>(0,10,20);
	for(int i=0;i<200;i++)
		points_to_fit->AddValue(RandomGauss(1.0,5.0));
	DifferentialRandomMutations<>
		fit(make_shared<Mul<Func3<Gaussian,Arg<0>,Par<0>,Par<1>>,Par<2>>>(),points_to_fit,THREADS_COUNT);
	fit.SetFilter(make_shared<Above>(ParamSet(INFINITY,0,0)));
	fit.Init(30,make_shared<GenerateByGauss>()<<make_pair(4,2)<<make_pair(1,1)<<make_pair(500,500));
	printf("Parameter count: %i\n",fit.ParamCount());
	printf("Population size: %i\n",fit.PopulationSize());
	while(!fit.AbsoluteOptimalityExitCondition(0.000001)){
		fit.Iterate();
		printf("%f <= chi^2 <= %f     \r",fit.Optimality(),fit.Optimality(fit.PopulationSize()-1));
	}
	printf("%i iterations                     \n",fit.iteration_count());
	printf("chi^2 = %f\n",fit.Optimality());
	for(int i=0; i<fit.ParamCount();i++)
		printf("par%i=%f\t",i,fit[i]);
	printf("\n");
	
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
		out.open("output.txt");
		if(out.is_open()){
			for(double x=0; x<=10; x+=0.1)
				out<<x<<" "<<fit(ParamSet(x))<<"\n";
			out.close();
		}
		ofstream script;
		script.open(".plotscript.gp");
		if(script.is_open()){
			script << "plot ";
			script << "\"output.data.txt\" using 1:2:($1-$3):($1+$3):($2-$4):($2+$4) with xyerrorbars title \"data\"";
			script << ",\\\n";
			script << "\"output.txt\" w l title \"total fit\"";
			script << "\n";
			script << "\npause -1";
			script.close();
		}
		system("gnuplot .plotscript.gp");
	}
	return 0;
}

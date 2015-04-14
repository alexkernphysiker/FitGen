#include <iostream>
#include <fstream>
#include <fit.h>
#include <initialconditions.h>
#include <math_h/functions.h>
using namespace std;
using namespace Genetic;
int main(int argcnt, char **arg){
	double left=0;
	double right=10;
	unsigned int bins=2;
	int count=500;
	auto distribution=make_shared<Distribution1D>(left,right,int(right-left)*bins);
	printf("Filling...\n");
	for(int i=0;i<count;i++)
		distribution->Fill(RandomGauss((right-left)/10.0)+(right+left)/2.0);
	printf("Prepare fitting...\n");
	Fit<DifferentialMutations<>> fit(distribution,[](ParamSet&X,ParamSet&P){return Gaussian(X[0],P[0],P[1])*P[2];},ChiSquareWithXError);
	fit.SetFilter([](ParamSet& P){return (P[1]>0)&&(P[2]>0);});
	fit.Init(30,make_shared<Initialiser>()
		<<[left,right](){return RandomUniformlyR(left,right);}
		<<[left,right](){return RandomUniformlyR(0.0,right-left);}
		<<[count,bins](){return RandomGauss(double(count/bins),double(count/bins));}
	);
	printf("Parameter count: %i\n",fit.ParamCount());
	printf("Population size: %i\n",fit.PopulationSize());
	printf("Fitting...\n");
	while(!fit.AbsoluteOptimalityExitCondition(0.0000001)){
		fit.Iterate();
		printf("%i iterations;   %f<=chi^2<=%f         \r",fit.iteration_count(),fit.Optimality(),fit.Optimality(fit.PopulationSize()-1));
	}
	printf("\nParameters:\n");
	for(double param:fit)
		printf("\t%f",param);
	printf("\nErrors:\n");
	for(double p:fit.GetParamParabolicErrors(parEq(fit.ParamCount(),0.01)))
		printf("\t%f",p);
	printf("\n");
	{//plot calculation results
		ofstream data;
		data.open("output.data.txt");
		if(data.is_open()){
			for(auto p:*distribution)
				data<<p.X[0]<<" "<<p.y<<" "<<p.WX[0]<<" "<<p.wy<<"\n";
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

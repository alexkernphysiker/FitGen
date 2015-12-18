// this file is distributed under 
// MIT license
#include <iostream>
#include <fstream>
#include <random>
#include <math_h/functions.h>
#include <Genetic/fit.h>
#include <Genetic/initialconditions.h>
using namespace std;
using namespace Genetic;
int main(){
	double left=0;
	double right=10;
	unsigned int bins=2;
	int count=5000;
	auto distribution=make_shared<Distribution1D>(left,right,int(right-left)*bins);
	RANDOM engine;
	normal_distribution<double> gauss((right+left)/2.0,(right-left)/10.0);
	for(int i=0;i<count;i++)
		distribution->Fill(gauss(engine));
	printf("Prepare fitting...\n");
	Fit<DifferentialMutations<>,ChiSquareWithXError> fit(distribution,[](const ParamSet&X,const ParamSet&P){return Gaussian(X[0],P[0],P[1])*P[2];});
	fit.SetFilter([](const ParamSet&P){return (P[1]>0)&&(P[2]>0);});
	fit.Init(30,make_shared<GenerateUniform>()<<make_pair(left,right)<<make_pair(0,right-left)<<make_pair(0,2.0*count/bins),engine);
	cout<<"Population:"<<fit.PopulationSize()<<endl;
	cout<<"Parameters:"<<fit.ParamCount()<<endl;
	while(!fit.AbsoluteOptimalityExitCondition(0.00000001)){
		fit.Iterate(engine);
		cout<<fit.iteration_count()<<" iterations; "<<fit.Optimality()<<"<S<"<<fit.Optimality(fit.PopulationSize()-1)<<"        \r";
	}
	cout<<endl<<"Fit parameters:"<<endl;
	for(double x:fit)cout<<x<<"\t";
	cout<<endl<<"Fit parameters dispersions:"<<endl;
	for(double x:fit.ParamDispersion())cout<<x<<"\t";
	cout<<endl;
	Plotter::Instance().SetOutput(".","distribution");
	PlotFit1D<decltype(fit)>().Fit("Fit distribution","Generated distribution",fit,0.1)<<"set xlabel 'argument value'\nset ylabel 'counts'";
	return 0;
}

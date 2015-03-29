#include <iostream>
#include <fstream>
#include <fitpoints.h>
#include <genetic.h>
#include <paramfunc.h>
#include <initialconditions.h>
using namespace std;
using namespace Fit;
int main(int argcnt, char **arg){
	double left=0;
	double right=10;
	unsigned int bins=2;
	int count=500;
	
	auto points_to_fit=make_shared<Distribution1D<ChiSquareWithXError>>(left,right,int(right-left)*bins);
	for(int i=0;i<count;i++)
		points_to_fit->Fill(RandomGauss((right-left)/10.0)+(right+left)/2.0);
	
	DifferentialRandomMutations<> fit(make_shared<ParameterFunction<>>(
			[](ParamSet& X,ParamSet& P){return Gaussian(X[0],P[0],P[1])*P[2];},
			[](ParamSet& P){return (P[1]>0)&&(P[2]>0);}
		),points_to_fit,THREADS_COUNT);
	fit.Init(30,make_shared<Initialiser>()
		<<[left,right](){return RandomUniformly(left,right);}
		<<[left,right](){return RandomUniformly(0.0,right-left);}
		<<[count,bins](){return RandomGauss(double(count/bins),double(count/bins)*0.5);}
	);
	while(!fit.ConcentratedInOnePoint())
		fit.Iterate();
	printf("Done in %i iterations\n",fit.iteration_count());
	printf("chi^2 = %f\n",fit.Optimality());
	for(int i=0; i<fit.ParamCount();i++)
		printf("par%i=%f\n",i,fit[i]);
	
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

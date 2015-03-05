#include <math.h>
#include "fitpoints.h"
#include "fitexception.h"
namespace Fit{
	FitPointsAbstract::~FitPointsAbstract(){
		m_data.clear();
		m_data_w.clear();
		m_y.clear();
		m_w.clear();
	}
	FitPointsAbstract &FitPointsAbstract::Add(ParamSet x, ParamSet x_w, double y, double weight){
		if(weight<=0)
			throw new FitException("AproPoints: wieght must be a positive value");
		m_data.push_back(x);
		m_data_w.push_back(x_w);
		m_y.push_back(y);
		m_w.push_back(weight);
		return *this;
	}
	FitPointsAbstract &FitPointsAbstract::Add(ParamSet x, double y, double weight){
		return Add(x,ParamSet(),y,weight);
	}
	int FitPointsAbstract::Count(){return m_data.size();}
	ParamSet FitPointsAbstract::X(int i){
		if((i<0)||(i>=m_data.size()))
			throw new FitException("AproPoints: attempt to get point out of range");
		return m_data[i];
	}
	ParamSet FitPointsAbstract::X_w(int i){
		if((i<0)||(i>=m_data_w.size()))
			throw new FitException("AproPoints: attempt to get point out of range");
		return m_data_w[i];
	}
	double FitPointsAbstract::Y(int i){
		if((i<0)||(i>=m_y.size()))
			throw new FitException("AproPoints: getting point out of range");
		return m_y[i];
	}
	double FitPointsAbstract::W(int i){
		if((i<0)||(i>=m_w.size()))
			throw new FitException("AproPoints: getting point out of range");
		return m_w[i];
	}
	
	double SquareDiff::operator ()(ParamSet params, IParamFunc &func){
		double res=0;
		for(int i=0; i<Count();i++)
			res+=::pow(Y(i)-func(X(i),params),2)*W(i);
		return res;
	}
	
	double chi_2::operator ()(ParamSet params, IParamFunc &func){
		double z=Count()-params.Count();
		if(z<=0)throw new FitException("wrong conditions for calculating xi^2");
		double res=0;
		for(int i=0; i<Count();i++)
			res+=::pow((Y(i)-func(X(i),params))/W(i),2);
		return res/z;
	}
	
	double chi_2_wx::operator ()(ParamSet params, IParamFunc &func){
		double z=Count()-params.Count();
		if(z<=0)throw new FitException("wrong conditions for calculating xi^2");
		double res=0;
		for(int i=0; i<Count();i++){
			double w=W(i);	w*=w;
			for(int j=0; (j<X(i).Count())&&(j<X_w(i).Count());j++){
				ParamSet x1=X(i);
				ParamSet x2=X(i);
				x1.Set(j,X(i)[j]+X_w(i)[j]);
				x2.Set(j,X(i)[j]-X_w(i)[j]);
				w+=::pow(0.5*(func(x1,params)-func(x2,params)),2);
			}
			res+=::pow((Y(i)-func(X(i),params)),2)/w;
		}
		return res/z;
	}
}

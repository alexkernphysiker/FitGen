#include <iostream>
#include <math.h>
#include "fitpoints.h"
#include "fitexception.h"
using namespace std;
namespace Fit{
	FitPoints::Point::Point(){}
	FitPoints::Point::Point(const Point& src){
		X=src.X;
		WX=src.WX;
		y=src.y;
		wy=src.wy;
	}
	FitPoints::Point& FitPoints::Point::operator=(const Point& src){
		X=src.X;
		WX=src.WX;
		y=src.y;
		wy=src.wy;
		return *this;
	}
	FitPoints::FitPoints(){}
	FitPoints::~FitPoints(){}
	int FitPoints::count(){
		return m_data.size();
	}
	FitPoints& FitPoints::operator<<(Point point){
		m_data.push_back(point);
		return *this;
	}
	shared_ptr<FitPoints> operator<<(shared_ptr<FitPoints>src, FitPoints::Point p){
		src->operator<<(p);
		return src;
	}
	shared_ptr< FitPoints > operator<<(shared_ptr<FitPoints> src, pair<double,double> p){
		FitPoints::Point P;
		P.X<<p.first;
		P.y=p.second;
		P.wy=1;
		return src<<P;
	}
	FitPoints::Point& FitPoints::operator[](int i){
		return m_data[i];
	}
	FitPoints::iterator FitPoints::begin(){
		return m_data.begin();
	}
	FitPoints::const_iterator FitPoints::cbegin() const{
		return m_data.cbegin();
	}
	FitPoints::iterator FitPoints::end(){
		return m_data.end();
	}
	FitPoints::const_iterator FitPoints::cend() const{
		return m_data.cend();
	}
	shared_ptr<FitPoints> SelectFitPoints(shared_ptr<FitPoints> src, shared_ptr<IParamCheck> condition){
		auto res=make_shared<FitPoints>();
		for(FitPoints::Point p:(*src))
			if(condition->CorrectParams(p.X))
				res<<p;
			return res;
	}
	shared_ptr<FitPoints> SelectFitPoints(shared_ptr<FitPoints> src, function<bool(double)> Ycond){
		auto res=make_shared<FitPoints>();
		for(FitPoints::Point p:(*src))
			if(Ycond(p.y))
				res<<p;
			return res;
	}
	shared_ptr<FitPoints> SelectFitPoints(shared_ptr<FitPoints> src, shared_ptr<IParamCheck> condition,function<bool(double)> Ycond){
		auto res=make_shared<FitPoints>();
		for(FitPoints::Point p:(*src))
			if(condition->CorrectParams(p.X))
				res<<p;
			return res;
	}
	Distribution1D::Distribution1D(double min, double max, int bins):FitPoints(){
		if(max<min)
			throw new FitException("Wrong distribution parameters");
		if(bins<1)
			throw new FitException("Wrong distribution parameters");
		double binwidth=(max-min)/double(bins);
		double halfwidth=binwidth/2.0;
		for(double x=min+halfwidth;x<max;x+=binwidth){
			Point P;
			P.X<<x;
			P.WX<<halfwidth;
			P.y=0;
			P.wy=1;
			operator<<(P);
		}
	}
	void Distribution1D::Fill(double x){
		for(Point&p:(*this))
			if((x>=(p.X[0]-p.WX[0]))&&(x<(p.X[0]+p.WX[0]))){
				p.y=p.y+1.0;
				p.wy=sqrt(p.y);
			}
	}
	OptimalityForPoints::OptimalityForPoints(std::shared_ptr< FitPoints > p, OptimalityForPoints::Coefficient c, OptimalityForPoints::Summand s){
		points=p;
		C=c;
		S=s;
	}
	OptimalityForPoints::~OptimalityForPoints(){}
	double OptimalityForPoints::operator()(ParamSet& P, IParamFunc& F){
		double res=0;
		for(FitPoints::Point p:*points)
			res+=S(p,P,F);
		return res*C(P,F);
	}
	shared_ptr<OptimalityForPoints> SumSquareDiff(shared_ptr<FitPoints> points){
		OptimalityForPoints::Coefficient c=[](ParamSet&,IParamFunc&){
			return 1.0;
		};
		OptimalityForPoints::Summand s=[](FitPoints::Point&p,ParamSet&P,IParamFunc&F){
			return pow(p.y-F(p.X,P),2);
		};
		return make_shared<OptimalityForPoints>(points,c,s);
	}
	shared_ptr<OptimalityForPoints> SumWeightedSquareDiff(shared_ptr<FitPoints> points){
		OptimalityForPoints::Coefficient c=[points](ParamSet&,IParamFunc&){
			double z=0;
			for(FitPoints::Point p:*points)
				z+=p.wy;
			return 1.0/z;
		};
		OptimalityForPoints::Summand s=[](FitPoints::Point&p,ParamSet&P,IParamFunc&F){
			return pow(p.y-F(p.X,P),2)*p.wy;
		};
		return make_shared<OptimalityForPoints>(points,c,s);
	}
	shared_ptr<OptimalityForPoints> ChiSquare(shared_ptr<FitPoints> points){
		OptimalityForPoints::Coefficient c=[points](ParamSet&P,IParamFunc&){
			double z=points->count()-P.Count();
			if(z<=0)
				throw new FitException("wrong conditions for calculating xi^2: there must be at least one degree of freedom");
			return 1.0/z;
		};
		OptimalityForPoints::Summand s=[](FitPoints::Point&p,ParamSet&P,IParamFunc&F){
			return pow((p.y-F(p.X,P))/p.wy,2);
		};
		return make_shared<OptimalityForPoints>(points,c,s);
	}
	shared_ptr<OptimalityForPoints> ChiSquareWithXError(shared_ptr<FitPoints> points){
		OptimalityForPoints::Coefficient c=[points](ParamSet&P,IParamFunc&){
			double z=points->count()-P.Count();
			if(z<=0)
				throw new FitException("wrong conditions for calculating xi^2: there must be at least one degree of freedom");
			return 1.0/z;
		};
		OptimalityForPoints::Summand s=[](FitPoints::Point&p,ParamSet&P,IParamFunc&F){
			double w=pow(p.wy,2);
			for(int j=0; (j<p.X.Count())&&(j<p.WX.Count());j++){
				ParamSet x1=p.X;
				ParamSet x2=p.X;
				x1.Set(j,p.X[j]+p.WX[j]);
				x2.Set(j,p.X[j]-p.WX[j]);
				w+=pow(0.5*(F(x1,P)-F(x2,P)),2);
			}
			return pow((p.y-F(p.X,P)),2)/w;
		};
		return make_shared<OptimalityForPoints>(points,c,s);
	}
}

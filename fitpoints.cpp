#include <iostream>
#include <math.h>
#include "fitpoints.h"
#include "fitexception.h"
using namespace std;
namespace Fit{
	FitPoints::DataPoint::DataPoint(){}
	FitPoints::DataPoint::DataPoint(const DataPoint& src){
		X=src.X;
		WX=src.WX;
		y=src.y;
		wy=src.wy;
	}
	FitPoints::DataPoint& FitPoints::DataPoint::operator=(const DataPoint& src){
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
	FitPoints& FitPoints::operator<<(DataPoint point){
		m_data.push_back(point);
		return *this;
	}
	shared_ptr<FitPoints> operator<<(shared_ptr<FitPoints>src, FitPoints::DataPoint p){
		src->operator<<(p);
		return src;
	}
	shared_ptr< FitPoints > operator<<(shared_ptr<FitPoints> src, pair<double,double> p){
		FitPoints::DataPoint P;
		P.X<<p.first;
		P.y=p.second;
		P.wy=1;
		return src<<P;
	}
	FitPoints::DataPoint& FitPoints::operator[](int i){
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
		for(FitPoints::DataPoint p:(*src))
			if(condition->CorrectParams(p.X))
				res<<p;
			return res;
	}
	shared_ptr<FitPoints> SelectFitPoints(shared_ptr<FitPoints> src, function<bool(double)> Ycond){
		auto res=make_shared<FitPoints>();
		for(FitPoints::DataPoint p:(*src))
			if(Ycond(p.y))
				res<<p;
			return res;
	}
	shared_ptr<FitPoints> SelectFitPoints(shared_ptr<FitPoints> src, shared_ptr<IParamCheck> condition,function<bool(double)> Ycond){
		auto res=make_shared<FitPoints>();
		for(FitPoints::DataPoint p:(*src))
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
			DataPoint P;
			P.X<<x;
			P.WX<<halfwidth;
			P.y=0;
			P.wy=1;
			operator<<(P);
		}
	}
	void Distribution1D::Fill(double x){
		for(DataPoint&p:(*this))
			if((x>=(p.X[0]-p.WX[0]))&&(x<(p.X[0]+p.WX[0]))){
				p.y=p.y+1.0;
				p.wy=sqrt(p.y);
			}
	}
	OptimalityForPoints::~OptimalityForPoints(){}
	class _square_diff:public OptimalityForPoints{
	public:
		_square_diff(shared_ptr<FitPoints> p){
			points=p;
		}
		virtual ~_square_diff(){}
		virtual double operator()(ParamSet&P,IParamFunc&F)override{
			double res=0;
			for(FitPoints::DataPoint p:*points)
				res+=pow(p.y-F(p.X,P),2)*p.wy;
			return res;
		}
	};
	shared_ptr<OptimalityForPoints> SquareDiff(shared_ptr<FitPoints> points){
		return make_shared<_square_diff>(points);
	}
	class _chi_square:public OptimalityForPoints{
	public:
		_chi_square(shared_ptr<FitPoints> p){
			points=p;
		}
		virtual ~_chi_square(){}
		virtual double operator()(ParamSet&P,IParamFunc&F)override{
			double z=points->count()-P.Count();
			if(z<=0)
				throw new FitException("wrong conditions for calculating xi^2: there must be at least one degree of freedom");
			double res=0;
			for(FitPoints::DataPoint p:*points){
				res+=pow((p.y-F(p.X,P))/p.wy,2);
			}
			return res/z;
		}
	};
	shared_ptr<OptimalityForPoints> ChiSquare(shared_ptr<FitPoints> points){
		return make_shared<_chi_square>(points);
	}
	class _chi_square_wxe:public OptimalityForPoints{
	public:
		_chi_square_wxe(shared_ptr<FitPoints> p){
			points=p;
		}
		virtual ~_chi_square_wxe(){}
		virtual double operator()(ParamSet&P,IParamFunc&F)override{
			double z=points->count()-P.Count();
			if(z<=0)
				throw new FitException("wrong conditions for calculating xi^2: there must be at least one degree of freedom");
			double res=0;
			for(FitPoints::DataPoint p:*points){
				double w=pow(p.wy,2);
				for(int j=0; (j<p.X.Count())&&(j<p.WX.Count());j++){
					ParamSet x1=p.X;
					ParamSet x2=p.X;
					x1.Set(j,p.X[j]+p.WX[j]);
					x2.Set(j,p.X[j]-p.WX[j]);
					w+=pow(0.5*(F(x1,P)-F(x2,P)),2);
				}
				res+=pow((p.y-F(p.X,P)),2)/w;
			}
			return res/z;
		}
	};
	shared_ptr<OptimalityForPoints> ChiSquareWithXError(shared_ptr<FitPoints> points){
		return make_shared<_chi_square_wxe>(points);
	}
}

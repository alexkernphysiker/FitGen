// this file is distributed under 
// GPL v 3.0 license
#include <math.h>
#include "fit.h"
#include "genetic_exception.h"
using namespace std;
namespace Genetic{
	ParameterFunction::ParameterFunction(function<double(ParamSet&&,ParamSet&&)> f){
		func=f;
	}
	ParameterFunction::~ParameterFunction(){}
	double ParameterFunction::operator()(ParamSet&&X, ParamSet&&P){
		return func(static_cast<ParamSet&&>(X),static_cast<ParamSet&&>(P));
	}
	
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
		if((i<0)||(i>=count()))
			throw GeneticException("Range check error when getting an element from FitPoints");
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
			if(condition->operator()(static_cast<ParamSet&&>(p.X)))
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
			if(condition->operator()(static_cast<ParamSet&&>(p.X))&&Ycond(p.y))
				res<<p;
			return res;
	}
	Distribution1D::Distribution1D(double min, double max, int bins):FitPoints(){
		if((max<min)||(bins<1))
			throw GeneticException("Wrong constructor parameters for Distribution1D");
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
	
	OptimalityForPoints::OptimalityForPoints(
		std::shared_ptr< FitPoints > p, 
		shared_ptr<IParamFunc> f,
		OptimalityForPoints::Coefficient c, 
		OptimalityForPoints::Summand s
	){
		points=p;
		func=f;
		C=c;
		S=s;
	}
	OptimalityForPoints::~OptimalityForPoints(){}
	double OptimalityForPoints::operator()(ParamSet&&P){
		double res=0;
		for(FitPoints::Point p:*points)
			res+=S(p,P,*func);
		return res*C(P,*func);
	}
	shared_ptr<IOptimalityFunction> SumSquareDiff(shared_ptr<FitPoints> points, shared_ptr<IParamFunc> f){
		OptimalityForPoints::Coefficient c=[](ParamSet&,IParamFunc&){
			return 1.0;
		};
		OptimalityForPoints::Summand s=[](FitPoints::Point&p,ParamSet&P,IParamFunc&F){
			return pow(p.y-F(static_cast<ParamSet&&>(p.X),static_cast<ParamSet&&>(P)),2);
		};
		return make_shared<OptimalityForPoints>(points,f,c,s);
	}
	
	shared_ptr<IOptimalityFunction> SumWeightedSquareDiff(shared_ptr<FitPoints> points, shared_ptr<IParamFunc> f){
		OptimalityForPoints::Coefficient c=[points](ParamSet&,IParamFunc&){
			double z=0;
			for(FitPoints::Point p:*points)
				z+=p.wy;
			return 1.0/z;
		};
		OptimalityForPoints::Summand s=[](FitPoints::Point&p,ParamSet&P,IParamFunc&F){
			return pow(p.y-F(static_cast<ParamSet&&>(p.X),static_cast<ParamSet&&>(P)),2)*p.wy;
		};
		return make_shared<OptimalityForPoints>(points,f,c,s);
	}
	shared_ptr<IOptimalityFunction> ChiSquare(shared_ptr<FitPoints> points, shared_ptr<IParamFunc> f){
		OptimalityForPoints::Coefficient c=[points](ParamSet&P,IParamFunc&){
			double z=points->count()-P.Count();
			if(z<=0)
				throw GeneticException("wrong conditions for calculating xi^2: there must be at least one degree of freedom");
			return 1.0/z;
		};
		OptimalityForPoints::Summand s=[](FitPoints::Point&p,ParamSet&P,IParamFunc&F){
			return pow((p.y-F(static_cast<ParamSet&&>(p.X),static_cast<ParamSet&&>(P)))/p.wy,2);
		};
		return make_shared<OptimalityForPoints>(points,f,c,s);
	}
	shared_ptr<IOptimalityFunction> ChiSquareWithXError(shared_ptr<FitPoints> points, shared_ptr<IParamFunc> f){
		OptimalityForPoints::Coefficient c=[points](ParamSet&P,IParamFunc&){
			double z=points->count()-P.Count();
			if(z<=0)
				throw GeneticException("wrong conditions for calculating xi^2: there must be at least one degree of freedom");
			return 1.0/z;
		};
		OptimalityForPoints::Summand s=[](FitPoints::Point&p,ParamSet&P,IParamFunc&F){
			double w=pow(p.wy,2);
			for(int j=0; (j<p.X.Count())&&(j<p.WX.Count());j++){
				ParamSet x1=p.X;
				ParamSet x2=p.X;
				x1.Set(j,p.X[j]+p.WX[j]);
				x2.Set(j,p.X[j]-p.WX[j]);
				w+=pow(0.5*(F(static_cast<ParamSet&&>(x1),static_cast<ParamSet&&>(P))-F(static_cast<ParamSet&&>(x2),static_cast<ParamSet&&>(P))),2);
			}
			return pow((p.y-F(static_cast<ParamSet&&>(p.X),static_cast<ParamSet&&>(P))),2)/w;
		};
		return make_shared<OptimalityForPoints>(points,f,c,s);
	}
	
	OptimalityForPointsWithFuncError::OptimalityForPointsWithFuncError(
		shared_ptr< FitPoints > p, 
		shared_ptr< IParamFunc > f, 
		shared_ptr< IParamFunc > e, 
		Coefficient c, Summand s
	){
		points=p;
		func=f;
		error=e;
		C=c;
		S=s;
	}
	OptimalityForPointsWithFuncError::~OptimalityForPointsWithFuncError(){}
	double OptimalityForPointsWithFuncError::operator()(ParamSet&&P){
		double res=0;
		for(FitPoints::Point p:*points)
			res+=S(p,P,*func,*error);
		return res*C(P,*func,*error);
	}
	shared_ptr<IOptimalityFunction> ChiSquare(shared_ptr<FitPoints> points,shared_ptr<IParamFunc> f,shared_ptr<IParamFunc> e){
		OptimalityForPointsWithFuncError::Coefficient c=[points](ParamSet&P,IParamFunc&,IParamFunc&){
			double z=points->count()-P.Count();
			if(z<=0)
				throw GeneticException("wrong conditions for calculating xi^2: there must be at least one degree of freedom");
			return 1.0/z;
		};
		OptimalityForPointsWithFuncError::Summand s=[](FitPoints::Point&p,ParamSet&P,IParamFunc&F,IParamFunc&E){
			return pow((p.y-F(static_cast<ParamSet&&>(p.X),static_cast<ParamSet&&>(P)))/(p.wy+E(static_cast<ParamSet&&>(p.X),static_cast<ParamSet&&>(P))),2);
		};
		return make_shared<OptimalityForPointsWithFuncError>(points,f,e,c,s);
	}
	shared_ptr<IOptimalityFunction> ChiSquareWithXError(shared_ptr<FitPoints> points,shared_ptr<IParamFunc> f,shared_ptr<IParamFunc> e){
		OptimalityForPointsWithFuncError::Coefficient c=[points](ParamSet&P,IParamFunc&,IParamFunc&){
			double z=points->count()-P.Count();
			if(z<=0)
				throw GeneticException("wrong conditions for calculating xi^2: there must be at least one degree of freedom");
			return 1.0/z;
		};
		OptimalityForPointsWithFuncError::Summand s=[](FitPoints::Point&p,ParamSet&P,IParamFunc&F,IParamFunc&E){
			double w=pow(p.wy+E(static_cast<ParamSet&&>(p.X),static_cast<ParamSet&&>(P)),2);
			for(int j=0; (j<p.X.Count())&&(j<p.WX.Count());j++){
				ParamSet x1=p.X;
				ParamSet x2=p.X;
				x1.Set(j,p.X[j]+p.WX[j]);
				x2.Set(j,p.X[j]-p.WX[j]);
				w+=pow(0.5*(F(static_cast<ParamSet&&>(x1),static_cast<ParamSet&&>(P))-F(static_cast<ParamSet&&>(x2),static_cast<ParamSet&&>(P))),2);
			}
			return pow((p.y-F(static_cast<ParamSet&&>(p.X),static_cast<ParamSet&&>(P))),2)/w;
		};
		return make_shared<OptimalityForPointsWithFuncError>(points,f,e,c,s);
	}
	Parabolic::Parabolic(){}
	Parabolic::~Parabolic(){}
	double Parabolic::GetParamParabolicError(double delta, int i){
		if(delta<=0)
			throw new GeneticException("Error in parabolic error calculation: delta cannot be zero or negative");
		double s=Optimality();
		ParamSet ab=Parameters();
		ParamSet be=ab;
		ab.Set(i,ab[i]+delta);
		be.Set(i,be[i]-delta);
		double sa=OptimalityCalculator()->operator()(static_cast<ParamSet&&>(ab));
		double sb=OptimalityCalculator()->operator()(static_cast<ParamSet&&>(be));
		double da=(sa-s)/delta;
		double db=(s-sb)/delta;
		double dd=(da-db)/delta;
		if(dd<=0)
			return INFINITY;
		else
			return sqrt(2.0/dd);
	}
	ParamSet Parabolic::GetParamParabolicErrors(ParamSet&&delta){
		ParamSet res;
		for(int i=0,n=AbstractGenetic::ParamCount();i<n;i++)
			res<<GetParamParabolicError(delta[i],i);
		return res;
	}
}

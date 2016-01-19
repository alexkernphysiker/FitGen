// this file is distributed under 
// MIT license
#include <math.h>
#include <gnuplot_wrap.h>
#include <math_h/error.h>
#include <Genetic/fit.h>
using namespace std;
namespace Genetic{
	ParameterFunction::ParameterFunction(function<double(const ParamSet&,const ParamSet&)> f){func=f;}
	ParameterFunction::~ParameterFunction(){}
	double ParameterFunction::operator()(const ParamSet&X,const ParamSet&P)const{return func(X,P);}
	
	FitPoints::Point::Point(const ParamSet&x,const ParamSet&wx, double y_, double wy_){
		if(x.size()==0)
			throw Error<Point>("Empty parameter set");
		if(x.size()!=wx.size())
			throw Error<Point>("Point parameter weights size mismatch");
		__X=x;__WX=wx;
		__y=y_;__wy=wy_;
	}
	FitPoints::Point::Point(const ParamSet&x,double y_, double wy_):Point(x,parEq(x.size(),0),y_,wy_){}
	FitPoints::Point::Point(const ParamSet&x,double y_):Point(x,y_,1.0){}
	FitPoints::Point::Point(ParamSet&&x,ParamSet&&wx,double y_,double wy_):Point(x,wx,y_,wy_){}
	FitPoints::Point::Point(ParamSet&&x,double y_, double wy_):Point(x,y_,wy_){}
	FitPoints::Point::Point(ParamSet&&x,double y_):Point(x,y_){}
	FitPoints::Point::Point(const Point& src):Point(src.__X,src.__WX,src.__y,src.__wy){}
	ParamSet&FitPoints::Point::X()const{return const_cast<ParamSet&>(__X);}
	ParamSet&FitPoints::Point::WX()const{return const_cast<ParamSet&>(__WX);}
	double FitPoints::Point::y() const{return __y;}
	double FitPoints::Point::wy() const{return __wy;}
	double&FitPoints::Point::y_modify(){return __y;}
	double&FitPoints::Point::wy_modify(){return __wy;}
	
	FitPoints::FitPoints(){}
	FitPoints::~FitPoints(){}
	size_t FitPoints::size()const{return m_data.size();}
	size_t FitPoints::dimensions() const{
		if(size()==0)
			throw Error<FitPoints>("No dimensions in empty FitPoints");
		return m_data[0].X().size();
	}
	ParamSet&FitPoints::min()const{return const_cast<ParamSet&>(m_min);}
	ParamSet&FitPoints::max()const{return const_cast<ParamSet&>(m_max);}
	double FitPoints::Ymax() const{return ymax;}
	double FitPoints::Ymin() const{return ymin;}
	FitPoints& FitPoints::operator<<(Point&&point){
		if(size()>0){
			if(point.X().size()!=dimensions())
				throw Error<FitPoints>("This point has different dimensions");
			for(size_t i=0;i<dimensions();i++){
				if(point.X()[i]<m_min[i])m_min[i]=point.X()[i];
				if(point.X()[i]>m_max[i])m_max[i]=point.X()[i];
				if(point.y()<ymin)ymin=point.y();
				if(point.y()>ymax)ymax=point.y();
			}
		}else{
			m_min=m_max=point.X();
			ymin=ymax=point.y();
		}
		m_data.push_back(point);
		return *this;
	}
	shared_ptr<FitPoints> operator<<(shared_ptr<FitPoints>src,Point&&p){
		src->operator<<(static_cast<FitPoints::Point&&>(p));
		return src;
	}
	shared_ptr<FitPoints> operator<<(shared_ptr<FitPoints> src, pair<double,double>&&p){
		return src<<FitPoints::Point({p.first},p.second);
	}
	shared_ptr<FitPoints> operator<<(shared_ptr<FitPoints> src, shared_ptr<FitPoints> data){
		for(Point&P:(*data))src<<static_cast<Point&&>(P);
		return src;
	}

	FitPoints::Point&&FitPoints::operator[](size_t i)const{
		if(i>=size())
			throw Error<FitPoints>("Range check error when getting an element from FitPoints");
		return const_cast<Point&&>(m_data[i]);
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
		for(FitPoints::Point&p:(*src))
			if(condition->operator()(static_cast<ParamSet&&>(p.X())))
				res<<static_cast<FitPoints::Point&&>(p);
		return res;
	}
	shared_ptr<FitPoints> SelectFitPoints(shared_ptr<FitPoints> src, function<bool(double)> Ycond){
		auto res=make_shared<FitPoints>();
		for(FitPoints::Point&p:(*src))
			if(Ycond(p.y()))
				res<<static_cast<FitPoints::Point&&>(p);
		return res;
	}
	shared_ptr<FitPoints> SelectFitPoints(shared_ptr<FitPoints> src, shared_ptr<IParamCheck> condition,function<bool(double)> Ycond){
		auto res=make_shared<FitPoints>();
		for(FitPoints::Point&p:(*src))
			if(condition->operator()(static_cast<ParamSet&&>(p.X()))&&Ycond(p.y()))
				res<<static_cast<FitPoints::Point&&>(p);
		return res;
	}
	Distribution1D::Distribution1D(double min, double max, int bins):FitPoints(){
		if((max<min)||(bins<1))
			throw Error<Distribution1D>("Wrong constructor parameters for Distribution1D");
		double binwidth=(max-min)/double(bins);
		double halfwidth=binwidth/2.0;
		for(double x=min+halfwidth;x<max;x+=binwidth)
			operator<<(Point({x},{halfwidth},0.,1.));
	}
	void Distribution1D::Fill(double x){
		for(Point&p:(*this))
			if((x>=(p.X()[0]-p.WX()[0]))&&(x<(p.X()[0]+p.WX()[0]))){
				p.y_modify()+=1.0;
				p.wy_modify()=sqrt(p.y());
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
	shared_ptr<FitPoints> OptimalityForPoints::Points()const{return points;}
	double OptimalityForPoints::operator()(const ParamSet&P)const{
		double res=0;
		for(FitPoints::Point p:*points)
			res+=S(p,P,*func);
		return res*C(P,*func);
	}
	shared_ptr<OptimalityForPoints> SumSquareDiff(shared_ptr<FitPoints> points, shared_ptr<IParamFunc> f){
		OptimalityForPoints::Coefficient c=[](const ParamSet&,const IParamFunc&){
			return 1.0;
		};
		OptimalityForPoints::Summand s=[](const FitPoints::Point&p,const ParamSet&P,const IParamFunc&F){
			return pow(p.y()-F(p.X(),P),2);
		};
		return make_shared<OptimalityForPoints>(points,f,c,s);
	}
	
	shared_ptr<OptimalityForPoints> SumWeightedSquareDiff(shared_ptr<FitPoints> points, shared_ptr<IParamFunc> f){
		OptimalityForPoints::Coefficient c=[points](const ParamSet&,const IParamFunc&){
			double z=0;
			for(FitPoints::Point p:*points)
				z+=p.wy();
			return 1.0/z;
		};
		OptimalityForPoints::Summand s=[](const FitPoints::Point&p,const ParamSet&P,const IParamFunc&F){
			return pow(p.y()-F(p.X(),P),2)*p.wy();
		};
		return make_shared<OptimalityForPoints>(points,f,c,s);
	}
	shared_ptr<OptimalityForPoints> ChiSquare(shared_ptr<FitPoints> points, shared_ptr<IParamFunc> f){
		OptimalityForPoints::Coefficient c=[points](const ParamSet&P,const IParamFunc&){
			double z=points->size()-P.size();
			if(z<=0)
				throw Error<IOptimalityFunction>("wrong conditions for calculating xi^2: there must be at least one degree of freedom");
			return 1.0/z;
		};
		OptimalityForPoints::Summand s=[](const FitPoints::Point&p,const ParamSet&P,const IParamFunc&F){
			return pow((p.y()-F(p.X(),P))/p.wy(),2);
		};
		return make_shared<OptimalityForPoints>(points,f,c,s);
	}
	shared_ptr<OptimalityForPoints> ChiSquareWithXError(shared_ptr<FitPoints> points, shared_ptr<IParamFunc> f){
		OptimalityForPoints::Coefficient c=[points](const ParamSet&P,const IParamFunc&){
			double z=points->size()-P.size();
			if(z<=0)
				throw Error<IOptimalityFunction>("wrong conditions for calculating xi^2: there must be at least one degree of freedom");
			return 1.0/z;
		};
		OptimalityForPoints::Summand s=[](const FitPoints::Point&p,const ParamSet&P,const IParamFunc&F){
			double w=pow(p.wy(),2);
			for(size_t j=0; j<p.X().size();j++){
				ParamSet x1=p.X();
				ParamSet x2=p.X();
				x1[j]+=p.WX()[j];
				x2[j]-=p.WX()[j];
				w+=pow(0.5*(F(x1,P)-F(x2,P)),2);
			}
			return pow((p.y()-F(p.X(),P)),2)/w;
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
	shared_ptr< FitPoints > OptimalityForPointsWithFuncError::Points()const{return points;}
	double OptimalityForPointsWithFuncError::operator()(const ParamSet&P)const{
		double res=0;
		for(FitPoints::Point p:*points)
			res+=S(p,P,*func,*error);
		return res*C(P,*func,*error);
	}
	shared_ptr<OptimalityForPointsWithFuncError> ChiSquare(shared_ptr<FitPoints> points,shared_ptr<IParamFunc> f,shared_ptr<IParamFunc> e){
		OptimalityForPointsWithFuncError::Coefficient c=[points](const ParamSet&P,const IParamFunc&,const IParamFunc&){
			double z=points->size()-P.size();
			if(z<=0)
				throw Error<IOptimalityFunction>("wrong conditions for calculating xi^2: there must be at least one degree of freedom");
			return 1.0/z;
		};
		OptimalityForPointsWithFuncError::Summand s=[](const FitPoints::Point&p,const ParamSet&P,const IParamFunc&F,const IParamFunc&E){
			return pow((p.y()-F(p.X(),P))/(p.wy()+E(p.X(),P)),2);
		};
		return make_shared<OptimalityForPointsWithFuncError>(points,f,e,c,s);
	}
	shared_ptr<OptimalityForPointsWithFuncError> ChiSquareWithXError(shared_ptr<FitPoints> points,shared_ptr<IParamFunc> f,shared_ptr<IParamFunc> e){
		OptimalityForPointsWithFuncError::Coefficient c=[points](const ParamSet&P,const IParamFunc&,const IParamFunc&){
			double z=points->size()-P.size();
			if(z<=0)
				throw Error<IOptimalityFunction>("wrong conditions for calculating xi^2: there must be at least one degree of freedom");
			return 1.0/z;
		};
		OptimalityForPointsWithFuncError::Summand s=[](const FitPoints::Point&p,const ParamSet&P,const IParamFunc&F,const IParamFunc&E){
			double w=pow(p.wy()+E(p.X(),P),2);
			for(size_t j=0; j<p.X().size();j++){
				ParamSet x1=p.X();
				ParamSet x2=p.X();
				x1[j]+=p.WX()[j];
				x2[j]-=p.WX()[j];
				w+=pow(0.5*(F(x1,P)-F(x2,P)),2);
			}
			return pow((p.y()-F(p.X(),P)),2)/w;
		};
		return make_shared<OptimalityForPointsWithFuncError>(points,f,e,c,s);
	}
	Parabolic::Parabolic(){}
	Parabolic::~Parabolic(){}
	double Parabolic::GetParamParabolicError(double delta, int i)const{
		if(delta<=0)
			throw new Error<Parabolic>("Error in parabolic error calculation: delta cannot be zero or negative");
		double s=Optimality();
		ParamSet ab=Parameters();
		ParamSet be=ab;
		ab[i]+=delta;
		be[i]-=delta;
		double sa=OptimalityCalculator()->operator()(ab);
		double sb=OptimalityCalculator()->operator()(be);
		double dd=(sa-2.0*s+sb)/pow(delta,2);
		if(dd<=0)
			return INFINITY;
		else
			return sqrt(2.0/dd);
	}
	ParamSet Parabolic::GetParamParabolicErrors(ParamSet&&delta)const{
		ParamSet res;
		for(int i=0,n=AbstractGenetic::ParamCount();i<n;i++)
			res<<GetParamParabolicError(delta[i],i);
		return res;
	}
	PlotPoints1D::PlotPoints1D():Plot<double>(){}
	PlotPoints1D::~PlotPoints1D(){}
	PlotPoints1D&PlotPoints1D::Points(string&&name,shared_ptr<FitPoints> points,size_t param_index){
		OutputPlot(static_cast<std::string&&>(name),[points,param_index](ofstream&out){
			for(auto p:*points)
				out<<p.X()[param_index]<<" "<<p.y()<<" "<<p.WX()[param_index]<<" "<<p.wy()<<"\n";
		},"using 1:2:($1-$3):($1+$3):($2-$4):($2+$4) with xyerrorbars");
		return *this;
	}
	PlotPoints1D&PlotPoints1D::PointsWithoutErrors(string&&name,shared_ptr<FitPoints> points, size_t param_index){
		OutputPlot(static_cast<std::string&&>(name),[points,param_index](ofstream&out){
			for(auto p:*points)
				out<<p.X()[param_index]<<" "<<p.y()<<"\n";
		},"using 1:2");
		return *this;
	}
	
}

// this file is distributed under 
// MIT license
#include <math.h>
#include <climits>
#include <gnuplot_wrap.h>
#include <math_h/error.h>
#include <Genetic/fit.h>
namespace Genetic{
	using namespace std;
	using namespace MathTemplates;
	using namespace GnuplotWrap;
	ParameterFunction::ParameterFunction(const function<double(const ParamSet&,const ParamSet&)> f){func=f;}
	ParameterFunction::~ParameterFunction(){}
	double ParameterFunction::operator()(const ParamSet&X,const ParamSet&P)const{return func(X,P);}
	
	FitPoints::Point::Point(const ParamSet&x,const ParamSet&wx, const double y_, const double wy_){
		if(x.size()==0)
			throw Exception<Point>("Empty parameter set");
		if(x.size()!=wx.size())
			throw Exception<Point>("Point parameter weights size mismatch");
		__X=x;__WX=wx;
		__y=y_;__wy=wy_;
	}
	FitPoints::Point::Point(const ParamSet&x,const double y_, const double wy_):Point(x,parEq(x.size(),0),y_,wy_){}
	FitPoints::Point::Point(const ParamSet&x,const double y_):Point(x,y_,1.0){}
	FitPoints::Point::Point(const ParamSet&&x,const ParamSet&&wx,const double y_,const double wy_):Point(x,wx,y_,wy_){}
	FitPoints::Point::Point(const ParamSet&&x,const double y_, const double wy_):Point(x,y_,wy_){}
	FitPoints::Point::Point(const ParamSet&&x,const double y_):Point(x,y_){}
	FitPoints::Point::Point(const Point& src):Point(src.__X,src.__WX,src.__y,src.__wy){}
	const ParamSet&FitPoints::Point::X()const{return __X;}
	const ParamSet&FitPoints::Point::WX()const{return __WX;}
	double FitPoints::Point::y() const{return __y;}
	double FitPoints::Point::wy() const{return __wy;}
	double&FitPoints::Point::var_y(){return __y;}
	double&FitPoints::Point::var_wy(){return __wy;}
	
	FitPoints::FitPoints(){}
	FitPoints::FitPoints(const SortedPoints<double>& h){
		for(const point<double>&p:h)
			operator<<(Point({p.X()},{0.0},p.Y(),0.0));
	}
	FitPoints::FitPoints(const BiSortedPoints<double>&d){
		d.FullCycle([this](const point3d<double>&p){
			operator<<(Point({p.X(),p.Y()},{0.0,0.0},p.Z(),0.0));
		});
	}
	FitPoints::FitPoints(const SortedPoints<value<double>>& h){
		for(const point<value<double>>&p:h)
			operator<<(Point({p.X().val()},{p.X().delta()},p.Y().val(),p.Y().delta()));
	}
	FitPoints::FitPoints(const BiSortedPoints<value<double>>&d){
		d.FullCycle([this](const point3d<value<double>>&p){
			operator<<(Point({p.X().val(),p.Y().val()},{p.X().delta(),p.Y().delta()},p.Z().val(),p.Z().delta()));
		});
	}
	FitPoints::~FitPoints(){}
	size_t FitPoints::size()const{return m_data.size();}
	size_t FitPoints::dimensions() const{
		if(size()==0)
			throw Exception<FitPoints>("No dimensions in empty FitPoints");
		return m_data[0].X().size();
	}
	const ParamSet&FitPoints::min()const{return m_min;}
	const ParamSet&FitPoints::max()const{return m_max;}
	double FitPoints::Ymax() const{return ymax;}
	double FitPoints::Ymin() const{return ymin;}
	FitPoints& FitPoints::operator<<(const Point&point){
		if(size()>0){
			if(point.X().size()!=dimensions())
				throw Exception<FitPoints>("This point has different dimensions");
			for(size_t i=0;i<dimensions();i++){
				if(point.X()[i]<m_min[i])m_min(i)=point.X()[i];
				if(point.X()[i]>m_max[i])m_max(i)=point.X()[i];
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
	const SortedPoints<value<double>> FitPoints::Hist1(const size_t parameter_index) const{
		SortedPoints<value<double>> data;
		for(const Point&P:m_data)
			data<<point<value<double>>(
				value<double>(P.X()[parameter_index],P.WX()[parameter_index]),
						   value<double>(P.y(),P.wy())
			);
		return data;
	}
	const SortedPoints<value<double>> FitPoints::Hist1(const size_t parameter_index_x,const size_t parameter_index_y) const{
		SortedPoints<value<double>> data;
		for(const Point&P:m_data)
			data<<point<value<double>>(
				value<double>(P.X()[parameter_index_x],P.WX()[parameter_index_x]),
						   value<double>(P.X()[parameter_index_y],P.WX()[parameter_index_y])
			);
		return data;
	}
	const SortedPoints<double> FitPoints::Line(const size_t parameter_index)const{
		SortedPoints<double> data;
		for(const Point&P:m_data)
			data<<point<double>(P.X()[parameter_index],P.y());
		return data;
	}
	const SortedPoints<double> FitPoints::Line(const size_t parameter_index_x,const size_t parameter_index_y) const{
		SortedPoints<double> data;
		for(const Point&P:m_data)
			data<<point<double>(P.X()[parameter_index_x],P.X()[parameter_index_y]);
		return data;
	}
	

	shared_ptr<FitPoints> operator<<(shared_ptr<FitPoints>src,const Point&p){
		src->operator<<(p);
		return src;
	}
	shared_ptr<FitPoints> operator<<(shared_ptr<FitPoints>src,const Point&&p){
		src->operator<<(p);
		return src;
	}
	shared_ptr< FitPoints > operator<<(shared_ptr< FitPoints > src, const point<double>& p){
		return src<<FitPoints::Point({p.X()},p.Y());
	}
	shared_ptr<FitPoints> operator<<(shared_ptr<FitPoints> src,const  point<double>&&p){
		return src<<FitPoints::Point({p.X()},p.Y());
	}
	shared_ptr<FitPoints> operator<<(shared_ptr<FitPoints> src,const  shared_ptr<FitPoints> data){
		for(const Point&P:(*data))src<<P;
		return src;
	}

	const FitPoints::Point&FitPoints::operator[](const size_t i)const{
		if(i>=size())
			throw Exception<FitPoints>("Range check error when getting an element from FitPoints");
		return m_data[i];
	}
	FitPoints::const_iterator FitPoints::begin()const {
		return m_data.begin();
	}
	FitPoints::const_iterator FitPoints::cbegin() const{
		return m_data.cbegin();
	}
	FitPoints::const_iterator FitPoints::end()const {
		return m_data.end();
	}
	FitPoints::const_iterator FitPoints::cend() const{
		return m_data.cend();
	}
	shared_ptr<FitPoints> SelectFitPoints(shared_ptr<FitPoints> src,const shared_ptr<IParamCheck> condition){
		auto res=make_shared<FitPoints>();
		for(const FitPoints::Point&p:(*src))
			if(condition->operator()(p.X()))res<<p;
		return res;
	}
	shared_ptr<FitPoints> SelectFitPoints(shared_ptr<FitPoints> src, const function<bool(double)> Ycond){
		auto res=make_shared<FitPoints>();
		for(const FitPoints::Point&p:(*src))
			if(Ycond(p.y()))
				res<<p;
		return res;
	}
	shared_ptr<FitPoints> SelectFitPoints(shared_ptr<FitPoints> src, const shared_ptr<IParamCheck> condition,const function<bool(double)> Ycond){
		auto res=make_shared<FitPoints>();
		for(const FitPoints::Point&p:(*src))
			if(condition->operator()(p.X())&&Ycond(p.y()))
				res<<p;
		return res;
	}
	OptimalityForPoints::OptimalityForPoints(
		const std::shared_ptr< FitPoints > p, 
		const shared_ptr<IParamFunc> f,
		const OptimalityForPoints::Coefficient c, 
		const OptimalityForPoints::Summand s
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
	shared_ptr<OptimalityForPoints> SumSquareDiff(const shared_ptr<FitPoints> points, const shared_ptr<IParamFunc> f){
		OptimalityForPoints::Coefficient c=[](const ParamSet&,const IParamFunc&){
			return 1.0;
		};
		OptimalityForPoints::Summand s=[](const FitPoints::Point&p,const ParamSet&P,const IParamFunc&F){
			return pow(p.y()-F(p.X(),P),2);
		};
		return make_shared<OptimalityForPoints>(points,f,c,s);
	}
	
	shared_ptr<OptimalityForPoints> SumWeightedSquareDiff(const shared_ptr<FitPoints> points, const shared_ptr<IParamFunc> f){
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
	shared_ptr<OptimalityForPoints> ChiSquare(const shared_ptr<FitPoints> points, const shared_ptr<IParamFunc> f){
		OptimalityForPoints::Coefficient c=[points](const ParamSet&P,const IParamFunc&){
			double z=points->size()-P.size();
			if(z<=0)
				throw Exception<IOptimalityFunction>("wrong conditions for calculating xi^2: there must be at least one degree of freedom");
			return 1.0/z;
		};
		OptimalityForPoints::Summand s=[](const FitPoints::Point&p,const ParamSet&P,const IParamFunc&F){
			return pow((p.y()-F(p.X(),P))/p.wy(),2);
		};
		return make_shared<OptimalityForPoints>(points,f,c,s);
	}
	shared_ptr<OptimalityForPoints> ChiSquareWithXError(const shared_ptr<FitPoints> points, const shared_ptr<IParamFunc> f){
		OptimalityForPoints::Coefficient c=[points](const ParamSet&P,const IParamFunc&){
			double z=points->size()-P.size();
			if(z<=0)
				throw Exception<IOptimalityFunction>("wrong conditions for calculating xi^2: there must be at least one degree of freedom");
			return 1.0/z;
		};
		OptimalityForPoints::Summand s=[](const FitPoints::Point&p,const ParamSet&P,const IParamFunc&F){
			double w=pow(p.wy(),2);
			for(size_t j=0; j<p.X().size();j++){
				ParamSet x1=p.X();
				ParamSet x2=p.X();
				x1(j)+=p.WX()[j];
				x2(j)-=p.WX()[j];
				w+=pow(0.5*(F(x1,P)-F(x2,P)),2);
			}
			return pow((p.y()-F(p.X(),P)),2)/w;
		};
		return make_shared<OptimalityForPoints>(points,f,c,s);
	}
	
	OptimalityForPointsWithFuncError::OptimalityForPointsWithFuncError(const shared_ptr< FitPoints > p, const Func f,const Coefficient c, const Summand s){
		points=p;
		m_func=f;
		C=c;
		S=s;
	}
	OptimalityForPointsWithFuncError::~OptimalityForPointsWithFuncError(){}
	shared_ptr< FitPoints > OptimalityForPointsWithFuncError::Points()const{return points;}
	double OptimalityForPointsWithFuncError::operator()(const ParamSet&P)const{
		double res=0;
		for(FitPoints::Point p:*points)
			res+=S(p,P,m_func);
		return res*C(P,m_func);
	}
	shared_ptr<OptimalityForPointsWithFuncError> ChiSquare(const shared_ptr<FitPoints> points,const OptimalityForPointsWithFuncError::Func f){
		OptimalityForPointsWithFuncError::Coefficient c=[points](const ParamSet&P,const OptimalityForPointsWithFuncError::Func&){
			double z=points->size()-P.size();
			if(z<=0)
				throw Exception<IOptimalityFunction>("wrong conditions for calculating xi^2: there must be at least one degree of freedom");
			return 1.0/z;
		};
		OptimalityForPointsWithFuncError::Summand s=[](const FitPoints::Point&p,const ParamSet&P,const OptimalityForPointsWithFuncError::Func&F){
			return pow((p.y()-F(p.X(),P).val())/(p.wy()+F(p.X(),P).delta()),2);
		};
		return make_shared<OptimalityForPointsWithFuncError>(points,f,c,s);
	}
	shared_ptr<OptimalityForPointsWithFuncError> ChiSquareWithXError(const shared_ptr<FitPoints> points,const OptimalityForPointsWithFuncError::Func f){
		OptimalityForPointsWithFuncError::Coefficient c=[points](const ParamSet&P,const OptimalityForPointsWithFuncError::Func&){
			double z=points->size()-P.size();
			if(z<=0)
				throw Exception<IOptimalityFunction>("wrong conditions for calculating xi^2: there must be at least one degree of freedom");
			return 1.0/z;
		};
		OptimalityForPointsWithFuncError::Summand s=[](const FitPoints::Point&p,const ParamSet&P,const OptimalityForPointsWithFuncError::Func&F){
			double w=pow(p.wy()+F(p.X(),P).delta(),2);
			for(size_t j=0; j<p.X().size();j++){
				ParamSet x1=p.X();
				ParamSet x2=p.X();
				x1(j)+=p.WX()[j];
				x2(j)-=p.WX()[j];
				w+=pow(0.5*(F(x1,P).val()-F(x2,P).val()),2);
			}
			return pow((p.y()-F(p.X(),P).val()),2)/w;
		};
		return make_shared<OptimalityForPointsWithFuncError>(points,f,c,s);
	}
	Parabolic::Parabolic(){
		m_uncertainty_cache=make_shared<vector<value<double>>>();
		m_iter_number=make_shared<unsigned long long int>(ULLONG_MAX);
	}
	Parabolic::~Parabolic(){}
	double Parabolic::GetParamParabolicError(const double delta, const size_t i)const{
		if(delta<=0)
			throw new Exception<Parabolic>("Exception in parabolic error calculation: delta cannot be zero or negative");
		double s=Optimality();
		ParamSet ab=Parameters();
		ParamSet be=ab;
		ab(i)+=delta;
		be(i)-=delta;
		double sa=OptimalityCalculator()->operator()(ab);
		double sb=OptimalityCalculator()->operator()(be);
		double dd=(sa-2.0*s+sb)/pow(delta,2);
		if(dd<=0)
			return INFINITY;
		else
			return sqrt(2.0/dd);
	}
	Parabolic& Parabolic::SetUncertaintyCalcDeltas(const ParamSet& P){
		m_delta=P;
		(*m_iter_number)=ULLONG_MAX;
		return *this;
	}
	Parabolic& Parabolic::SetUncertaintyCalcDeltas(const ParamSet&& P){
		return SetUncertaintyCalcDeltas(P);
	}
	const vector<value<double>>&Parabolic::ParametersWithUncertainties()const{
		if(iteration_count()!=(*m_iter_number)){
			m_uncertainty_cache->clear();
			for(size_t i=0;(i<ParamCount())&&(i<m_delta.size());i++){
				m_uncertainty_cache->push_back(value<double>(Parameters()[i],GetParamParabolicError(m_delta[i],i)));
			}
			(*m_iter_number)=iteration_count();
		}
		return *m_uncertainty_cache;
	}
}

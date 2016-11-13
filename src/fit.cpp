// this file is distributed under 
// MIT license
#include <math.h>
#include <climits>
#include <math_h/error.h>
#include <Genetic/fit.h>
namespace Genetic{
	using namespace std;
	using namespace MathTemplates;
	ParameterFunction::ParameterFunction(const function<double(const ParamSet&,const ParamSet&)> f){func=f;}
	ParameterFunction::~ParameterFunction(){}
	double ParameterFunction::operator()(const ParamSet&X,const ParamSet&P)const{return func(X,P);}
	
	FitPoints::Point::~Point(){}
	FitPoints::Point::Point(const FitPoints::Point&src):yy(src.yy){
		for(const auto&x:src.XX)XX.push_back(x);
	}
	FitPoints::Point::Point(const vector<value<double>>&x,const value<double>&y_):yy(y_){
		for(const auto&p:x)XX.push_back(p);
		if(XX.size()==0)
			throw Exception<FitPoints::Point>("Cannot create Point with zero arguments");
	}
	FitPoints::Point::Point(const vector<value<double>>&&x,const value<double>&y_):Point(x,y_){}
	FitPoints::Point::Point(const vector<value<double>>&x,const value<double>&&y_):Point(x,y_){}
	FitPoints::Point::Point(const vector<value<double>>&&x,const value<double>&&y_):Point(x,y_){}
	const vector<value<double>>&FitPoints::Point::X() const{return XX;}
	const ParamSet FitPoints::Point::x() const{
		ParamSet res;
		for(const auto&x:X())
			res<<x.val();
		return res;
	}
	const value<double>&FitPoints::Point::y() const{return yy;}
	value<double>&FitPoints::Point::var_y(){return yy;}
	
	
	FitPoints::FitPoints(){}
	FitPoints::FitPoints(const SortedPoints<double>& h){
		for(const point<double>&p:h)
			operator<<(Point({p.X()},p.Y()));
	}
	FitPoints::FitPoints(const BiSortedPoints<double>&d){
		d.FullCycle([this](const point3d<double>&p){
			operator<<(Point({p.X(),p.Y()},p.Z()));
		});
	}
	FitPoints::FitPoints(const SortedPoints<value<double>>& h){
		for(const point<value<double>>&p:h)
			operator<<(Point({p.X()},p.Y()));
	}
	FitPoints::FitPoints(const BiSortedPoints<value<double>>&d){
		d.FullCycle([this](const point3d<value<double>>&p){
			operator<<(Point({p.X(),p.Y()},p.Z()));
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
	const double&FitPoints::Ymax() const{return ymax;}
	const double&FitPoints::Ymin() const{return ymin;}
	FitPoints& FitPoints::operator<<(const Point&point){
		if(size()>0){
			if(point.X().size()!=dimensions())
				throw Exception<FitPoints>("This point has different dimensions");
			for(size_t i=0;i<dimensions();i++){
				if(point.X()[i].val()<m_min[i])m_min(i)=point.X()[i].val();
				if(point.X()[i].val()>m_max[i])m_max(i)=point.X()[i].val();
				if(point.y().val()<ymin)ymin=point.y().val();
				if(point.y().val()>ymax)ymax=point.y().val();
			}
		}else{
			m_min=m_max={};
			for(const auto&x:point.X()){
				m_min << x.val();
				m_max << x.val();
			}
			ymin=ymax=point.y().val();
		}
		m_data.push_back(point);
		return *this;
	}
	const SortedPoints<value<double>> FitPoints::Hist1(const size_t parameter_index) const{
		SortedPoints<value<double>> data;
		for(const Point&P:m_data)
			data<<point<value<double>>(P.X()[parameter_index],P.y());
		return data;
	}
	const SortedPoints<value<double>> FitPoints::Hist1(const size_t parameter_index_x,const size_t parameter_index_y) const{
		SortedPoints<value<double>> data;
		for(const Point&P:m_data)
			data<<point<value<double>>(P.X()[parameter_index_x],P.X()[parameter_index_y]);
		return data;
	}
	const SortedPoints<double> FitPoints::Line(const size_t parameter_index)const{
		SortedPoints<double> data;
		for(const Point&P:m_data)
			data<<point<double>(P.X()[parameter_index].val(),P.y().val());
		return data;
	}
	const SortedPoints<double> FitPoints::Line(const size_t parameter_index_x,const size_t parameter_index_y) const{
		SortedPoints<double> data;
		for(const Point&P:m_data)
			data<<point<double>(P.X()[parameter_index_x].val(),P.X()[parameter_index_y].val());
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
	FitPoints::const_iterator FitPoints::end()const {
		return m_data.end();
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
		OptimalityForPoints::Coefficient c=[](const ParamSet&,const IParamFunc&){return 1.0;};
		OptimalityForPoints::Summand s=[](const FitPoints::Point&p,const ParamSet&P,const IParamFunc&F){
			return pow(p.y().val()-F(p.x(),P),2);
		};
		return make_shared<OptimalityForPoints>(points,f,c,s);
	}
	
	shared_ptr<OptimalityForPoints> ChiSquare(const shared_ptr<FitPoints> points, const shared_ptr<IParamFunc> f){
		OptimalityForPoints::Coefficient c=[points](const ParamSet&,const IParamFunc&){return 1.0;};
		OptimalityForPoints::Summand s=[](const FitPoints::Point&p,const ParamSet&P,const IParamFunc&F){
			return p.y().NumCompare(F(p.x(),P));
		};
		return make_shared<OptimalityForPoints>(points,f,c,s);
	}
	shared_ptr<OptimalityForPoints> ChiSquareWithXError(const shared_ptr<FitPoints> points, const shared_ptr<IParamFunc> f){
		OptimalityForPoints::Coefficient c=[points](const ParamSet&,const IParamFunc&){return 1.0;};
		OptimalityForPoints::Summand s=[](const FitPoints::Point&p,const ParamSet&P,const IParamFunc&F){
			double w=pow(p.y().uncertainty(),2);
			for(size_t j=0; j<p.X().size();j++){
				ParamSet x1=p.x(),x2=x1;
				x1(j)+=p.X()[j].uncertainty();
				x2(j)-=p.X()[j].uncertainty();
				w+=pow(0.5*(F(x1,P)-F(x2,P)),2);
			}
			return pow((p.y().val()-F(p.x(),P)),2)/w;
		};
		return make_shared<OptimalityForPoints>(points,f,c,s);
	}
	
}

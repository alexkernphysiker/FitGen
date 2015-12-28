// this file is distributed under 
// MIT license
#ifndef ____lrEPWamH___
#define ____lrEPWamH___
#include <functional>
#include "../gnuplot_wrap.h"
#include "abstract.h"
#include "genetic.h"
#include "paramfunc.h"
namespace Genetic{
	using namespace std;
	class ParameterFunction:public IParamFunc{
	public:
		ParameterFunction(paramFunc f);
		virtual ~ParameterFunction();
		virtual double operator()(const ParamSet&X,const ParamSet&P)const override;
	private:
		paramFunc func;
	};
	class FitPoints{
	public:
		class Point{
		public:
			Point(const Point &src);
			Point(const ParamSet&x,double y_);
			Point(const ParamSet&x,double y_,double wy_);
			Point(const ParamSet&x,const ParamSet&wx,double y_,double wy_);
			Point(ParamSet&&x,double y_);
			Point(ParamSet&&x,double y_,double wy_);
			Point(ParamSet&&x,ParamSet&&wx,double y_,double wy_);
			ParamSet&X()const;
			ParamSet&WX()const;
			double y()const;
			double wy()const;
			double&Y();
			double&WY();
		private:
			ParamSet __X;
			ParamSet __WX;
			double __y;
			double __wy;
		};
		FitPoints();
		virtual ~FitPoints();
		FitPoints&operator<<(Point&&point);
		Point&&operator[](size_t i)const;
		size_t size()const;
		size_t dimensions()const;
		ParamSet&min()const;
		ParamSet&max()const;
		typedef vector<Point>::iterator iterator;
		typedef vector<Point>::const_iterator const_iterator;
		iterator begin();
		const_iterator cbegin()const;
		iterator end();
		const_iterator cend() const;
	private:
		vector<Point> m_data;
		ParamSet m_min,m_max;
	};
	typedef FitPoints::Point Point;
	shared_ptr<FitPoints> operator<<(shared_ptr<FitPoints> src,Point&&p);
	shared_ptr<FitPoints> operator<<(shared_ptr<FitPoints> src,pair<double,double>&&p);
	shared_ptr<FitPoints> SelectFitPoints(shared_ptr<FitPoints> src,shared_ptr<IParamCheck> condition);
	shared_ptr<FitPoints> SelectFitPoints(shared_ptr<FitPoints> src,function<bool(double)> Ycond);
	shared_ptr<FitPoints> SelectFitPoints(shared_ptr<FitPoints> src,shared_ptr<IParamCheck> condition,function<bool(double)> Ycond);
	template<class IndexerX,class IndexerY=IndexerX>
	shared_ptr<FitPoints> FitPointsXY(int from,int to,IndexerX X,IndexerY Y){
		auto res=make_shared<FitPoints>();
		for(int i=from;i<=to;i++)
			res<<FitPoints::Point({X[i]},Y[i]);
		return res;
	}
	template<class IndexerX,class IndexerY=IndexerX,class IndexerWY=IndexerY>
	shared_ptr<FitPoints> FitPointsXYdY(int from,int to,IndexerX X,IndexerY Y,IndexerWY WY){
		auto res=make_shared<FitPoints>();
		for(int i=from;i<=to;i++)
			res<<FitPoints::Point({X[i]},Y[i],WY[i]);
		return res;
	}
	template<class IndexerX,class IndexerWX=IndexerX,class IndexerY=IndexerX,class IndexerWY=IndexerY>
	shared_ptr<FitPoints> FitPointsXdXYdY(int from,int to,IndexerX X,IndexerWX WX,IndexerY Y,IndexerWY WY){
		auto res=make_shared<FitPoints>();
		for(int i=from;i<=to;i++)
			res<<FitPoints::Point({X[i]},{WX[i]},Y[i],WY[i]);
		return res;
	}
	class Distribution1D:public FitPoints{
	public:
		Distribution1D(double min, double max, int bins);
		void Fill(double x);
	};
	
	class OptimalityForPoints:public IOptimalityFunction{
	public:
		typedef function<double(const ParamSet&,const IParamFunc&)> Coefficient;
		typedef function<double(const FitPoints::Point&,const ParamSet&,const IParamFunc&)> Summand;
		OptimalityForPoints(shared_ptr<FitPoints> p, shared_ptr<IParamFunc> f,Coefficient c,Summand s);
		virtual ~OptimalityForPoints();
		virtual double operator()(const ParamSet&P)const override;
		shared_ptr<FitPoints> Points()const;
	protected:
		shared_ptr<FitPoints> points;
		shared_ptr<IParamFunc> func;
		Coefficient C;
		Summand S;
	};
	class Parabolic:public virtual AbstractGenetic{
	protected:
		Parabolic();
	public:
		virtual ~Parabolic();
		double GetParamParabolicError(double delta,int i)const;
		ParamSet GetParamParabolicErrors(ParamSet&&delta)const;
	};
	shared_ptr<OptimalityForPoints> SumSquareDiff(shared_ptr<FitPoints> points, shared_ptr<IParamFunc> f);
	shared_ptr<OptimalityForPoints> SumWeightedSquareDiff(shared_ptr<FitPoints> points, shared_ptr<IParamFunc> f);
	shared_ptr<OptimalityForPoints> ChiSquare(shared_ptr<FitPoints> points, shared_ptr<IParamFunc> f);
	shared_ptr<OptimalityForPoints> ChiSquareWithXError(shared_ptr<FitPoints> points, shared_ptr<IParamFunc> f);
	template<class GENETIC,shared_ptr<OptimalityForPoints> OptimalityAlgorithm(shared_ptr<FitPoints>,shared_ptr<IParamFunc>)>
	class Fit:public virtual GENETIC,public virtual Parabolic{
	private:
		shared_ptr<IParamFunc> m_func;
	protected:
		Fit(shared_ptr<IParamFunc> f):GENETIC(),Parabolic(){
			m_func=f;
		}
	public:
		Fit(
			shared_ptr<FitPoints> points, 
			shared_ptr<IParamFunc> f
		):AbstractGenetic(OptimalityAlgorithm(points,f)){
			m_func=f;
		}
		Fit(shared_ptr<FitPoints> points, paramFunc f):Fit(points,make_shared<ParameterFunction>(f)){}
		virtual ~Fit(){}
		double operator()(ParamSet&&X)const{return m_func->operator()(X,AbstractGenetic::Parameters());}
		shared_ptr<IParamFunc> Func()const{return m_func;}
		shared_ptr<FitPoints> Points()const{return dynamic_pointer_cast<OptimalityForPoints>(AbstractGenetic::OptimalityCalculator())->Points();}
	};
	template<class GENETIC,class FUNC,shared_ptr<OptimalityForPoints> OptimalityAlgorithm(shared_ptr<FitPoints>,shared_ptr<IParamFunc>)>
	class FitFunction:public virtual Fit<GENETIC,OptimalityAlgorithm>{
	public:
		typedef FUNC functype;
		FitFunction(shared_ptr<FitPoints> points):
			AbstractGenetic(OptimalityAlgorithm(points,make_shared<FUNC>())),
			Fit<GENETIC,OptimalityAlgorithm>(make_shared<FUNC>()){}
		virtual ~FitFunction(){}
	};
	
	class OptimalityForPointsWithFuncError:public IOptimalityFunction{
	public:
		typedef function<double(const ParamSet&,const IParamFunc&,const IParamFunc&)> Coefficient;
		typedef function<double(const FitPoints::Point&,const ParamSet&,const IParamFunc&,const IParamFunc&)> Summand;
		OptimalityForPointsWithFuncError(shared_ptr<FitPoints> p,shared_ptr<IParamFunc> f,shared_ptr<IParamFunc> e,Coefficient c,Summand s);
		virtual ~OptimalityForPointsWithFuncError();
		virtual double operator()(const ParamSet&P)const override;
		shared_ptr<FitPoints> Points()const;
	protected:
		shared_ptr<FitPoints> points;
		shared_ptr<IParamFunc> func;
		shared_ptr<IParamFunc> error;
		Coefficient C;
		Summand S;
	};
	shared_ptr<OptimalityForPointsWithFuncError> ChiSquare(shared_ptr<FitPoints> points, shared_ptr<IParamFunc> f,shared_ptr<IParamFunc> e);
	shared_ptr<OptimalityForPointsWithFuncError> ChiSquareWithXError(shared_ptr<FitPoints> points, shared_ptr<IParamFunc> f,shared_ptr<IParamFunc> e);
	template<class GENETIC,shared_ptr<OptimalityForPointsWithFuncError> OptimalityAlgorithm(shared_ptr<FitPoints>,shared_ptr<IParamFunc>,shared_ptr<IParamFunc>)>
	class FitFunctionWithError:public virtual GENETIC,public virtual Parabolic{
	private:
		shared_ptr<IParamFunc> m_func;
	public:
		FitFunctionWithError(shared_ptr<FitPoints> points,shared_ptr<IParamFunc> f,shared_ptr<IParamFunc> e)
			:AbstractGenetic(OptimalityAlgorithm(points,f,e)),GENETIC(),Parabolic(){
			m_func=f;
		}
		FitFunctionWithError(shared_ptr<FitPoints> points,paramFunc f,paramFunc e):
			FitFunctionWithError(points,make_shared<ParameterFunction>(f),make_shared<ParameterFunction>(e)){}
		virtual ~FitFunctionWithError(){}
		double operator()(ParamSet&&X)const{return m_func->operator()(X,AbstractGenetic::Parameters());}
		shared_ptr<FitPoints> Points()const{return dynamic_pointer_cast<OptimalityForPointsWithFuncError>(AbstractGenetic::OptimalityCalculator())->Points();}
	};
	class PlotPoints1D:public Plot<double>{
	public:
		PlotPoints1D();
		virtual ~PlotPoints1D();
		PlotPoints1D& Points(std::string&&name,shared_ptr<FitPoints> points,size_t param_index=0);
		PlotPoints1D& PointsWithoutErrors(std::string&&name,shared_ptr<FitPoints> points,size_t param_index=0);
	};
	template<class FIT>
	class PlotFit1D:public PlotPoints1D{
	private:
		double min,max;
	public:
		PlotFit1D():PlotPoints1D(){
			min=+INFINITY;
			max=-INFINITY;
		}
		virtual ~PlotFit1D(){}
		PlotFit1D& Points(std::string&&name,shared_ptr<FitPoints> points,size_t param_index=0){
			PlotPoints1D::Points(static_cast<std::string&&>(name),points,param_index);
			for(FitPoints::Point p:*points){
				if(p.X()[param_index]<min)
					min=p.X()[param_index];
				if(p.X()[param_index]>max)
					max=p.X()[param_index];
			}
			return *this;
		}
		PlotFit1D& PointsWithoutErrors(std::string&&name,shared_ptr<FitPoints> points,size_t param_index=0){
			PlotPoints1D::PointsWithoutErrors(static_cast<std::string&&>(name),points,param_index);
			for(FitPoints::Point p:*points){
				if(p.X()[param_index]<min)
					min=p.X()[param_index];
				if(p.X()[param_index]>max)
					max=p.X()[param_index];
			}
			return *this;
		}
		PlotFit1D& Fit(std::string&&name,const FIT&fit,double step){
			if(max<min)throw math_h_error<PlotFit1D>("No FitPoints instance initialized the ranges on the plot");
			Plot<double>::Line(static_cast<std::string&&>(name),[&fit](double x){return fit({x});},min,max,step);
			return *this;
		}
		PlotFit1D& Fit(std::string&&name,std::string&&name_func,const FIT&fit,double step){
			Points(static_cast<std::string&&>(name_func),fit.Points());
			Plot<double>::Line(static_cast<std::string&&>(name),[&fit](double x){return fit({x});},min,max,step);
			return *this;
		}
		PlotFit1D& ParamFunc(std::string&&name,IParamFunc&&func,const FIT&fit,double step){
			if(max<min)throw math_h_error<PlotFit1D>("No FitPoints instance initialized the ranges on the plot");
			Plot<double>::Line(static_cast<std::string&&>(name),[&func,&fit](double x){return func({x},fit.Parameters());},min,max,step);
			return *this;
		}
		PlotFit1D& ParamFunc(std::string&&name,paramFunc func,const FIT&fit,double step){
			if(max<min)throw math_h_error<PlotFit1D>("No FitPoints instance initialized the ranges on the plot");
			Plot<double>::Line(static_cast<std::string&&>(name),[&func,&fit](double x){return func({x},fit.Parameters());},min,max,step);
			return *this;
		}
	};
}
#endif

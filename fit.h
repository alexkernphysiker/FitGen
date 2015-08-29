// this file is distributed under 
// GPL v 3.0 license
#ifndef ____lrEPWamH___
#define ____lrEPWamH___
#include <functional>
#include "math_h/gnuplot/gnuplot.h"
#include "genetic_exception.h"
#include "abstract.h"
#include "genetic.h"
namespace Genetic{
	using namespace std;
	class IParamFunc{
	public:
		virtual ~IParamFunc(){}
		virtual double operator()(ParamSet&&X,ParamSet&&P)=0;
		inline double F(ParamSet&&X,ParamSet&P){return operator()(static_cast<ParamSet&&>(X),static_cast<ParamSet&&>(P));}
		inline double F(ParamSet&X,ParamSet&&P){return operator()(static_cast<ParamSet&&>(X),static_cast<ParamSet&&>(P));}
		inline double F(ParamSet&X,ParamSet&P){return operator()(static_cast<ParamSet&&>(X),static_cast<ParamSet&&>(P));}
	};
	typedef function<double(ParamSet&&,ParamSet&&)> paramFunc;
	class ParameterFunction:public IParamFunc{
	public:
		ParameterFunction(paramFunc f);
		virtual ~ParameterFunction();
		virtual double operator()(ParamSet&&X,ParamSet&&P) override;
	private:
		paramFunc func;
	};
	
	class FitPoints{
	public:
		struct Point{
			Point();
			Point(const Point &src);
			Point &operator=(const Point &src);
			ParamSet X;
			ParamSet WX;
			double y;
			double wy;
		};
		FitPoints();
		virtual ~FitPoints();
		FitPoints &operator<<(Point point);
		Point &operator[](int i);
		int count();
		typedef vector<Point>::iterator iterator;
		typedef vector<Point>::const_iterator const_iterator;
		iterator begin();
		const_iterator cbegin()const;
		iterator end();
		const_iterator cend() const;
	private:
		vector<Point> m_data;
	};
	shared_ptr<FitPoints> operator<<(shared_ptr<FitPoints> src,FitPoints::Point p);
	shared_ptr<FitPoints> operator<<(shared_ptr<FitPoints> src,pair<double,double> p);
	shared_ptr<FitPoints> operator<<(shared_ptr<FitPoints> src,pair<double,pair<double,double>> p);
	shared_ptr<FitPoints> operator<<(shared_ptr<FitPoints> src,pair<pair<double,double>,pair<double,double>> p);
	shared_ptr<FitPoints> SelectFitPoints(shared_ptr<FitPoints> src,shared_ptr<IParamCheck> condition);
	shared_ptr<FitPoints> SelectFitPoints(shared_ptr<FitPoints> src,function<bool(double)> Ycond);
	shared_ptr<FitPoints> SelectFitPoints(shared_ptr<FitPoints> src,shared_ptr<IParamCheck> condition,function<bool(double)> Ycond);
	template<class IndexerX,class IndexerY=IndexerX>
	shared_ptr<FitPoints> FitPointsXY(int from,int to,IndexerX X,IndexerY Y){
		auto res=make_shared<FitPoints>();
		for(int i=from;i<=to;i++){
			FitPoints::Point P;
			P.X<<X[i];
			P.y=Y[i];
			res<<P;
		}
		return res;
	}
	template<class IndexerX,class IndexerY=IndexerX,class IndexerWY=IndexerY>
	shared_ptr<FitPoints> FitPointsXYdY(int from,int to,IndexerX X,IndexerY Y,IndexerWY WY){
		auto res=make_shared<FitPoints>();
		for(int i=from;i<=to;i++){
			FitPoints::Point P;
			P.X<<X[i];
			P.y=Y[i];
			P.wy=WY[i];
			res<<P;
		}
		return res;
	}
	template<class IndexerX,class IndexerWX=IndexerX,class IndexerY=IndexerX,class IndexerWY=IndexerY>
	shared_ptr<FitPoints> FitPointsXdXYdY(int from,int to,IndexerX X,IndexerWX WX,IndexerY Y,IndexerWY WY){
		auto res=make_shared<FitPoints>();
		for(int i=from;i<=to;i++){
			FitPoints::Point P;
			P.X<<X[i];
			P.WX<<WX[i];
			P.y=Y[i];
			P.wy=WY[i];
			res<<P;
		}
		return res;
	}
	class Distribution1D:public FitPoints{
	public:
		Distribution1D(double min, double max, int bins);
		void Fill(double x);
	};
	
	class OptimalityForPoints:public IOptimalityFunction{
	public:
		typedef function<double(ParamSet&,IParamFunc&)> Coefficient;
		typedef function<double(FitPoints::Point&,ParamSet&,IParamFunc&)> Summand;
		OptimalityForPoints(shared_ptr<FitPoints> p, shared_ptr<IParamFunc> f,Coefficient c,Summand s);
		virtual ~OptimalityForPoints();
		virtual double operator()(ParamSet&&P)override;
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
		double GetParamParabolicError(double delta,int i);
		ParamSet GetParamParabolicErrors(ParamSet&&delta);
	};
	shared_ptr<IOptimalityFunction> SumSquareDiff(shared_ptr<FitPoints> points, shared_ptr<IParamFunc> f);
	shared_ptr<IOptimalityFunction> SumWeightedSquareDiff(shared_ptr<FitPoints> points, shared_ptr<IParamFunc> f);
	shared_ptr<IOptimalityFunction> ChiSquare(shared_ptr<FitPoints> points, shared_ptr<IParamFunc> f);
	shared_ptr<IOptimalityFunction> ChiSquareWithXError(shared_ptr<FitPoints> points, shared_ptr<IParamFunc> f);
	template<class GENETIC,shared_ptr<IOptimalityFunction> OptimalityAlgorithm(shared_ptr<FitPoints>,shared_ptr<IParamFunc>)>
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
		double operator()(ParamSet&&X){return m_func->operator()(static_cast<ParamSet&&>(X),AbstractGenetic::Parameters());}
	};
	template<class GENETIC,class FUNC,shared_ptr<IOptimalityFunction> OptimalityAlgorithm(shared_ptr<FitPoints>,shared_ptr<IParamFunc>)>
	class FitFunction:public virtual Fit<GENETIC,OptimalityAlgorithm>{
	public:
		typedef FUNC functype;
		FitFunction(shared_ptr<FitPoints> points):Fit<GENETIC,OptimalityAlgorithm>(make_shared<FUNC>()),
			AbstractGenetic(OptimalityAlgorithm(points,make_shared<FUNC>())){}
		virtual ~FitFunction(){}
	};
	
	class OptimalityForPointsWithFuncError:public IOptimalityFunction{
	public:
		typedef function<double(ParamSet&,IParamFunc&,IParamFunc&)> Coefficient;
		typedef function<double(FitPoints::Point&,ParamSet&,IParamFunc&,IParamFunc&)> Summand;
		OptimalityForPointsWithFuncError(shared_ptr<FitPoints> p,shared_ptr<IParamFunc> f,shared_ptr<IParamFunc> e,Coefficient c,Summand s);
		virtual ~OptimalityForPointsWithFuncError();
		virtual double operator()(ParamSet&&P)override;
	protected:
		shared_ptr<FitPoints> points;
		shared_ptr<IParamFunc> func;
		shared_ptr<IParamFunc> error;
		Coefficient C;
		Summand S;
	};
	shared_ptr<IOptimalityFunction> ChiSquare(shared_ptr<FitPoints> points, shared_ptr<IParamFunc> f,shared_ptr<IParamFunc> e);
	shared_ptr<IOptimalityFunction> ChiSquareWithXError(shared_ptr<FitPoints> points, shared_ptr<IParamFunc> f,shared_ptr<IParamFunc> e);
	template<class GENETIC,shared_ptr<IOptimalityFunction> OptimalityAlgorithm(shared_ptr<FitPoints>,shared_ptr<IParamFunc>,shared_ptr<IParamFunc>)>
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
		double operator()(ParamSet&&X){return m_func->operator()(static_cast<ParamSet&&>(X),AbstractGenetic::Parameters());}
	};
	class PlotPoints1D:public Plot<double>{
	public:
		PlotPoints1D();
		virtual ~PlotPoints1D();
		PlotPoints1D& Points(std::string name,shared_ptr<FitPoints> points,unsigned int param_index=0);
		PlotPoints1D& PointsWithoutErrors(std::string name,shared_ptr<FitPoints> points,unsigned int param_index=0);
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
		PlotFit1D& Points(std::string name,shared_ptr<FitPoints> points,unsigned int param_index=0){
			PlotPoints1D::Points(name,points,param_index);
			for(FitPoints::Point p:*points){
				if(p.X[param_index]<min)
					min=p.X[param_index];
				if(p.X[param_index]>max)
					max=p.X[param_index];
			}
			return *this;
		}
		PlotFit1D& PointsWithoutErrors(std::string name,shared_ptr<FitPoints> points,unsigned int param_index=0){
			PlotPoints1D::PointsWithoutErrors(name,points,param_index);
			for(FitPoints::Point p:*points){
				if(p.X[param_index]<min)
					min=p.X[param_index];
				if(p.X[param_index]>max)
					max=p.X[param_index];
			}
			return *this;
		}
		PlotFit1D& Fit(std::string name,FIT&fit,double step=0.1){
			if(max<min)throw GeneticException("No FitPoints instance initialized the ranges on the plot");
			Plot<double>::Line(name,[&fit](double x){return fit(ParamSet(x));},min,max,step);
			return *this;
		}
		PlotFit1D& ParamFunc(std::string name,IParamFunc&&func,FIT&fit,double step=0.1){
			if(max<min)throw GeneticException("No FitPoints instance initialized the ranges on the plot");
			Plot<double>::Line(name,[&func,&fit](double x){return func(ParamSet(x),fit.Parameters());},min,max,step);
			return *this;
		}
		PlotFit1D& ParamFunc(std::string name,paramFunc func,FIT&fit,double step=0.1){
			if(max<min)throw GeneticException("No FitPoints instance initialized the ranges on the plot");
			Plot<double>::Line(name,[&func,&fit](double x){return func(ParamSet(x),fit.Parameters());},min,max,step);
			return *this;
		}
	};
}
#endif

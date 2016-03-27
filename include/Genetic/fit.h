// this file is distributed under 
// MIT license
#ifndef ____lrEPWamH___
#define ____lrEPWamH___
#include <functional>
#include "../math_h/structures.h"
#include "abstract.h"
#include "genetic.h"
#include "paramfunc.h"
namespace Genetic{
	class ParameterFunction:public IParamFunc{
	public:
		ParameterFunction(const paramFunc f);
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
			Point(const ParamSet&x,const double y_);
			Point(const ParamSet&x,const double y_,const double wy_);
			Point(const ParamSet&x,const ParamSet&wx,const double y_,const double wy_);
			Point(const ParamSet&&x,const double y_);
			Point(const ParamSet&&x,const double y_,const double wy_);
			Point(const ParamSet&&x,const ParamSet&&wx,const double y_,const double wy_);
			const ParamSet&X()const;
			const ParamSet&WX()const;
			double y()const;
			double wy()const;
			double&var_y();
			double&var_wy();
		private:
			ParamSet __X;
			ParamSet __WX;
			double __y;
			double __wy;
		};
		FitPoints();
		FitPoints(const MathTemplates::SortedPoints<double>&h);
		FitPoints(const MathTemplates::BiSortedPoints<double>&d);
		FitPoints(const MathTemplates::SortedPoints<MathTemplates::value<double>>&h);
		FitPoints(const MathTemplates::BiSortedPoints<MathTemplates::value<double>>&d);
		virtual ~FitPoints();
		FitPoints&operator<<(const Point&point);
		const Point&operator[](const size_t i)const;
		size_t size()const;
		size_t dimensions()const;
		const ParamSet&min()const;
		const ParamSet&max()const;
		double Ymin()const;
		double Ymax()const;
		typedef std::vector<Point>::const_iterator const_iterator;
		const_iterator begin()const;
		const_iterator cbegin()const;
		const_iterator end()const;
		const_iterator cend() const;
		const MathTemplates::SortedPoints<MathTemplates::value<double>> Hist1(const size_t parameter_index_x)const;
		const MathTemplates::SortedPoints<MathTemplates::value<double>> Hist1(const size_t parameter_index_x,const size_t parameter_index_y)const;
		const MathTemplates::SortedPoints<double> Line(const size_t parameter_index_x)const;
		const MathTemplates::SortedPoints<double> Line(const size_t parameter_index_x,const size_t parameter_index_y)const;
	private:
		std::vector<Point> m_data;
		ParamSet m_min,m_max;
		double ymin,ymax;
	};
	typedef FitPoints::Point Point;
	std::shared_ptr<FitPoints> operator<<(std::shared_ptr<FitPoints> src,const Point&p);
	std::shared_ptr<FitPoints> operator<<(std::shared_ptr<FitPoints> src,const MathTemplates::point<double>&p);
	std::shared_ptr<FitPoints> operator<<(std::shared_ptr<FitPoints> src,const Point&&p);
	std::shared_ptr<FitPoints> operator<<(std::shared_ptr<FitPoints> src,const MathTemplates::point<double>&&p);
	std::shared_ptr<FitPoints> operator<<(std::shared_ptr<FitPoints> src,const std::shared_ptr<FitPoints>data);
	std::shared_ptr<FitPoints> SelectFitPoints(std::shared_ptr<FitPoints> src,const std::shared_ptr<IParamCheck> condition);
	std::shared_ptr<FitPoints> SelectFitPoints(std::shared_ptr<FitPoints> src,const std::function<bool(double)> Ycond);
	std::shared_ptr<FitPoints> SelectFitPoints(std::shared_ptr<FitPoints> src,const std::shared_ptr<IParamCheck> condition,const std::function<bool(double)> Ycond);
	template<class IndexerX,class IndexerY=IndexerX>
	std::shared_ptr<FitPoints> FitPointsXY(const int from,const int to,const IndexerX&X,const IndexerY&Y){
		auto res=std::make_shared<FitPoints>();
		for(int i=from;i<=to;i++)
			res<<FitPoints::Point({X[i]},Y[i]);
		return res;
	}
	template<class IndexerX,class IndexerY=IndexerX,class IndexerWY=IndexerY>
	std::shared_ptr<FitPoints> FitPointsXYdY(const int from,const int to,const IndexerX&X,const IndexerY&Y,const IndexerWY&WY){
		auto res=std::make_shared<FitPoints>();
		for(int i=from;i<=to;i++)
			res<<FitPoints::Point({X[i]},Y[i],WY[i]);
		return res;
	}
	template<class IndexerX,class IndexerWX=IndexerX,class IndexerY=IndexerX,class IndexerWY=IndexerY>
	std::shared_ptr<FitPoints> FitPointsXdXYdY(const int from,const int to,const IndexerX&X,const IndexerWX&WX,const IndexerY&Y,const IndexerWY&WY){
		auto res=std::make_shared<FitPoints>();
		for(int i=from;i<=to;i++)
			res<<FitPoints::Point({X[i]},{WX[i]},Y[i],WY[i]);
		return res;
	}
	
	class OptimalityForPoints:public IOptimalityFunction{
	public:
		typedef std::function<double(const ParamSet&,const IParamFunc&)> Coefficient;
		typedef std::function<double(const FitPoints::Point&,const ParamSet&,const IParamFunc&)> Summand;
		OptimalityForPoints(const std::shared_ptr<FitPoints> p,const  std::shared_ptr<IParamFunc> f,const Coefficient c,const Summand s);
		virtual ~OptimalityForPoints();
		virtual double operator()(const ParamSet&P)const override;
		std::shared_ptr<FitPoints> Points()const;
	protected:
		std::shared_ptr<FitPoints> points;
		std::shared_ptr<IParamFunc> func;
		Coefficient C;
		Summand S;
	};
	class Parabolic:public virtual AbstractGenetic{
	protected:
		Parabolic();
	public:
		virtual ~Parabolic();
		Parabolic&SetUncertaintyCalcDeltas(const ParamSet&P);
		Parabolic&SetUncertaintyCalcDeltas(const ParamSet&&P);
		const std::vector<MathTemplates::value<double>>&ParametersWithUncertainties()const;
	private:
		double GetParamParabolicError(const double delta, const size_t i)const;
		ParamSet m_delta;
		std::shared_ptr<std::vector<MathTemplates::value<double>>> m_uncertainty_cache;
		std::shared_ptr<unsigned long long int> m_iter_number;
	};
	std::shared_ptr<OptimalityForPoints> SumSquareDiff(const std::shared_ptr<FitPoints> points, const std::shared_ptr<IParamFunc> f);
	std::shared_ptr<OptimalityForPoints> SumWeightedSquareDiff(const std::shared_ptr<FitPoints> points, const std::shared_ptr<IParamFunc> f);
	std::shared_ptr<OptimalityForPoints> ChiSquare(const std::shared_ptr<FitPoints> points, const std::shared_ptr<IParamFunc> f);
	std::shared_ptr<OptimalityForPoints> ChiSquareWithXError(const std::shared_ptr<FitPoints> points, const std::shared_ptr<IParamFunc> f);
	template<class GENETIC,std::shared_ptr<OptimalityForPoints> OptimalityAlgorithm(const std::shared_ptr<FitPoints>,const std::shared_ptr<IParamFunc>)>
	class Fit:public virtual GENETIC,public virtual Parabolic{
	private:
		std::shared_ptr<IParamFunc> m_func;
	protected:
		Fit(std::shared_ptr<IParamFunc> f):GENETIC(),Parabolic(){
			m_func=f;
		}
	public:
		Fit(
			const std::shared_ptr<FitPoints> points, 
			const std::shared_ptr<IParamFunc> f
		):AbstractGenetic(OptimalityAlgorithm(points,f)){
			m_func=f;
		}
		Fit(const std::shared_ptr<FitPoints> points,const paramFunc f):Fit(points,std::make_shared<ParameterFunction>(f)){}
		virtual ~Fit(){}
		double operator()(const ParamSet&X)const{return m_func->operator()(X,AbstractGenetic::Parameters());}
		double operator()(const ParamSet&&X)const{return operator()(X);}
		std::shared_ptr<IParamFunc> Func()const{return m_func;}
		std::shared_ptr<FitPoints> Points()const{
			return std::dynamic_pointer_cast<OptimalityForPoints>(AbstractGenetic::OptimalityCalculator())->Points();
		}
	};
	template<class GENETIC,class FUNC,std::shared_ptr<OptimalityForPoints> OptimalityAlgorithm(const std::shared_ptr<FitPoints>,const std::shared_ptr<IParamFunc>)>
	class FitFunction:public virtual Fit<GENETIC,OptimalityAlgorithm>{
	public:
		typedef FUNC functype;
		FitFunction(const std::shared_ptr<FitPoints> points):
			AbstractGenetic(OptimalityAlgorithm(points,std::make_shared<FUNC>())),
			Fit<GENETIC,OptimalityAlgorithm>(std::make_shared<FUNC>()){}
		virtual ~FitFunction(){}
	};
	
	
	
	class OptimalityForPointsWithFuncError:public IOptimalityFunction{
	public:
		typedef std::function<MathTemplates::value<double>(const ParamSet&,const ParamSet&)> Func;
		typedef std::function<double(const ParamSet&,const Func&)> Coefficient;
		typedef std::function<double(const FitPoints::Point&,const ParamSet&,const Func&)> Summand;
		OptimalityForPointsWithFuncError(const std::shared_ptr<FitPoints> p,const Func f,const Coefficient c,const Summand s);
		virtual ~OptimalityForPointsWithFuncError();
		virtual double operator()(const ParamSet&P)const override;
		std::shared_ptr<FitPoints> Points()const;
	protected:
		std::shared_ptr<FitPoints> points;
		Func m_func;
		Coefficient C;
		Summand S;
	};
	std::shared_ptr<OptimalityForPointsWithFuncError> ChiSquare(const std::shared_ptr<FitPoints> points, const OptimalityForPointsWithFuncError::Func f);
	std::shared_ptr<OptimalityForPointsWithFuncError> ChiSquareWithXError(const std::shared_ptr<FitPoints> points, const OptimalityForPointsWithFuncError::Func f);
	template<class GENETIC,std::shared_ptr<OptimalityForPointsWithFuncError> OptimalityAlgorithm(const std::shared_ptr<FitPoints>, const OptimalityForPointsWithFuncError::Func f)>
	class FitFunctionWithError:public virtual GENETIC,public virtual Parabolic{
	private:
		OptimalityForPointsWithFuncError::Func m_func;
	public:
		typedef typename OptimalityForPointsWithFuncError::Func Func;
		FitFunctionWithError(std::shared_ptr<FitPoints> points, const Func f)
			:AbstractGenetic(OptimalityAlgorithm(points,f)),GENETIC(),Parabolic(){
			m_func=f;
		}
		virtual ~FitFunctionWithError(){}
		const MathTemplates::value<double> operator()(ParamSet&&X)const{return m_func(X,AbstractGenetic::Parameters());}
		std::shared_ptr<FitPoints> Points()const{
			return std::dynamic_pointer_cast<OptimalityForPointsWithFuncError>(AbstractGenetic::OptimalityCalculator())->Points();
		}
	};
}
#endif

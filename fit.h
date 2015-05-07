#ifndef ____lrEPWamH___
#define ____lrEPWamH___
#include <functional>
#include "abstract.h"
#include "genetic.h"
#include "genetic_exception.h"
namespace Genetic{
	using namespace std;
	class IParamFunc{
	public:
		virtual ~IParamFunc(){}
		virtual double operator()(ParamSet&X,ParamSet&P)=0;
	};
	class ParameterFunction:public IParamFunc{
	public:
		ParameterFunction(function<double(ParamSet &,ParamSet &)> f);
		virtual ~ParameterFunction();
		virtual double operator()(ParamSet&X,ParamSet&P) override;
	private:
		function<double(ParamSet &,ParamSet &)> func;
	};
	
	class FitPoints{
	public:
		struct Point{
		public:
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
		virtual double operator()(ParamSet&P)override;
	protected:
		shared_ptr<FitPoints> points;
		shared_ptr<IParamFunc> func;
		Coefficient C;
		Summand S;
	};
	template<class GENETIC>
	class Parabolic:public GENETIC{
	public:
		Parabolic(shared_ptr<IOptimalityFunction> opt):GENETIC(opt){}
		virtual ~Parabolic(){}
		double GetParamParabolicError(double delta,int i){
			if(delta<=0)
				throw new GeneticException("Error in parabolic error calculation: delta cannot be zero or negative");
			double s=AbstractGenetic::Optimality();
			ParamSet ab=AbstractGenetic::Parameters();
			ParamSet be=ab;
			ab.Set(i,ab[i]+delta);
			be.Set(i,be[i]-delta);
			double sa=AbstractGenetic::OptimalityCalculator()->operator()(ab);
			double sb=AbstractGenetic::OptimalityCalculator()->operator()(be);
			double da=(sa-s)/delta;
			double db=(s-sb)/delta;
			double dd=(da-db)/delta;
			if(dd<=0)
				return INFINITY;
			else
				return sqrt(2.0/dd);
		}
		ParamSet GetParamParabolicErrors(ParamSet delta){
			ParamSet res;
			for(int i=0,n=AbstractGenetic::ParamCount();i<n;i++)
				res<<GetParamParabolicError(delta[i],i);
			return res;
		}
	};
	shared_ptr<IOptimalityFunction> SumSquareDiff(shared_ptr<FitPoints> points, shared_ptr<IParamFunc> f);
	shared_ptr<IOptimalityFunction> SumWeightedSquareDiff(shared_ptr<FitPoints> points, shared_ptr<IParamFunc> f);
	shared_ptr<IOptimalityFunction> ChiSquare(shared_ptr<FitPoints> points, shared_ptr<IParamFunc> f);
	shared_ptr<IOptimalityFunction> ChiSquareWithXError(shared_ptr<FitPoints> points, shared_ptr<IParamFunc> f);
	template<class GENETIC,shared_ptr<IOptimalityFunction> OptimalityAlgorithm(shared_ptr<FitPoints>,shared_ptr<IParamFunc>)>
	class Fit:public Parabolic<GENETIC>{
	private:
		shared_ptr<IParamFunc> m_func;
	public:
		Fit(
			shared_ptr<FitPoints> points, 
			shared_ptr<IParamFunc> f
		):Parabolic<GENETIC>(OptimalityAlgorithm(points,f)){
			m_func=f;
		}
		Fit(
			shared_ptr<FitPoints> points,
			function<double(ParamSet &,ParamSet &)> f
		):Fit(points,make_shared<ParameterFunction>(f)){}
		virtual ~Fit(){}
		double operator()(ParamSet X){
			return m_func->operator()(X,AbstractGenetic::Parameters());
		}
	};
	template<class GENETIC,class FUNC,shared_ptr<IOptimalityFunction> OptimalityAlgorithm(shared_ptr<FitPoints>,shared_ptr<IParamFunc>)>
	class FitFunction:public Fit<GENETIC,OptimalityAlgorithm>{
	public:
		FitFunction(shared_ptr<FitPoints> points):Fit<GENETIC,OptimalityAlgorithm>(points,make_shared<FUNC>()){}
		virtual ~FitFunction(){}
	};
	
	class OptimalityForPointsWithFuncError:public IOptimalityFunction{
	public:
		typedef function<double(ParamSet&,IParamFunc&,IParamFunc&)> Coefficient;
		typedef function<double(FitPoints::Point&,ParamSet&,IParamFunc&,IParamFunc&)> Summand;
		OptimalityForPointsWithFuncError(shared_ptr<FitPoints> p,shared_ptr<IParamFunc> f,shared_ptr<IParamFunc> e,Coefficient c,Summand s);
		virtual ~OptimalityForPointsWithFuncError();
		virtual double operator()(ParamSet&P)override;
	protected:
		shared_ptr<FitPoints> points;
		shared_ptr<IParamFunc> func;
		shared_ptr<IParamFunc> error;
		Coefficient C;
		Summand S;
	};
	shared_ptr<IOptimalityFunction> ChiSquare(shared_ptr<FitPoints> points, shared_ptr<IParamFunc> f,shared_ptr<IParamFunc> e);
	shared_ptr<IOptimalityFunction> ChiSquareWithXError(shared_ptr<FitPoints> points, shared_ptr<IParamFunc> f,shared_ptr<IParamFunc> e);
	template<class GENETIC,shared_ptr<IOptimalityFunction> OptiimalityAlgorithm(shared_ptr<FitPoints>,shared_ptr<IParamFunc>,shared_ptr<IParamFunc>)>
	class FitFunctionWithError:public Parabolic<GENETIC>{
	private:
		shared_ptr<IParamFunc> m_func;
	public:
		FitFunctionWithError(shared_ptr<FitPoints> points,shared_ptr<IParamFunc> f,shared_ptr<IParamFunc> e)
			:Parabolic<GENETIC>(OptiimalityAlgorithm(points,f,e)){
			m_func=f;
		}
		FitFunctionWithError(
			shared_ptr<FitPoints> points,
			function<double(ParamSet &,ParamSet &)> f,
			function<double(ParamSet &,ParamSet &)> e
		):FitFunctionWithError(points,make_shared<ParameterFunction>(f),make_shared<ParameterFunction>(e)){}
		virtual ~FitFunctionWithError(){}
		double operator()(ParamSet X){
			return m_func->operator()(X,AbstractGenetic::Parameters());
		}
	};
}
#endif

#ifndef ____lrEPWamH___
#define ____lrEPWamH___
#include <functional>
#include "fit_gen.h"
#include "fitexception.h"
namespace Fit{
	using namespace std;
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
		OptimalityForPoints(shared_ptr<FitPoints> p,Coefficient c,Summand s);
		virtual ~OptimalityForPoints();
		virtual double operator()(ParamSet&P,IParamFunc&F)override;
	protected:
		shared_ptr<FitPoints> points;
		Coefficient C;
		Summand S;
	};
	shared_ptr<OptimalityForPoints> SumSquareDiff(shared_ptr<FitPoints> points);
	shared_ptr<OptimalityForPoints> SumWeightedSquareDiff(shared_ptr<FitPoints> points);
	shared_ptr<OptimalityForPoints> ChiSquare(shared_ptr<FitPoints> points);
	shared_ptr<OptimalityForPoints> ChiSquareWithXError(shared_ptr<FitPoints> points);
}
#endif

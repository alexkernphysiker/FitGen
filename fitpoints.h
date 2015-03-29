#ifndef ____lrEPWamH___
#define ____lrEPWamH___
#include "fit_gen.h"
#include "fitexception.h"
namespace Fit{
	using namespace std;
	class FitPointsAbstract:public IOptimalityFunction{
	public:
		virtual ~FitPointsAbstract();
		FitPointsAbstract &Add(ParamSet x, double y, double weight=1);
		FitPointsAbstract &Add(ParamSet x, ParamSet x_w, double y, double weight=1);
		virtual double operator()(ParamSet params, IParamFunc &func)override=0;
		int Count();
		ParamSet X(int i);
		ParamSet X_w(int i);
		double Y(int i);
		double W(int i);
	protected:
		vector<ParamSet> m_data;
		vector<ParamSet> m_data_w;
		vector<double> m_y;
		vector<double> m_w;
	};
	class SquareDiff:public FitPointsAbstract{
	public:
		SquareDiff(){}
		virtual ~SquareDiff(){}
		virtual double operator()(ParamSet params, IParamFunc &func)override;
	};
	class ChiSquare:public FitPointsAbstract{
	public:
		ChiSquare(){}
		virtual ~ChiSquare(){}
		virtual double operator()(ParamSet params, IParamFunc &func)override;
	};
	class ChiSquareWithXError:public FitPointsAbstract{
	public:
		ChiSquareWithXError(){}
		virtual ~ChiSquareWithXError(){}
		virtual double operator()(ParamSet params, IParamFunc &func)override;
	};
	
	template<class fitpoints>
	shared_ptr<fitpoints> SelectFitPoints(shared_ptr<FitPointsAbstract> data, IParamCheck &condition){
		shared_ptr<fitpoints> res=make_shared<fitpoints>();
		for(int i=0; i<data->Count();i++){
			if(condition.CorrectParams(data->X(i)))
				res->Add(data->X(i),data->Y(i),data->W(i));
		}
		return res;
	}
	template <class fitpoints, class func>
	shared_ptr<fitpoints> SelectFitPointsByY(shared_ptr<FitPointsAbstract> data, IParamCheck &condition, func Ycond){
		shared_ptr<fitpoints> res=make_shared<fitpoints>();
		for(int i=0; i<data->Count();i++){
			if(condition.CorrectParams(data->X(i)) && Ycond(data->Y(i)))
				res->Add(data->X(i),data->Y(i),data->W(i));
		}
		return res;
	}
	template <class fitpoints, class indexerx, class indexery = indexerx>
	shared_ptr<fitpoints> FitPointsXY(int from, int to, indexerx X, indexery Y){
		if(to<from)throw;
		auto res=make_shared<fitpoints>();
		for(int i=from; i<=to;i++)
			res->Add(ParamSet(X[i]),Y[i]);
		return res;
	}
	template <class fitpoints, class indexerx, class indexery = indexerx, class indexererr = indexery>
	shared_ptr<fitpoints> FitPointsXYdY(int from, int to, indexerx X, indexery Y,indexererr W){
		if(to<from)throw;
		auto res=make_shared<fitpoints>();
		for(int i=from; i<=to;i++)
			res->Add(ParamSet(X[i]),Y[i],W[i]);
		return res;
	}
	template <class fitpoints, class indexerx, class indexerx_err=indexerx, class indexery = indexerx, class indexererr = indexery>
	shared_ptr<fitpoints> FitPointsXdXYdY(int from, int to, indexerx X, indexerx_err X_w, indexery Y,indexererr W){
		if(to<from)throw;
		auto res=make_shared<fitpoints>();
		for(int i=from; i<=to;i++)
			res->Add(ParamSet()<<X[i],ParamSet()<<X_w[i],Y[i],W[i]);
		return res;
	}
	template <class fitpoints>
	class Distribution1D:public fitpoints{
		private:double m_min;double m_max;double binwidth;
	public:
		Distribution1D(double min, double max, int bins):fitpoints(){
			if(max<min)throw;if(bins<1)throw;
			m_max=max;m_min=min;
			binwidth=(max-min)/double(bins);
			double offs=binwidth/2.0;
			for(double x=min+offs;x<max;x+=binwidth)
				fitpoints::Add(ParamSet(x),ParamSet(offs),double(0),double(1));
		}
		void Fill(double x){
			int bin_pos=int((x-m_min)/binwidth);
			if((bin_pos>=0)&&(bin_pos<fitpoints::Count())){
				fitpoints::m_y[bin_pos]+=1;
				fitpoints::m_w[bin_pos]=sqrt(fitpoints::m_y[bin_pos]);
			}
		}
	};
}
#endif

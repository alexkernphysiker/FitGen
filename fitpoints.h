#ifndef __APROPOINTS__
#define __APROPOINTS__
#include "fit_gen.h"
#include "fitexception.h"
namespace Fit{
class FitPointsAbstract:public IOptimalityFunction{
public:
	virtual ~FitPointsAbstract();
	FitPointsAbstract &Add(ParamSet x, double y, double weight=1);
	FitPointsAbstract &Add(ParamSet x, ParamSet x_w, double y, double weight=1);
	virtual double operator()(ParamSet &params, IParamFunc &func)override=0;
	int Count();
	ParamSet &X(int i);// gets measured parameters for the points
	ParamSet &X_w(int i);// gets measured parameters errors for the points (if used)
	double Y(int i);// gets the value
	double W(int i);// gets the point's weight/error
protected:
	std::vector<ParamSet> m_data;
	std::vector<ParamSet> m_data_w;
	std::vector<double> m_y;
	std::vector<double> m_w;
};
class SquareDiff:public FitPointsAbstract{
public:
	SquareDiff(){}
	virtual ~SquareDiff(){}
	virtual double operator()(ParamSet &params, IParamFunc &func)override;
};
class xi_2:public FitPointsAbstract{
public:
	xi_2(){}
	virtual ~xi_2(){}
	virtual double operator()(ParamSet &params, IParamFunc &func)override;
};
class xi_2_wx:public FitPointsAbstract{
public:
	xi_2_wx(){}
	virtual ~xi_2_wx(){}
	virtual double operator()(ParamSet &params, IParamFunc &func)override;
};
template<class fitpoints,class fitpoints2=fitpoints>
std::shared_ptr<fitpoints> SelectFitPoints(std::shared_ptr<fitpoints2> data, IParamCheck &condition){
	std::shared_ptr<fitpoints> res=std::make_shared<fitpoints>();
	for(uint i=0; i<data->Count();i++){
		if(condition.CorrectParams(data->X(i)))
			res->Add(data->X(i),data->Y(i),data->W(i));
	}
	return res;
}
template <class fitpoints, class func,class fitpoints2=fitpoints>
std::shared_ptr<fitpoints> SelectFitPointsByY(std::shared_ptr<fitpoints2> data, IParamCheck &condition, func Ycond){
	std::shared_ptr<fitpoints> res=std::make_shared<fitpoints>();
	for(uint i=0; i<data->Count();i++){
		if(condition.CorrectParams(data->X(i)) && Ycond(data->Y(i)))
			res->Add(data->X(i),data->Y(i),data->W(i));
	}
	return res;
}
template <class fitpoints, class indexerx, class indexery = indexerx>
class FitPoints:public fitpoints{
public:
	FitPoints(int from, int to, indexerx X, indexery Y):fitpoints(){
		if(to<from)throw;
		for(int i=from; i<=to;i++)
			this->Add(ParamSet()<<X[i],Y[i]);
	}
};
template <class fitpoints, class indexerx, class indexery = indexerx, class indexererr = indexery>
class FitPointsWithErrors:public fitpoints{
public:
	FitPointsWithErrors(int from, int to, indexerx X, indexery Y,indexererr W):fitpoints(){
		if(to<from)throw;
		for(int i=from; i<=to;i++)
			this->Add(ParamSet()<<X[i],Y[i],W[i]);
	}
};
template <class fitpoints, class indexerx, class indexerx_err=indexerx, class indexery = indexerx, class indexererr = indexery>
class FitPointsWithXErrors:public fitpoints{
public:
	FitPointsWithXErrors(int from, int to, indexerx X, indexerx_err X_w, indexery Y,indexererr W):fitpoints(){
		if(to<from)throw;
		for(int i=from; i<=to;i++)
			this->Add(ParamSet()<<X[i],ParamSet()<<X_w[i],Y[i],W[i]);
	}
};
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
			fitpoints::Add(ParamSet()<<x,double(0));
	}
	void AddValue(double x){
		int bin_pos=int((x-m_min)/binwidth);
		if((bin_pos>=0)&&(bin_pos<fitpoints::Count()))
			fitpoints::m_y[bin_pos]+=1;
	}
};
}
#endif // __APROPOINTS__

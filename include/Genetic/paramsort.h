// this file is distributed under 
// MIT license
#ifndef PUWJZCZDMORMZODA_sort
#define PUWJZCZDMORMZODA_sort
#include <memory>
#include "../math_h/error.h"
#include "../gnuplot_wrap.h"
#include "paramset.h"
namespace Genetic{
	using namespace std;
	using namespace MathTemplates;
	using namespace GnuplotWrap;
	class BinningParam{
	public:
		BinningParam(size_t ind,const pair<double,double>&R,size_t cnt);
		BinningParam(size_t ind,pair<double,double>&&R,size_t cnt);
		BinningParam(const BinningParam&source);
		size_t param_index()const;
		size_t count()const;
		pair<double,double>&limits()const;
		double bin_width()const;
		double bin_center(size_t i)const;
	protected:
		bool FindBinIndex(const ParamSet& P,size_t&res)const;
	private:
		void CheckCorrectness()const;
		size_t index,cnt;
		pair<double,double> range;
	};
	class BinChecker:public BinningParam{
	public:
		BinChecker(size_t ind,const pair<double,double>&R,size_t cnt):BinningParam(ind,R,cnt){}
		BinChecker(size_t ind,pair<double,double>&&R,size_t cnt):BinningParam(ind,R,cnt){}
		BinChecker(const BinChecker&source):BinningParam(source){}
		bool FindBin(const ParamSet& P,size_t&res)const{return FindBinIndex(P,res);}
	};
	template<class ParamProcessor>
	class AbstractPerBinSeparator:public BinningParam{
	private:
		vector<shared_ptr<ParamProcessor>> m_data;
	protected:
		virtual shared_ptr<ParamProcessor> CreateParamProcessor(size_t for_pos)=0;
		virtual void ProcessParams(ParamProcessor&proc,const ParamSet&P)=0;
		AbstractPerBinSeparator(size_t ind,const pair<double,double>&R,size_t cnt):BinningParam(ind,R,cnt){}
		AbstractPerBinSeparator(size_t ind,pair<double,double>&&R,size_t cnt):AbstractPerBinSeparator(ind,R,cnt){}
		AbstractPerBinSeparator(const BinningParam& source):AbstractPerBinSeparator(source.param_index(),source.limits(),source.count()){}
		AbstractPerBinSeparator(BinningParam&&source):AbstractPerBinSeparator(source.param_index(),source.limits(),source.count()){}
		AbstractPerBinSeparator(const AbstractPerBinSeparator&source):BinningParam(source){m_data=source.m_data;}
		void Init(){
			for(size_t i=0;i<count();i++)
				m_data.push_back(CreateParamProcessor(i));
		}
	public:
		ParamProcessor&operator[](size_t index)const{
			if(m_data.size()==0)throw Exception<AbstractPerBinSeparator>("Not inited");
			if(index>=count())throw Exception<AbstractPerBinSeparator>("Range check error");
			return const_cast<ParamProcessor&>(*(m_data[index]));
		}
		AbstractPerBinSeparator&operator<<(const ParamSet&P){
			if(m_data.size()==0)throw Exception<AbstractPerBinSeparator>("Not inited");
			size_t index=0;
			if(FindBinIndex(P,index))
				ProcessParams(*(m_data[index]),P);
			return *this;
		}
		AbstractPerBinSeparator&operator<<(ParamSet&&P){return operator<<(P);}
	};
	
	class ParamsPerBins:public AbstractPerBinSeparator<vector<ParamSet>>{
	public:
		ParamsPerBins(size_t ind,const pair<double,double>&R,size_t cnt);
		ParamsPerBins(size_t ind,pair<double,double>&&R,size_t cnt);
		ParamsPerBins(const BinningParam&src);
		ParamsPerBins(BinningParam&&src);
	protected:
		virtual shared_ptr<vector<ParamSet>> CreateParamProcessor(size_t for_pos)override;
		virtual void ProcessParams(vector<ParamSet>&proc,const ParamSet&P)override;
	};
	
	template<unsigned int dimensions>
	class ParamsPerBinsCounter:public AbstractPerBinSeparator<ParamsPerBinsCounter<(dimensions-1)>>{
	public:
		friend class ParamsPerBinsCounter<dimensions+1>;
		typedef function<void(const ParamSet&,const unsigned long)> Delegate;
	private:
		vector<BinningParam>*m_binning_addr;
		vector<BinningParam> m_binning;
	protected:
		virtual shared_ptr<ParamsPerBinsCounter<(dimensions-1)>> CreateParamProcessor(size_t)override{
			return make_shared<ParamsPerBinsCounter<(dimensions-1)>>(m_binning_addr);
		}
		virtual void ProcessParams(ParamsPerBinsCounter<(dimensions-1)>&proc,const ParamSet&P)override{proc<<P;}
		void Full_Cycle(Delegate func,ParamSet&P)const{
			for(size_t i=0;i<AbstractPerBinSeparator<ParamsPerBinsCounter<(dimensions-1)>>::count();i++){
				P[AbstractPerBinSeparator<ParamsPerBinsCounter<(dimensions-1)>>::param_index()]=
					AbstractPerBinSeparator<ParamsPerBinsCounter<(dimensions-1)>>::bin_center(i);
				AbstractPerBinSeparator<ParamsPerBinsCounter<(dimensions-1)>>::operator[](i).Full_Cycle(func,P);
			}
		}
	public:
		ParamsPerBinsCounter<dimensions>(vector<BinningParam>*binning)
			:AbstractPerBinSeparator<ParamsPerBinsCounter<(dimensions-1)>>((*binning)[binning->size()-dimensions]){
				m_binning_addr=binning;
				AbstractPerBinSeparator<ParamsPerBinsCounter<(dimensions-1)>>::Init();
		}
		ParamsPerBinsCounter<dimensions>(const vector<BinningParam>&binning)
			:AbstractPerBinSeparator<ParamsPerBinsCounter<(dimensions-1)>>(binning[binning.size()-dimensions]){
				m_binning=binning;
				m_binning_addr=&m_binning;
				AbstractPerBinSeparator<ParamsPerBinsCounter<(dimensions-1)>>::Init();
		}
		ParamsPerBinsCounter<dimensions>(vector<BinningParam>&&binning):ParamsPerBinsCounter<dimensions>(binning){}
		void FullCycle(Delegate func)const{
			size_t cnt=0;
			for(const BinningParam&bp: *m_binning_addr)
				if(cnt<bp.param_index())
					cnt=bp.param_index();
			ParamSet P=parZeros(cnt+1);
			Full_Cycle(func,P);
		}
	};
	template<>
	class ParamsPerBinsCounter<1>:public AbstractPerBinSeparator<unsigned long>{
	public:
		friend class ParamsPerBinsCounter<2>;
		typedef function<void(const ParamSet&,const unsigned long)> Delegate;
		ParamsPerBinsCounter(const vector<BinningParam>&binning);
		ParamsPerBinsCounter(vector<BinningParam>&&binning);
		ParamsPerBinsCounter(const vector<BinningParam>*binning);
		void FullCycle(Delegate func)const;
	protected:
		virtual shared_ptr<unsigned long> CreateParamProcessor(size_t)override;
		virtual void ProcessParams(unsigned long&proc,const ParamSet&)override;
		void Full_Cycle(Delegate func,ParamSet&P)const;
	};
	
	typedef PlotPoints<double,vector<pair<double,double>>> PlotEngine;
	class AbstractPlotStream{
	protected:
		AbstractPlotStream(string&&name,shared_ptr<PlotEngine>plot);
		AbstractPlotStream(string&&name);
		virtual ~AbstractPlotStream();
	public:
		AbstractPlotStream&operator<<(const ParamSet&P);
		AbstractPlotStream&operator<<(ParamSet&&P);
		shared_ptr<PlotEngine>Plot();
	protected:
		virtual void ProcessPoint(const ParamSet&P)=0;
		string&&Name()const;
	private:
		string m_name;
		shared_ptr<PlotEngine> m_plot;
	};
	class SimplePlotStream:public AbstractPlotStream{
	public:
		SimplePlotStream(string&& name,pair<size_t,size_t>&&indexes);
		SimplePlotStream(string&& name,pair<size_t,size_t>&&indexes,shared_ptr<PlotEngine> plot);
		virtual ~SimplePlotStream();
		typedef function<double(double)> func;
		SimplePlotStream&AddFunc(func f);
	protected:
		virtual void ProcessPoint(const ParamSet& P)override;
	private:
		pair<size_t,size_t> m_indexes;
		vector<pair<double,double>> m_data;
		pair<double,double> m_xrange;
		vector<func> m_funcs;
	};
};
#endif
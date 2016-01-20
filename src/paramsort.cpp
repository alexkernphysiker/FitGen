// this file is distributed under 
// MIT license
#include <Genetic/paramsort.h>
namespace Genetic{
	using namespace std;
	using namespace MathTemplates;
	using namespace GnuplotWrap;
	BinningParam::BinningParam(size_t ind, const pair< double, double >& R, size_t c){
		index=ind;cnt=c;
		range=R;
		CheckCorrectness();
	}
	BinningParam::BinningParam(size_t ind, pair< double, double >&& R, size_t cnt):BinningParam(ind,R,cnt){}
	BinningParam::BinningParam(const BinningParam& source):BinningParam(source.index,source.range,source.cnt){}
	pair<double,double>& BinningParam::limits()const{return const_cast<pair<double,double>&>(range);}
	size_t BinningParam::param_index() const{return index;}
	size_t BinningParam::count() const{return cnt;}
	void BinningParam::CheckCorrectness() const{
		if(range.second<=range.first)
			throw Exception<BinningParam>("wrong binning ranges");
		if(0==cnt)
			throw Exception<BinningParam>("there cannot be zero bins");
	}
	double BinningParam::bin_width() const{
		if(0==count())throw Exception<BinningParam>("count==0");
		return (range.second-range.first)/double(count());
	}
	double BinningParam::bin_center(size_t i) const{
		if(i>=count())throw Exception<BinningParam>("bin range check error");
		return range.first+bin_width()*(double(i)+0.5);
	}
	bool BinningParam::FindBinIndex(const ParamSet& P, size_t& res) const{
		double x=P[index],delta=bin_width()/2.0;
		if(delta<=0)throw Exception<BinningParam>("delta<=0");
		for(size_t i=0,n=count();i<n;i++){
			double pos=bin_center(i);
			if((x>=(pos-delta))&&(x<(pos+delta))){
				res=i;
				return true;
			}
		}
		return false;
	}
	
	ParamsPerBins::ParamsPerBins(size_t ind, const pair<double,double>& R, size_t cnt):AbstractPerBinSeparator<vector<ParamSet>>(ind,R,cnt){Init();}
	ParamsPerBins::ParamsPerBins(size_t ind, pair< double, double >&& R, size_t cnt):ParamsPerBins(ind,R,cnt){}
	ParamsPerBins::ParamsPerBins(const BinningParam& src):AbstractPerBinSeparator<vector<ParamSet>>(src){Init();}
	ParamsPerBins::ParamsPerBins(BinningParam&&src):ParamsPerBins(src){}
	
	shared_ptr<vector<ParamSet>> ParamsPerBins::CreateParamProcessor(size_t){return make_shared<vector<ParamSet>>();}
	void ParamsPerBins::ProcessParams(vector<ParamSet>& proc, const ParamSet& P){proc.push_back(P);}
	
	ParamsPerBinsCounter<1>::ParamsPerBinsCounter(const vector<BinningParam>& binning):AbstractPerBinSeparator<unsigned long>(binning[binning.size()-1]){Init();}
	ParamsPerBinsCounter<1>::ParamsPerBinsCounter(vector< BinningParam >&& binning):ParamsPerBinsCounter(binning){}
	ParamsPerBinsCounter<1>::ParamsPerBinsCounter(const vector<BinningParam>*binning):ParamsPerBinsCounter(*binning){}
	shared_ptr<unsigned long> ParamsPerBinsCounter<1>::CreateParamProcessor(size_t){return make_shared<unsigned long>(0);}
	void ParamsPerBinsCounter<1>::ProcessParams(long unsigned int& proc, const ParamSet&){proc++;}
	void ParamsPerBinsCounter<1>::Full_Cycle(ParamsPerBinsCounter<int(1)>::Delegate func,ParamSet&P)const{
		for(size_t i=0;i<AbstractPerBinSeparator<unsigned long>::count();i++){
			double x=AbstractPerBinSeparator<unsigned long>::bin_center(i);
			P<<x;
			func(P,AbstractPerBinSeparator<unsigned long>::operator[](i));
			P>>x;
		}
	}
	void ParamsPerBinsCounter<1>::FullCycle(Genetic::ParamsPerBinsCounter<int(1)>::Delegate func) const{
		ParamSet P;
		Full_Cycle(func,P);
	}

	
	AbstractPlotStream::AbstractPlotStream(string&& name, shared_ptr<PlotEngine> plot){
		m_plot=plot;
		m_name=name;
	}
	AbstractPlotStream::AbstractPlotStream(string&& name):
	AbstractPlotStream(static_cast<string&&>(name),make_shared<PlotEngine>()){}
	string&&AbstractPlotStream::Name() const{return const_cast<string&&>(m_name);}
	shared_ptr<PlotEngine> AbstractPlotStream::Plot(){return m_plot;}
	AbstractPlotStream::~AbstractPlotStream(){}
	AbstractPlotStream& AbstractPlotStream::operator<<(const ParamSet& P){
		ProcessPoint(P);
		return *this;
	}
	AbstractPlotStream& AbstractPlotStream::operator<<(ParamSet&& P){return operator<<(P);}

	SimplePlotStream::SimplePlotStream(string&&name, pair<size_t,size_t>&&indexes,shared_ptr<PlotEngine>plot)
		:AbstractPlotStream(static_cast<string&&>(name), plot){m_indexes=indexes;m_xrange=make_pair(+INFINITY,-INFINITY);}
	SimplePlotStream::SimplePlotStream(std::string&& name,pair<size_t,size_t>&& indexes)
		:AbstractPlotStream(static_cast<string&&>(name)){m_indexes=indexes;m_xrange=make_pair(+INFINITY,-INFINITY);}
	SimplePlotStream& SimplePlotStream::AddFunc(SimplePlotStream::func f){
		m_funcs.push_back(f);
		return *this;
	}
	SimplePlotStream::~SimplePlotStream(){
		if(m_data.size()>0)
			Plot()->WithoutErrors(Name(),m_data);
		size_t cnt=1;
		for(func F:m_funcs){
			Plot()->Line(Name()+" Line "+to_string(cnt),F,m_xrange.first,m_xrange.second,(m_xrange.second-m_xrange.first)/100.0);
			cnt++;
		}
	}
	void SimplePlotStream::ProcessPoint(const ParamSet& P){
		auto p=make_pair(P[m_indexes.first],P[m_indexes.second]);
		m_data.push_back(p);
		if(m_xrange.first>p.first)m_xrange.first=p.first;
		if(m_xrange.second<p.second)m_xrange.second=p.second;
	}
};
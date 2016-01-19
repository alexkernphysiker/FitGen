// this file is distributed under 
// MIT license
#include <Genetic/paramsort.h>
namespace Genetic{
	using namespace std;
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
			throw math_h_error<BinningParam>("wrong binning ranges");
		if(0==cnt)
			throw math_h_error<BinningParam>("there cannot be zero bins");
	}
	double BinningParam::bin_width() const{
		if(0==count())throw math_h_error<BinningParam>("count==0");
		return (range.second-range.first)/double(count());
	}
	double BinningParam::bin_center(size_t i) const{
		if(i>=count())throw math_h_error<BinningParam>("bin range check error");
		return range.first+bin_width()*(double(i)+0.5);
	}
	bool BinningParam::FindBinIndex(const ParamSet& P, size_t& res) const{
		double x=P[index],delta=bin_width()/2.0;
		if(delta<=0)throw math_h_error<BinningParam>("delta<=0");
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
	shared_ptr<vector<ParamSet>> ParamsPerBins::CreateParamProcessor(size_t){return make_shared<vector<ParamSet>>();}
	void ParamsPerBins::ProcessParams(vector<ParamSet>& proc, const ParamSet& P){proc.push_back(P);}
	
	ParamsPerBinsCounter<1>::ParamsPerBinsCounter(const vector<BinningParam>& binning):AbstractPerBinSeparator<unsigned long>(binning[binning.size()-1]){Init();}
	ParamsPerBinsCounter<1>::ParamsPerBinsCounter(vector< BinningParam >&& binning):ParamsPerBinsCounter(binning){}
	ParamsPerBinsCounter<1>::ParamsPerBinsCounter(const vector<BinningParam>*binning):ParamsPerBinsCounter(*binning){}
	shared_ptr<unsigned long> ParamsPerBinsCounter<1>::CreateParamProcessor(size_t){return make_shared<unsigned long>(0);}
	void ParamsPerBinsCounter<1>::ProcessParams(long unsigned int& proc, const ParamSet&){proc++;}

};
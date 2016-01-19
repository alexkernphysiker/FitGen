// this file is distributed under 
// MIT license
#include <math.h>
#include <gtest/gtest.h>
#include <math_h/exception_math_h.h>
#include <Genetic/paramsort.h>
using namespace Genetic;
using namespace std;
TEST(BinningParam,Throwing){
	EXPECT_NO_THROW(BinningParam(0,make_pair(0,1),1));
	EXPECT_THROW(BinningParam(0,make_pair(0,1),0),math_h_error<BinningParam>);
	EXPECT_THROW(BinningParam(0,make_pair(0,0),0),math_h_error<BinningParam>);
	EXPECT_THROW(BinningParam(0,make_pair(1,0),0),math_h_error<BinningParam>);
	EXPECT_THROW(BinningParam(0,make_pair(0,0),1),math_h_error<BinningParam>);
	EXPECT_THROW(BinningParam(0,make_pair(1,0),1),math_h_error<BinningParam>);
}
TEST(BinningParam,Base){
	BinChecker b(0,make_pair(0,1),10);
	ASSERT_EQ(0,b.param_index());
	ASSERT_EQ(10,b.count());
	EXPECT_EQ(0.1,b.bin_width());
	for(size_t i=0;i<10;i++)
		EXPECT_TRUE(pow(double(i)/10+0.05-b.bin_center(i),2)<0.00000001);
	size_t i=0;
	ParamSet P={0.72};
	EXPECT_EQ(true,b.FindBin(P,i));
	EXPECT_EQ(7,i);
	P[0]=0.49;
	EXPECT_EQ(true,b.FindBin(P,i));
	EXPECT_EQ(4,i);
	P[0]=0.33;
	EXPECT_EQ(true,b.FindBin(P,i));
	EXPECT_EQ(3,i);
	P[0]=1.01;
	EXPECT_EQ(false,b.FindBin(P,i));
	P[0]=-0.01;
	EXPECT_EQ(false,b.FindBin(P,i));
}
TEST(ParamsPerBins,Base){
	ParamsPerBins Binner(0,make_pair(0,1),4);
	Binner<<ParamSet({-1})<<ParamSet({0.15})<<ParamSet({0.85})<<ParamSet({2})<<ParamSet({0.45})<<ParamSet({0.55});
	ASSERT_EQ(4,Binner.count());
	EXPECT_EQ(1,Binner[0].size());
	EXPECT_EQ(1,Binner[1].size());
	EXPECT_EQ(1,Binner[2].size());
	EXPECT_EQ(1,Binner[3].size());
	EXPECT_EQ(1,Binner[0][0].size());
	EXPECT_EQ(1,Binner[1][0].size());
	EXPECT_EQ(1,Binner[2][0].size());
	EXPECT_EQ(1,Binner[3][0].size());
	EXPECT_EQ(0.15,Binner[0][0][0]);
	EXPECT_EQ(0.45,Binner[1][0][0]);
	EXPECT_EQ(0.55,Binner[2][0][0]);
	EXPECT_EQ(0.85,Binner[3][0][0]);
}
TEST(ParamsPerBinsCounter,Base){
	ParamsPerBinsCounter<1> Binner({BinningParam(0,make_pair(0,1),4)});
	Binner<<ParamSet({-1})<<ParamSet({0.15})<<ParamSet({0.85})<<ParamSet({2})<<ParamSet({0.45})<<ParamSet({0.55});
	ASSERT_EQ(4,Binner.count());
	EXPECT_EQ(1,Binner[0]);
	EXPECT_EQ(1,Binner[1]);
	EXPECT_EQ(1,Binner[2]);
	EXPECT_EQ(1,Binner[3]);
}
TEST(ParamsPerBinsCounter2,Base){
	ParamsPerBinsCounter<2> Binner({BinningParam(0,make_pair(0,1),2),BinningParam(1,make_pair(1,2),2)});
	Binner<<ParamSet({0.2,1.2})<<ParamSet({0.8,1.2})<<ParamSet({0.2,1.8})<<ParamSet({0.8,1.8})
		<<ParamSet({0.2,0})<<ParamSet({0.2,3})<<ParamSet({-1,1.2})<<ParamSet({2,1.2})
		<<ParamSet({-1,0})<<ParamSet({-1,3})<<ParamSet({3,0})<<ParamSet({3,3});
	ASSERT_EQ(2,Binner.count());
	ASSERT_EQ(2,Binner[0].count());
	ASSERT_EQ(2,Binner[1].count());
	EXPECT_EQ(1,Binner[0][0]);
	EXPECT_EQ(1,Binner[0][1]);
	EXPECT_EQ(1,Binner[1][0]);
	EXPECT_EQ(1,Binner[1][1]);
}

// this file is distributed under 
// GPL v 3.0 license
#include <gtest/gtest.h>
#include <exception_math_h.h>
#include <filter.h>
using namespace Genetic;
using namespace std;
TEST(Filter, BaseTest){
	Filter must_be_false([](const ParamSet&){return false;});
	Filter must_be_true([](const ParamSet&){return true;});
	ParamSet p;
	EXPECT_EQ(false,must_be_false(static_cast<const ParamSet&>(p)));
	EXPECT_EQ(true,must_be_true(static_cast<const ParamSet&>(p)));
}
TEST(Above,BasicTest){
	for(int count=0;count<10;count++){
		ParamSet Check=parZeros(10);
		{Above tester;
			for(int i=0;i<count;i++){
				EXPECT_EQ(&tester,&(tester<<1));
				EXPECT_EQ(false,tester(static_cast<const ParamSet&>(Check)));
			}
		}
		{Above tester;
			for(int i=0;i<count;i++){
				EXPECT_EQ(&tester,&(tester<<-1));
				EXPECT_EQ(true,tester(static_cast<const ParamSet&>(Check)));
			}
		}
		{Above tester;double c=1;
			for(int i=0;i<count;i++){
				EXPECT_EQ(&tester,&(tester<<c));
				EXPECT_EQ(false,tester(static_cast<const ParamSet&>(Check)));
				c=-c;
			}
		}
	}
}
TEST(Below,BasicTest){
	for(int count=0;count<10;count++){
		ParamSet Check=parZeros(10);
		{Below tester;
			for(int i=0;i<count;i++){
				EXPECT_EQ(&tester,&(tester<<-1));
				EXPECT_EQ(false,tester(static_cast<const ParamSet&>(Check)));
			}
		}
		{Below tester;
			for(int i=0;i<count;i++){
				EXPECT_EQ(&tester,&(tester<<1));
				EXPECT_EQ(true,tester(static_cast<const ParamSet&>(Check)));
			}
		}
		{Below tester;double c=-1;
			for(int i=0;i<count;i++){
				EXPECT_EQ(&tester,&(tester<<c));
				EXPECT_EQ(false,tester(static_cast<const ParamSet&>(Check)));
				c=-c;
			}
		}
	}
}
const int n=5;
shared_ptr<Filter> filters[]={
	make_shared<Filter>([](const ParamSet&){return true;}),
	make_shared<Filter>([](const ParamSet&){return true;}),
	make_shared<Filter>([](const ParamSet&){return true;}),
	make_shared<Filter>([](const ParamSet&){return true;}),
	make_shared<Filter>([](const ParamSet&){return true;})
};
template<class FilterMulti>void test_multi(){
	for(int count=0;count<n;count++){
		FilterMulti I;
		for(int i=0;i<count;i++)
			EXPECT_EQ(&I,&(I.Add(filters[i])));
		EXPECT_EQ(count,I.Count());
		for(int i=0;i<count;i++)
			EXPECT_EQ(filters[i].get(),&I.Get(i));
		EXPECT_THROW(I.Get(count),math_h_error<AbstractFilterMulti>);
		EXPECT_THROW(I.Get(-1),math_h_error<AbstractFilterMulti>);
	}
}
TEST(And,Add){test_multi<And>();}
TEST(Or,Add){test_multi<Or>();}
template<class FilterMulti>void test_add(){
	auto I=make_shared<FilterMulti>();
	EXPECT_EQ(0,I->Count());
	auto filter=[](const ParamSet&){return true;};
	EXPECT_EQ(I.get(),&(I->Add(filter)));
	EXPECT_EQ(1,I->Count());
	EXPECT_EQ(I.get(),(I<<filter).get());
	EXPECT_EQ(2,I->Count());
	EXPECT_EQ(I.get(),(I<<make_shared<Filter>(filter)).get());
	EXPECT_EQ(3,I->Count());
}
TEST(And,Add2){test_add<And>();}
TEST(Or,Add2){test_add<Or>();}
class ConstFilter:public IParamCheck{
private:
	bool value;
public:
	ConstFilter(bool v):value(v){}
	virtual ~ConstFilter(){}
	virtual bool operator()(const ParamSet&P)const override{return value;}
};
#define TRUE make_shared<ConstFilter>(true)
#define FALSE make_shared<ConstFilter>(false)
TEST(And,Work){
	ParamSet P;
	EXPECT_EQ(true,(make_shared<And>()<<TRUE<<TRUE)->operator()(static_cast<const ParamSet&>(P)));
	EXPECT_EQ(false,(make_shared<And>()<<FALSE<<TRUE)->operator()(static_cast<const ParamSet&>(P)));
	EXPECT_EQ(false,(make_shared<And>()<<TRUE<<FALSE)->operator()(static_cast<const ParamSet&>(P)));
	EXPECT_EQ(false,(make_shared<And>()<<FALSE<<FALSE)->operator()(static_cast<const ParamSet&>(P)));
}
TEST(Or,Work){
	ParamSet P;
	EXPECT_EQ(true,(make_shared<Or>()<<TRUE<<TRUE)->operator()(static_cast<const ParamSet&>(P)));
	EXPECT_EQ(true,(make_shared<Or>()<<FALSE<<TRUE)->operator()(static_cast<const ParamSet&>(P)));
	EXPECT_EQ(true,(make_shared<Or>()<<TRUE<<FALSE)->operator()(static_cast<const ParamSet&>(P)));
	EXPECT_EQ(false,(make_shared<Or>()<<FALSE<<FALSE)->operator()(static_cast<const ParamSet&>(P)));
}

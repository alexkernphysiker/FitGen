#include <gtest/gtest.h>
#include <filter.h>
#include <genetic_exception.h>
using namespace Genetic;
using namespace std;

TEST(Above,BasicTest){
	for(int count=0;count<10;count++){
		ParamSet Check=parZeros(10);
		{Above tester;
			for(int i=0;i<count;i++){
				EXPECT_EQ(&tester,&(tester<<1));
				EXPECT_EQ(false,tester(Check));
			}
		}
		{Above tester;
			for(int i=0;i<count;i++){
				EXPECT_EQ(&tester,&(tester<<-1));
				EXPECT_EQ(true,tester(Check));
			}
		}
		{Above tester;double c=1;
			for(int i=0;i<count;i++){
				EXPECT_EQ(&tester,&(tester<<c));
				EXPECT_EQ(false,tester(Check));
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
				EXPECT_EQ(false,tester(Check));
			}
		}
		{Below tester;
			for(int i=0;i<count;i++){
				EXPECT_EQ(&tester,&(tester<<1));
				EXPECT_EQ(true,tester(Check));
			}
		}
		{Below tester;double c=-1;
			for(int i=0;i<count;i++){
				EXPECT_EQ(&tester,&(tester<<c));
				EXPECT_EQ(false,tester(Check));
				c=-c;
			}
		}
	}
}
const int n=5;
shared_ptr<Filter> filters[]={
	make_shared<Filter>([](ParamSet&){return true;}),
	make_shared<Filter>([](ParamSet&){return true;}),
	make_shared<Filter>([](ParamSet&){return true;}),
	make_shared<Filter>([](ParamSet&){return true;}),
	make_shared<Filter>([](ParamSet&){return true;})
};
template<class FilterMulti>void test_multi(){
	for(int count=0;count<n;count++){
		FilterMulti I;
		for(int i=0;i<count;i++)
			EXPECT_EQ(&I,&(I.Add(filters[i])));
		EXPECT_EQ(count,I.Count());
		for(int i=0;i<count;i++)
			EXPECT_EQ(filters[i].get(),&I.Get(i));
		ASSERT_ANY_THROW(I.Get(count));
		ASSERT_ANY_THROW(I.Get(-1));
	}
}
TEST(And,Add){test_multi<And>();}
TEST(Or,Add){test_multi<Or>();}
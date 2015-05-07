#include <gtest/gtest.h>
#include <paramset.h>
#include <genetic_exception.h>
using namespace Genetic;
using namespace std;
TEST(ParamSet, Constructor0){
	ParamSet p;
	EXPECT_EQ(0,p.Count());
}
TEST(ParamSet, Constructor1){
	ParamSet p(1);
	EXPECT_EQ(1,p.Count());
	for(int i=1;i<=1;i++)
		EXPECT_EQ(i,p[i-1]);
}
TEST(ParamSet, Constructor2){
	ParamSet p(1,2);
	EXPECT_EQ(2,p.Count());
	for(int i=1;i<=2;i++)
		EXPECT_EQ(i,p[i-1]);
}
TEST(ParamSet, Constructor3){
	ParamSet p(1,2,3);
	EXPECT_EQ(3,p.Count());
	for(int i=1;i<=3;i++)
		EXPECT_EQ(i,p[i-1]);
}
TEST(ParamSet, Constructor4){
	ParamSet p(1,2,3,4);
	EXPECT_EQ(4,p.Count());
	for(int i=1;i<=4;i++)
		EXPECT_EQ(i,p[i-1]);
}
TEST(ParamSet, Constructor5){
	ParamSet p(1,2,3,4,5);
	EXPECT_EQ(5,p.Count());
	for(int i=1;i<=5;i++)
		EXPECT_EQ(i,p[i-1]);
}
TEST(ParamSet, Constructor6){
	ParamSet p(1,2,3,4,5,6);
	EXPECT_EQ(6,p.Count());
	for(int i=1;i<=6;i++)
		EXPECT_EQ(i,p[i-1]);
}
const int n=10;
double data[]={1,2,3,4,5,6,7,8,9,10};
void TestEq(ParamSet &p1,ParamSet &p2){
	EXPECT_EQ(p1.Count(),p2.Count());
	for(int i=0,n=p1.Count();i<n;i++)
		EXPECT_EQ(p1[i],p2[i]);
}
TEST(ParamSet, AddNumber){
	ParamSet p;
	for(int i=0;i<n;i++){
		EXPECT_EQ(&p,&(p<<data[i]));
		EXPECT_EQ(i+1,p.Count());
		EXPECT_EQ(data[i],p[p.Count()-1]);
	}
}
TEST(ParamSet, AddParamSet){
	for(int count=0;count<=n;count++){
		ParamSet P;
		for(int i=0;i<5;i++){
			ParamSet p;
			for(int i=0;i<count;i++)p<<data[i];
			int c=P.Count();
			EXPECT_EQ(&P,&(P<<p));
			EXPECT_EQ(c+count,P.Count());
			for(int i=0;i<count;i++)
				EXPECT_EQ(P[P.Count()-1-i],p[p.Count()-1-i]);
		}
	}
}
TEST(ParamSet, AddToItself){
	for(int count=0;count<=n;count++){
		ParamSet P;
		for(int i=0;i<count;i++)P<<data[i];
		int c=P.Count();
		EXPECT_EQ(&P,&(P<<P));
		EXPECT_EQ(c,P.Count());
	}
}
TEST(ParamSet, Copying){
	for(int count=0;count<=n;count++){
		ParamSet source;
		for(int i=0;i<count;i++)source<<data[i];
		ParamSet copy(source);
		TestEq(source,copy);
		ParamSet copy2;
		copy2=source;
		TestEq(source,copy2);
		copy=source;
		TestEq(source,copy);
	}
}
TEST(ParamSet, Set){
	for(int count=0;count<=n;count++){
		for(int si=0;si<count;si++){
			ParamSet source;
			for(int i=0;i<count;i++)source<<data[i];
			source.Set(si,0.0);
			EXPECT_EQ(count,source.Count());
			for(int i=0;i<count;i++){
				if(i==si)
					EXPECT_EQ(0.0,source[i]);
				else
					EXPECT_EQ(data[i],source[i]);
			}
		}
		ParamSet source;
		for(int i=0;i<count;i++)source<<data[i];
		ASSERT_ANY_THROW(source.Set(-1,0.0));
		ASSERT_ANY_THROW(source.Set(source.Count(),0.0));
	}
}
TEST(ParamSet,ParEQ){
	double val=1.2;
	for(int count=0;count<=n;count++){
		ParamSet P=parEq(count,val);
		EXPECT_EQ(count,P.Count());
		for(int i=0;i<count;i++)
			EXPECT_EQ(val,P[i]);
	}
}

#include <gtest/gtest.h>
#include <initialconditions.h>
#include <genetic_exception.h>
using namespace Genetic;
using namespace std;
TEST(Initialiser,Create){
	Initialiser I;
	EXPECT_EQ(0,I.Count());
}
const int n=5;
generator Funcs[]={[](){return 0;},[](){return 1;},[](){return 2;},[](){return 3;},[](){return 4;}};
TEST(Initialiser,Add){
	for(int count=1;count<n;count++){
		Initialiser I;
		for(int i=0;i<count;i++)I<<Funcs[i];
		EXPECT_EQ(count,I.Count());
		for(int i=0;i<count;i++)
			EXPECT_EQ(Funcs[i](),I[i]());
		//EXPECT_ANY_THROW(I[count]());
		//EXPECT_ANY_THROW(I[-1]());
	}
}
TEST(Initialiser,Generate){
	for(int count=1;count<5;count++){
		Initialiser I;
		for(int i=0;i<count;i++)I<<Funcs[i];
		ParamSet P=I.Generate();
		EXPECT_EQ(count,P.Count());
		for(int i=0;i<count;i++)EXPECT_EQ(Funcs[i](),P[i]);
	}
}


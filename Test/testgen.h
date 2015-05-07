#include <abstract.h>
class GeneticTest:public Genetic::AbstractGenetic{
public:
	GeneticTest(std::shared_ptr<Genetic::IOptimalityFunction> optimality):AbstractGenetic(optimality){}
	virtual ~GeneticTest(){}
};

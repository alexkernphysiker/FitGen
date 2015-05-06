#ifndef ____EuiGhCCk
#define ____EuiGhCCk
#include <exception>
#include <string>
#include <functional>
void SetGeneticErrorOutput(std::function<void(std::string)> func);
class GeneticException:public std::exception{
public:
	GeneticException(std::string msg);
	virtual ~GeneticException() throw(){}
	virtual const char* what() const throw();
private:
	std::string m_msg;
};
#endif
